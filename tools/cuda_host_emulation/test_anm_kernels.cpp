// Runs the CUDA ANM kernels (edge-based and neighbour-matrix) on the host and compares the force and torque on every particle,
// and the total energy (sum of the w components / 2, as MD_CUDABackend computes it), with the CPU interaction.
// usage: test_anm_kernels <interaction> <topology> <conf> <parfile|none> [kb kt]
#define HOST_EMULATION
#include "cuda_runtime_api.h"
#include "Interactions/DNANMInteraction.h"
#include "Interactions/RNANMInteraction.h"
#include "Interactions/NANMInteraction.h"
#include "Particles/ANMParticle.h"
#include "Particles/RNANucleotide.h"
#include "Boxes/CubicBox.h"
#include "Utilities/Logger.h"
#include "Utilities/ConfigInfo.h"
#include "Particles/Molecule.h"
#include "CUDA/CUDAUtils.h"
#include "CUDA/cuda_utils/CUDABox.h"
#include "CUDA/Lists/CUDABaseList.h"
#include "CUDA/Interactions/CUDA_DRH.cuh"
#include "CUDA/Interactions/CUDA_NA_constants.cuh"
#include "CUDA/Interactions/CUDA_ANM.cuh"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <type_traits>

#include "copy_model.inc"

template<class Base>
class Testable: public Base {
	friend struct NAConstants;
public:
	explicit Testable(bool torsion) : Base(torsion) {}
	const ANMModel &model() { return this->_anm; }
};

template<class I> using has_dna = std::is_base_of<DNA2Interaction, I>;
template<class I> using has_rna = std::is_base_of<RNA2Interaction, I>;
template<class I> using has_hybrid = std::is_base_of<DRHInteraction, I>;
template<class I> void init_dna(I &inter, int N, cuda_drh::HybridModel &hm, cuda_drh::KernelParameters &params, std::true_type) {
	NAConstants::dna(inter, N);
	NAConstants::dna_sites(hm);
	NAConstants::dna_kernel_parameters(inter, params);
}
template<class I> void init_dna(I &, int, cuda_drh::HybridModel &, cuda_drh::KernelParameters &, std::false_type) {}
template<class I> void init_rna(I &inter, int N, cuda_drh::HybridModel &hm, cuda_drh::KernelParameters &params, std::true_type) {
	NAConstants::rna(inter, N);
	NAConstants::rna_sites(hm, inter.RNAInteraction::get_model());
	NAConstants::rna_kernel_parameters(inter, params);
}
template<class I> void init_rna(I &, int, cuda_drh::HybridModel &, cuda_drh::KernelParameters &, std::false_type) {}
template<class I> void init_hybrid(I &inter, int N, const cuda_drh::HybridModel &, std::true_type) { NAConstants::hybrid(inter, N, true, 1); }
template<class I> void init_hybrid(I &, int N, const cuda_drh::HybridModel &hm, std::false_type) { NAConstants::copy_hybrid_model(hm, N, true, 1); }

static cuda_anm::RepulsiveLJ dev(const ANMModel::RepulsiveLJ &lj) {
	return { lj.sigma, lj.rstar, lj.b, lj.rcut, lj.sqr_rcut, lj.stiffness };
}

static void launch(int n, const std::function<void()> &kernel) {
	blockDim = dim3(1, 1, 1);
	gridDim = dim3(n, 1, 1);
	threadIdx.x = threadIdx.y = threadIdx.z = 0;
	for(int i = 0; i < n; i++) {
		blockIdx.x = i;
		blockIdx.y = blockIdx.z = 0;
		kernel();
	}
}

template<class Base>
int run(const std::string &interaction, bool torsion, const std::string &top, const std::string &conf, const std::string &par, const std::string &extra) {
	static std::vector<BaseParticle *> config_particles;
	static std::vector<std::shared_ptr<Molecule>> config_mols;
	ConfigInfo::init(&config_particles, &config_mols);
	Testable<Base> inter(torsion);
	input_file inp;
	std::string s = "topology = " + top + "\nparfile = " + par + "\nT = 300K\nsalt_concentration = 1.0\n" + extra;
	inp.init_from_string(s);
	inter.get_settings(inp);
	inter.init();

	std::ifstream tf(top);
	int N;
	tf >> N;
	std::vector<BaseParticle *> particles(N);
	int N_strands;
	inter.read_topology(&N_strands, particles);
	config_particles = particles;

	// configuration
	std::ifstream cf(conf);
	std::string line;
	std::getline(cf, line);
	std::getline(cf, line);
	std::stringstream bs(line.substr(line.find('=') + 1));
	double L;
	bs >> L;
	std::getline(cf, line);
	CubicBox box;
	box.init(L, L, L);
	inter.set_box(&box);
	for(int i = 0; i < N; i++) {
		std::getline(cf, line);
		std::stringstream ls(line);
		LR_vector pos, a1, a3;
		ls >> pos.x >> pos.y >> pos.z >> a1.x >> a1.y >> a1.z >> a3.x >> a3.y >> a3.z;
		a3.normalize();
		a1 -= a3 * (a1 * a3);
		a1.normalize();
		LR_vector a2 = a3.cross(a1);
		BaseParticle *p = particles[i];
		p->pos = pos;
		p->orientationT = LR_matrix(a1, a2, a3);
		p->orientation = p->orientationT.get_transpose();
		p->set_positions();
		p->force = p->torque = LR_vector(0, 0, 0);
	}

	// CPU forces and torques, all pairs
	double cpu_energy = 0.;
	for(int i = 0; i < N; i++) {
		for(int j = i + 1; j < N; j++) {
			cpu_energy += inter.pair_interaction(particles[i], particles[j], true, true);
		}
	}

	// GPU data
	cuda_drh::HybridModel hm;
	memset(&hm, 0, sizeof(hm));
	cuda_drh::KernelParameters params;
	memset(&params, 0, sizeof(params));
	params.average = true;
	init_dna(inter, N, hm, params, has_dna<Base>());
	init_rna(inter, N, hm, params, has_rna<Base>());
	init_hybrid(inter, N, hm, has_hybrid<Base>());

	const ANMModel &anm = inter.model();
	cuda_anm::ANMParameters ap;
	ap.protein = dev(anm.protein_lj());
	ap.backbone = dev(anm.backbone_lj());
	ap.base = dev(anm.base_lj());
	ap.protein_nucleotide_sqr_rcut = anm.protein_nucleotide_sqr_rcut();
	ap.torsion = anm.torsion();
	ap.per_pair_stiffness = anm.per_pair_stiffness();
	ap.k_bend = anm.global_k_bend();
	ap.k_tor = anm.global_k_tor();
	cuda_anm::anm_params = ap;
	cuda_anm::MD_N[0] = N;
	cuda_anm::MD_n_forces[0] = 1;

	std::vector<int> species(N);
	std::vector<c_number4> poss(N);
	std::vector<GPU_quat> orientations(N);
	std::vector<LR_bonds> bonds(N);
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		species[i] = dynamic_cast<ANMParticle *>(p) ? 2 : (dynamic_cast<RNANucleotide *>(p) ? 1 : 0);
		poss[i] = make_c_number4(p->pos.x, p->pos.y, p->pos.z, 0);
		const unsigned int mask22 = 0x003FFFFF;
		poss[i].w = GpuUtils::int_as_float((p->btype << 22) | (mask22 & p->index));
		bonds[i].n3 = (p->n3 == P_VIRTUAL) ? P_INVALID : p->n3->index;
		bonds[i].n5 = (p->n5 == P_VIRTUAL) ? P_INVALID : p->n5->index;
		LR_matrix &o = p->orientation;
		GPU_quat &q = orientations[i];
		number trace = o.v1.x + o.v2.y + o.v3.z;
		if(trace > 0) {
			number s = 0.5 / sqrt(trace + 1.0);
			q.w = 0.25 / s; q.x = (o.v3.y - o.v2.z) * s; q.y = (o.v1.z - o.v3.x) * s; q.z = (o.v2.x - o.v1.y) * s;
		}
		else if((o.v1.x > o.v2.y) && (o.v1.x > o.v3.z)) {
			number s = 0.5 / sqrt(1.0 + o.v1.x - o.v2.y - o.v3.z);
			q.w = (o.v3.y - o.v2.z) * s; q.x = 0.25 / s; q.y = (o.v1.y + o.v2.x) * s; q.z = (o.v1.z + o.v3.x) * s;
		}
		else if(o.v2.y > o.v3.z) {
			number s = 0.5 / sqrt(1.0 + o.v2.y - o.v1.x - o.v3.z);
			q.w = (o.v1.z - o.v3.x) * s; q.x = (o.v1.y + o.v2.x) * s; q.y = 0.25 / s; q.z = (o.v2.z + o.v3.y) * s;
		}
		else {
			number s = 0.5 / sqrt(1.0 + o.v3.z - o.v1.x - o.v2.y);
			q.w = (o.v2.x - o.v1.y) * s; q.x = (o.v1.z + o.v3.x) * s; q.y = (o.v2.z + o.v3.y) * s; q.z = 0.25 / s;
		}
	}

	// springs and angular parameters, as in CUDAANMInteraction::_init_arrays
	std::vector<std::vector<std::pair<int, const ANMModel::Spring *>>> per_particle(N);
	for(auto &spring : anm.springs()) {
		per_particle[spring.first.first].emplace_back(spring.first.second, &spring.second);
		per_particle[spring.first.second].emplace_back(spring.first.first, &spring.second);
	}
	std::vector<int> start(N + 1), neigh;
	std::vector<c_number> k, r0;
	for(int i = 0; i < N; i++) {
		std::sort(per_particle[i].begin(), per_particle[i].end(), [](const std::pair<int, const ANMModel::Spring *> &a, const std::pair<int, const ANMModel::Spring *> &b) { return a.first < b.first; });
		start[i] = neigh.size();
		for(auto &n : per_particle[i]) {
			neigh.push_back(n.first);
			k.push_back(n.second->k);
			r0.push_back(n.second->r0);
		}
	}
	start[N] = neigh.size();
	neigh.push_back(0); k.push_back(0); r0.push_back(0);
	std::vector<c_number> cosines(4 * N, 0), angk(2 * N, 0);
	for(auto &a : anm.angular_parameters()) {
		for(int c = 0; c < 4; c++) cosines[4 * a.first + c] = a.second.cosines[c];
		angk[2 * a.first] = a.second.k_bend;
		angk[2 * a.first + 1] = a.second.k_tor;
	}

	CUDABox cbox;
	cbox.change_sides(L, L, L);

	// edge list: all the pairs within the cutoff that are not n3/n5 neighbours
	number rcut = inter.get_rcut();
	std::vector<edge_bond> edges;
	for(int i = 0; i < N; i++) {
		for(int j = i + 1; j < N; j++) {
			if(bonds[i].n3 == j || bonds[i].n5 == j) continue;
			LR_vector r = box.min_image(particles[i]->pos, particles[j]->pos);
			if(r.norm() < SQR(rcut)) {
				edge_bond e;
				e.from = i;
				e.to = j;
				edges.push_back(e);
			}
		}
	}

	auto compare = [&](const char *name, const std::vector<c_number4> &F, const std::vector<c_number4> &T) {
		double fmax[3] = {0, 0, 0}, tmax[3] = {0, 0, 0}, fscale[3] = {0, 0, 0}, tscale[3] = {0, 0, 0};
		int worst[3] = {-1, -1, -1};
		for(int i = 0; i < N; i++) {
			BaseParticle *p = particles[i];
			int sp = species[i];
			double df = std::max({fabs(F[i].x - p->force.x), fabs(F[i].y - p->force.y), fabs(F[i].z - p->force.z)});
			double dt = std::max({fabs(T[i].x - p->torque.x), fabs(T[i].y - p->torque.y), fabs(T[i].z - p->torque.z)});
			if(df > fmax[sp]) { fmax[sp] = df; worst[sp] = i; }
			tmax[sp] = std::max(tmax[sp], dt);
			fscale[sp] = std::max(fscale[sp], (double) p->force.module());
			tscale[sp] = std::max(tscale[sp], (double) p->torque.module());
		}
		// as MD_CUDABackend (CUDA_print_energy): every pair is counted twice in the w components
		double gpu_energy = 0.;
		for(int i = 0; i < N; i++) gpu_energy += F[i].w;
		gpu_energy /= 2.;
		printf("  %-8s energy       GPU %.10e  CPU %.10e  diff %.3e (per particle %.3e)\n", name, gpu_energy, cpu_energy, gpu_energy - cpu_energy, (gpu_energy - cpu_energy) / N);
		const char *names[3] = {"DNA", "RNA", "amino acids"};
		for(int sp = 0; sp < 3; sp++) {
			if(fscale[sp] == 0 && tscale[sp] == 0) continue;
			printf("  %-8s %-12s max |dF| %.3e (max |F| %.3e)  max |dT| %.3e (max |T| %.3e)  worst particle %d\n", name, names[sp], fmax[sp], fscale[sp], tmax[sp], tscale[sp], worst[sp]);
		}
	};

	printf("%s, %d particles, %d springs, %d edges\n", interaction.c_str(), N, (int) anm.springs().size(), (int) edges.size());

	// edge-based kernels
	std::vector<c_number4> F(N, make_c_number4(0, 0, 0, 0)), T(N, make_c_number4(0, 0, 0, 0));
	launch(edges.size(), [&]() {
		cuda_anm::anm_forces_edge_nonbonded(poss.data(), orientations.data(), F.data(), T.data(), edges.data(), edges.size(), bonds.data(), species.data(), params, &cbox);
	});
	launch(N, [&]() {
		cuda_anm::anm_forces_edge_bonded(poss.data(), orientations.data(), F.data(), T.data(), bonds.data(), species.data(), params, start.data(), neigh.data(), k.data(), r0.data(), cosines.data(), angk.data(), &cbox);
	});
	compare("edge", F, T);

	// neighbour-matrix kernel with all the particles as neighbours
	std::vector<c_number4> F2(N, make_c_number4(0, 0, 0, 0)), T2(N, make_c_number4(0, 0, 0, 0));
	launch(N, [&]() {
		cuda_anm::anm_forces(poss.data(), orientations.data(), F2.data(), T2.data(), nullptr, nullptr, bonds.data(), species.data(), params, start.data(), neigh.data(), k.data(), r0.data(), cosines.data(), angk.data(), &cbox);
	});
	compare("matrix", F2, T2);
	return 0;
}

int main(int argc, char *argv[]) {
	Logger::init();
	Logger::instance()->disable_log();
	std::string inter = argv[1];
	std::string extra;
	if(argc > 5) extra = std::string("bending_k = ") + argv[5] + "\ntorsion_k = " + argv[6] + "\n";
	if(inter == "DNANM") return run<DNANMInteraction>(inter, false, argv[2], argv[3], argv[4], extra);
	if(inter == "DNACT") return run<DNANMInteraction>(inter, true, argv[2], argv[3], argv[4], extra);
	if(inter == "RNANM") return run<RNANMInteraction>(inter, false, argv[2], argv[3], argv[4], extra);
	if(inter == "RNACT") return run<RNANMInteraction>(inter, true, argv[2], argv[3], argv[4], extra);
	if(inter == "NANM") return run<NANMInteraction>(inter, false, argv[2], argv[3], argv[4], extra);
	if(inter == "NACT") return run<NANMInteraction>(inter, true, argv[2], argv[3], argv[4], extra);
	return 1;
}
