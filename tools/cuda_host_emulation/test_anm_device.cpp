// Host emulation of the CUDA ANM device functions, compared with the CPU ANMModel
#define HOST_EMULATION
#include "cuda_runtime_api.h"
#include "Interactions/DRHInteraction.h"
#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "CUDA/Interactions/CUDA_DRH.cuh"
#include "CUDA/Interactions/CUDA_ANM.cuh"
#include "Interactions/ANMModel.h"
#include "Particles/ANMParticle.h"
#include "Particles/DNANucleotide.h"
#include "model.h"
#include "Utilities/Logger.h"
#include <cstdio>
#include <fstream>
#include <random>

static std::mt19937 rng(42);
static double uni(double a, double b) { return std::uniform_real_distribution<double>(a, b)(rng); }

static LR_matrix random_rotation() {
	LR_vector axis(uni(-1, 1), uni(-1, 1), uni(-1, 1));
	axis.normalize();
	double t = uni(0, 2 * M_PI);
	double c = cos(t), s = sin(t);
	return LR_matrix(c + axis.x * axis.x * (1 - c), axis.x * axis.y * (1 - c) - axis.z * s, axis.x * axis.z * (1 - c) + axis.y * s,
			axis.y * axis.x * (1 - c) + axis.z * s, c + axis.y * axis.y * (1 - c), axis.y * axis.z * (1 - c) - axis.x * s,
			axis.z * axis.x * (1 - c) - axis.y * s, axis.z * axis.y * (1 - c) + axis.x * s, c + axis.z * axis.z * (1 - c));
}

static void set_orientation(BaseParticle *p, const LR_matrix &R) {
	p->orientation = R;
	p->orientationT = R.get_transpose();
	p->set_positions();
	p->force = LR_vector(0, 0, 0);
	p->torque = LR_vector(0, 0, 0);
}

static c_number4 c4(const LR_vector &v) { return make_c_number4(v.x, v.y, v.z, 0); }
static double diff(const c_number4 &a, const LR_vector &b) { return std::max({fabs(a.x - b.x), fabs(a.y - b.y), fabs(a.z - b.z)}); }

static cuda_anm::RepulsiveLJ dev(const ANMModel::RepulsiveLJ &lj) {
	return { lj.sigma, lj.rstar, lj.b, lj.rcut, lj.sqr_rcut, lj.stiffness };
}

int main() {
	Logger::init();
	double worst = 0;
	auto report = [&worst](const char *what, double d) {
		if(d > worst) worst = d;
		if(d > 1e-9) printf("MISMATCH %s: %g\n", what, d);
	};

	// --- protein-protein excluded volume
	ANMModel anm(true);
	cuda_anm::anm_params.protein = dev(anm.protein_lj());
	cuda_anm::anm_params.backbone = dev(anm.backbone_lj());
	cuda_anm::anm_params.base = dev(anm.base_lj());
	cuda_anm::anm_params.torsion = true;
	ANMParticle p(true), q(true);
	p.index = 0;
	q.index = 1;
	for(int t = 0; t < 2000; t++) {
		set_orientation(&p, random_rotation());
		set_orientation(&q, random_rotation());
		LR_vector r(uni(-1, 1), uni(-1, 1), uni(-1, 1));
		r.normalize();
		r *= uni(0.3, 0.36);
		number E = anm.excluded_volume(&p, &q, r, true);
		c_number4 F = cuda_anm::_anm_repulsive_lj(c4(r), cuda_anm::anm_params.protein);
		if(t < 3) printf("LJ E cpu %.15g gpu %.15g\n", E, F.w);
		report("protein LJ energy", 1e-3 * fabs(F.w - E) / std::max(1.0, fabs(E)));
		report("protein LJ force", diff(F, q.force) / std::max(1.0, q.force.module()));
	}

	// --- bending and torsion
	{
		std::ofstream par("/tmp/anm_devtest.par");
		par << "2\n0 1 0.38 s 1.0 0.3 -0.2 0.7 0.5\n";
	}
	input_file inp;
	inp.init_from_string("parfile = /tmp/anm_devtest.par\nbending_k = 37\ntorsion_k = 23\n");
	anm.get_settings(inp, true);
	anm.parse_parameter_file(2);
	const auto &ang = anm.angular_parameters().at(0);
	c_number cosines[4] = { ang.cosines[0], ang.cosines[1], ang.cosines[2], ang.cosines[3] };
	for(int t = 0; t < 2000; t++) {
		set_orientation(&p, random_rotation());
		set_orientation(&q, random_rotation());
		LR_vector r(uni(-1, 1), uni(-1, 1), uni(-1, 1));
		r.normalize();
		r *= uni(0.3, 0.5);
		number E = anm.angular(&p, &q, r, true);
		c_number4 a1 = c4(p.orientationT.v1), a2 = c4(p.orientationT.v2), a3 = c4(p.orientationT.v3);
		c_number4 b1 = c4(q.orientationT.v1), b2 = c4(q.orientationT.v2), b3 = c4(q.orientationT.v3);
		c_number4 Fq, Tp, Tq;
		cuda_anm::_anm_angular(c4(r), a1, a3, b1, b3, cosines, 37, 23, Fq, Tp, Tq);
		if(t < 3) printf("ang E cpu %.12g gpu %.12g | Fq cpu %.9g %.9g %.9g gpu %.9g %.9g %.9g\n", E, Fq.w, q.force.x, q.force.y, q.force.z, Fq.x, Fq.y, Fq.z);
		report("angular energy", 1e-3 * fabs(Fq.w - E) / std::max(1.0, fabs(E)));
		report("angular force q", diff(Fq, q.force));
		report("angular force p", diff(-Fq, p.force));
		// the CPU torques are in the frame of each particle
		report("angular torque p", diff(_vectors_transpose_c_number4_product(a1, a2, a3, Tp), p.torque));
		report("angular torque q", diff(_vectors_transpose_c_number4_product(b1, b2, b3, Tq), q.torque));
	}

	// --- amino acid - DNA nucleotide excluded volume
	cuda_drh::HybridModel hm;
	memset(&hm, 0, sizeof(hm));
	hm.dna_back1 = POS_MM_BACK1;
	hm.dna_back2 = POS_MM_BACK2;
	hm.dna_base = POS_BASE;
	cuda_drh::hybrid_model = hm;
	DNANucleotide nuc(true);
	int contacts = 0;
	for(int t = 0; t < 20000; t++) {
		set_orientation(&p, random_rotation());
		set_orientation(&nuc, random_rotation());
		LR_vector r(uni(-1, 1), uni(-1, 1), uni(-1, 1));
		r.normalize();
		r *= uni(0.0, 1.2);
		number E = anm.nucleotide_excluded_volume(&p, &nuc, r, true);
		c_number4 n1 = c4(nuc.orientationT.v1), n2 = c4(nuc.orientationT.v2), n3 = c4(nuc.orientationT.v3);
		c_number4 a1 = c4(p.orientationT.v1), a2 = c4(p.orientationT.v2), a3 = c4(p.orientationT.v3);
		c_number4 Fp, Tn, Tp;
		cuda_anm::_anm_nucleotide_excluded_volume(c4(r), n1, n2, n3, false, Fp, Tn, Tp);
		if(E != 0) contacts++;
		report("nucleotide exc energy", 1e-3 * fabs(Fp.w - E) / std::max(1.0, fabs(E)));
		report("nucleotide exc force protein", diff(Fp, p.force) / std::max(1.0, p.force.module()));
		report("nucleotide exc force nucleotide", diff(-Fp, nuc.force) / std::max(1.0, nuc.force.module()));
		report("nucleotide exc torque nucleotide", diff(_vectors_transpose_c_number4_product(n1, n2, n3, Tn), nuc.torque) / std::max(1.0, nuc.torque.module()));
		// the force is central, so this torque is round-off noise proportional to the force
		report("nucleotide exc torque protein", diff(_vectors_transpose_c_number4_product(a1, a2, a3, Tp), p.torque) / std::max(1.0, p.force.module()));
	}
	printf("nucleotide contacts tested: %d\n", contacts);
	printf("largest relative difference: %g\n", worst);
	return worst > 1e-9;
}
