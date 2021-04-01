#pragma once
namespace STD {
    class CalcDensity {
        kernel_t kernel;
    public:
        void operator()(const EPI::Dens *const ep_i, const PS::S32 Nip, const EPJ::Dens *const ep_j, const PS::S32 Njp,
                        RESULT::Dens *const dens) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                const EPI::Dens &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const EPJ::Dens &jth = ep_j[j];
                    const PS::F64vec dr = jth.pos - ith.pos;
                    dens[i].dens += jth.mass * kernel.W(dr, ith.smth);
                }
#ifdef FLAG_GI
                dens[i].dens = std::max(100.0, dens[i].dens);
#endif
                dens[i].smth = PARAM::SMTH * pow(ith.mass / dens[i].dens, 1.0 / (PS::F64) (PARAM::Dim));
            }
        }
    };

    void CalcPressure(PS::ParticleSystem<STD::RealPtcl> &sph_system) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            sph_system[i].pres = sph_system[i].EoS->Pressure(sph_system[i].dens, sph_system[i].eng);
            sph_system[i].snds = sph_system[i].EoS->SoundSpeed(sph_system[i].dens, sph_system[i].eng);
        }
    }

    class CalcDerivative {
        kernel_t kernel;
    public:
        void operator()(const EPI::Drvt *ep_i, const PS::S32 Nip, const EPJ::Drvt *ep_j, const PS::S32 Njp,
                        RESULT::Drvt *const drvt) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                const EPI::Drvt &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const EPJ::Drvt &jth = ep_j[j];
                    const PS::F64vec dr = ith.pos - jth.pos;
                    const PS::F64vec dv = ith.vel - jth.vel;
                    drvt[i].div_v += -jth.mass * dv * kernel.gradW(dr, ith.smth);
                    drvt[i].rot_v += -jth.mass * dv ^ kernel.gradW(dr, ith.smth);
                    drvt[i].grad_smth -= jth.mass / ith.smth *
                                         (PARAM::Dim * kernel.W(dr, ith.smth) + dr * kernel.gradW(dr, ith.smth));
                }
                drvt[i].grad_smth = 1.0 / (1.0 + ith.smth * drvt[i].grad_smth / (PARAM::Dim * ith.dens));
                drvt[i].div_v /= ith.dens;
                drvt[i].rot_v /= ith.dens;
            }
        }
    };

    class CalcHydroForce {
        const kernel_t kernel;
    public:
        void
        operator()(const EPI::Hydro *const ep_i, const PS::S32 Nip, const EPJ::Hydro *const ep_j, const PS::S32 Njp,
                   RESULT::Hydro *const hydro) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                PS::F64 v_sig_max = 0.0;
                const EPI::Hydro &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const EPJ::Hydro &jth = ep_j[j];
                    const PS::F64vec dr = ith.pos - jth.pos;
                    const PS::F64vec dv = ith.vel - jth.vel;
                    const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
                    const PS::F64 v_sig = ith.snds + jth.snds - 3.0 * w_ij;
                    v_sig_max = std::max(v_sig_max, v_sig);
                    PS::F64 AV = -PARAM::AV_STRENGTH * 0.5 * v_sig * w_ij / (0.5 * (ith.dens + jth.dens)) * 0.5 *
                                 (ith.Bal + jth.Bal);
                    if (PARAM::FLAG_R00 == true) {
                        AV *= 0.5 * (ith.AVa + jth.AVa);
                    }
#if 1
                    const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ith.smth) * ith.grad_smth +
                                                    kernel.gradW(dr, jth.smth) * jth.grad_smth);
                    hydro[i].acc -= jth.mass *
                                    (ith.grad_smth * ith.pres / (ith.dens * ith.dens) * kernel.gradW(dr, ith.smth) +
                                     jth.grad_smth * jth.pres / (jth.dens * jth.dens) * kernel.gradW(dr, jth.smth) +
                                     AV * gradW);
                    hydro[i].eng_dot +=
                            jth.mass * (ith.grad_smth * ith.pres / (ith.dens * ith.dens) + 0.5 * AV) * dv * gradW;
#else
                    const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ith.smth) + kernel.gradW(dr, jth.smth));
                    hydro[i].acc     -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW;
                    hydro[i].eng_dot += jth.mass * (ith.pres / (ith.dens * ith.dens) + 0.5 * AV) * dv * gradW;
#endif
                }
                hydro[i].dt = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
            }
        }
    };

    void CalcMomentum(PS::ParticleSystem<STD::RealPtcl> &sph_system) {
        PS::F64vec momentum = 0; // vector#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            momentum += sph_system[i].vel; // specific momentum, i.e. normalized to mass
        }
        momentum = PS::Comm::getSum(momentum);
        if (PS::Comm::getRank() == 0) {
            printf("momentum x: %.16e\n", momentum.x);
            printf("momentum y: %.16e\n", momentum.y);
            printf("momentum z: %.16e\n", momentum.z);
            printf("momentum total: %.16e\n", momentum.x + momentum.y + momentum.z);
        }
    };

    template<class TPtclJ>
    class CalcGravityForce {
        static const double G;
        kernel_t kernel;
    public:
        void operator()(const EPI::Grav *const __restrict ep_i, const PS::S32 Nip, const TPtclJ *const __restrict ep_j,
                        const PS::S32 Njp, RESULT::Grav *const grav) {
            for (PS::S32 i = 0; i < Nip; ++i) {  // loop through ith particles
                const EPI::Grav &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) { // loop through jth particles
                    const TPtclJ &jth = ep_j[j];
                    const EPI::Grav &jth2 = ep_i[j];
                    const PS::F64vec dr = ith.pos - jth.pos;
                    const PS::F64 dr2 = dr * dr;
                    PS::F64 mj_smth = 0.0;
                    const PS::F64 abs_dist = sqrt(dr2);
                    if (abs_dist < (2 * ith.smth)) {
                        mj_smth = jth.mass;
                    } else {
                        mj_smth = (4.0 / 3.0) * 3.1415 * jth2.dens * math::pow3(abs_dist);
                    }
                    std::cout << jth2.dens << std::endl;
                    const PS::F64 dr_inv = 1.0 / sqrt(dr2 + 1e-6);
                    const PS::F64 m_dr3_inv = mj_smth * math::pow3(dr_inv);
                    grav[i].acc -= G * m_dr3_inv * dr;
                    grav[i].pot -= G * mj_smth * dr_inv;
                }
            }
        }
    };

//    template<class TPtclJ>
//    class CalcGravityForce {
//        static const double G;
//    public:
//        void operator()(const EPI::Grav *const __restrict ep_i, const PS::S32 Nip, const TPtclJ *const __restrict ep_j,
//                        const PS::S32 Njp, RESULT::Grav *const grav) {
//            for (PS::S32 i = 0; i < Nip; ++i) {
//                const EPI::Grav &ith = ep_i[i];
//                for (PS::S32 j = 0; j < Njp; ++j) {
//                    const TPtclJ &jth = ep_j[j];
//                    const PS::F64vec dr = ith.pos - jth.pos;
//                    const PS::F64 dr2 = dr * dr;
//                    const PS::F64 dr_inv = 1.0 / sqrt(dr2 + 1e-6);
//                    const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
//                    grav[i].acc -= G * m_dr3_inv * dr;
//                    grav[i].pot -= G * jth.mass * dr_inv;
//                }
//            }
//        }
//    };

    template<class TPtclJ>
    const double CalcGravityForce<TPtclJ>::G = 6.67e-11;
}

namespace DI {
    class CalcDensity {
        kernel_t kernel;
    public:
        void operator()(const EPI::Dens *const ep_i, const PS::S32 Nip, const EPJ::Dens *const ep_j, const PS::S32 Njp,
                        RESULT::Dens *const dens) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                const EPI::Dens &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const EPJ::Dens &jth = ep_j[j];
                    const PS::F64vec dr = jth.pos - ith.pos;
                    dens[i].dens_smth += jth.mass * kernel.W(dr, ith.smth);
                    dens[i].pres_smth += jth.pV * kernel.W(dr, ith.smth);
                }
                dens[i].dens_smth = std::max(100.0, dens[i].dens_smth);
                dens[i].smth = PARAM::SMTH * pow(ith.mass / dens[i].dens_smth, 1.0 / (PS::F64) (PARAM::Dim));
            }
        }
    };

    void CalcPressure(PS::ParticleSystem<DI::RealPtcl> &sph_system) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            sph_system[i].snds = sph_system[i].EoS->SoundSpeed(sph_system[i].dens, sph_system[i].eng);
        }
    }

    void CalcMomentum(PS::ParticleSystem<STD::RealPtcl> &sph_system) {
        PS::F64vec momentum = 0; // vector#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            momentum += sph_system[i].vel; // specific momentum, i.e. normalized to mass
        }
        momentum = PS::Comm::getSum(momentum);
        if (PS::Comm::getRank() == 0) {
            printf("momentum x: %.16e\n", momentum.x);
            printf("momentum y: %.16e\n", momentum.y);
            printf("momentum z: %.16e\n", momentum.z);
            printf("momentum total: %.16e\n", momentum.x + momentum.y + momentum.z);
        }
    };

    class CalcDerivative {
        kernel_t kernel;
    public:
        void operator()(const EPI::Drvt *ep_i, const PS::S32 Nip, const EPJ::Drvt *ep_j, const PS::S32 Njp,
                        RESULT::Drvt *const drvt) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                const EPI::Drvt &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const EPJ::Drvt &jth = ep_j[j];
                    const PS::F64vec dr = ith.pos - jth.pos;
                    const PS::F64vec dv = ith.vel - jth.vel;
                    drvt[i].div_v += -jth.pV * dv * kernel.gradW(dr, ith.smth);
                    drvt[i].rot_v += -jth.pV * dv ^ kernel.gradW(dr, ith.smth);
                    drvt[i].grad_smth -= jth.mass / ith.smth *
                                         (PARAM::Dim * kernel.W(dr, ith.smth) + dr * kernel.gradW(dr, ith.smth));
                }
                drvt[i].grad_smth = 1.0 / (1.0 + ith.smth * drvt[i].grad_smth / (PARAM::Dim * ith.dens_smth));
                drvt[i].div_v *= 1.0 / ith.pres_smth;
                drvt[i].rot_v *= 1.0 / ith.pres_smth;
            }
        }
    };

    class CalcHydroForce {
        const kernel_t kernel;
    public:
        void
        operator()(const EPI::Hydro *const ep_i, const PS::S32 Nip, const EPJ::Hydro *const ep_j, const PS::S32 Njp,
                   RESULT::Hydro *const hydro) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                PS::F64 v_sig_max = 0.0;
                const EPI::Hydro &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const EPJ::Hydro &jth = ep_j[j];
                    const PS::F64vec dr = ith.pos - jth.pos;
                    const PS::F64vec dv = ith.vel - jth.vel;
                    const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
                    const PS::F64 v_sig = ith.snds + jth.snds - 3.0 * w_ij;
                    v_sig_max = std::max(v_sig_max, v_sig);
                    PS::F64 AV =
                            -0.5 * v_sig * w_ij / (0.5 * (ith.dens_smth + jth.dens_smth)) * 0.5 * (ith.Bal + jth.Bal);
                    if (PARAM::FLAG_R00 == true) {
                        AV *= 0.5 * (ith.AVa + jth.AVa);
                    }
#if 1
                    const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ith.smth) * ith.grad_smth +
                                                    kernel.gradW(dr, jth.smth) * jth.grad_smth);
                    hydro[i].acc -= ith.pV * jth.pV / ith.mass * (ith.grad_smth * kernel.gradW(dr, ith.smth) *
                                                                  pow(ith.pres_smth, 1.0 / PARAM::DISPH_POWER - 2.0) +
                                                                  jth.grad_smth * kernel.gradW(dr, jth.smth) *
                                                                  pow(jth.pres_smth, 1.0 / PARAM::DISPH_POWER - 2.0)) +
                                    jth.mass * AV * gradW;
                    hydro[i].eng_dot += ith.pV * jth.pV / ith.mass * ith.grad_smth *
                                        pow(ith.pres_smth, 1.0 / PARAM::DISPH_POWER - 2.0) * (dv * gradW) +
                                        0.5 * jth.mass * AV * dv * gradW;
#else
                    //#warning untested...
                    const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ith.smth) + kernel.gradW(dr, jth.smth));
                    hydro[i].acc     -= ith.pV * jth.pV / ith.mass * (pow(ith.pres_smth, 1.0 / PARAM::DISPH_POWER - 2.0) + pow(jth.pres_smth, 1.0 / PARAM::DISPH_POWER - 2.0)) * gradW + jth.mass * AV * gradW;
                    hydro[i].eng_dot += ith.pV * jth.pV / ith.mass * pow(ith.pres_smth, 1.0 / PARAM::DISPH_POWER - 2.0) * (dv * gradW) + 0.5 * jth.mass * AV * (dv * gradW);
#endif
                }
                hydro[i].dt = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
            }
        }
    };

    template<class TPtclJ>
    class CalcGravityForce {
        static const double G;
    public:
        void operator()(const EPI::Grav *const __restrict ep_i, const PS::S32 Nip, const TPtclJ *const __restrict ep_j,
                        const PS::S32 Njp, RESULT::Grav *const grav) {
            for (PS::S32 i = 0; i < Nip; ++i) {
                const EPI::Grav &ith = ep_i[i];
                for (PS::S32 j = 0; j < Njp; ++j) {
                    const TPtclJ &jth = ep_j[j];
                    const PS::F64vec dr = ith.pos - jth.pos;  // the distance between two particles
                    const PS::F64 dr2 = dr * dr;  // dr^2
                    const PS::F64 dr_inv = 1.0 / sqrt(dr2 + 1.0e-6);  // 1 / (dr^2 + A^2)
                    const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv); // (m / (dr^2 + A^2)^3)
                    grav[i].acc -= G * m_dr3_inv * dr;
                    grav[i].pot -= G * jth.mass * dr_inv;
                }
            }
        }
    };

    template<class TPtclJ>
    const double CalcGravityForce<TPtclJ>::G = 6.67e-11;
}

