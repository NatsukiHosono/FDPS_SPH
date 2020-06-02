#pragma once

#include "EoS.h"

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
                dens[i].dens = std::max(5.0, dens[i].dens);
#endif
                dens[i].smth = PARAM::SMTH * pow(ith.mass / dens[i].dens, 1.0 / (PS::F64) (PARAM::Dim));
            }
        }
    };

    class AngularVelocity {
    private:
        static std::pair<PS::F64, PS::F64> velocity_xy(PS::F64 &position_x, PS::F64 &position_y, PS::F64 &angular_velocity, PS::F64 &dt) {
            PS::F64 v_x = -(position_x * angular_velocity * sin(angular_velocity * dt)) - (position_y * angular_velocity * cos(angular_velocity * dt));
            PS::F64 v_y = (position_x * angular_velocity * cos(angular_velocity * dt)) - (position_y * angular_velocity * sin(angular_velocity * dt));
            return std::make_pair(v_x, v_y);
        };
    public:
        static void add_angular_velocity_xy(PS::ParticleSystem<STD::RealPtcl> &sph_system, PS::F64 angular_velocity, PS::F64 &dt) {
            if (angular_velocity != 0) {
#pragma omp parallel for
                for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                    std::pair<PS::F64, PS::F64> v_x_v_y = velocity_xy(sph_system[i].pos.x, sph_system[i].pos.y,
                                                                      angular_velocity, dt);
                    sph_system[i].vel.x = v_x_v_y.first;
                    sph_system[i].vel.y = v_x_v_y.second;
                }
            }
        }
    };

    void CalcPressure(PS::ParticleSystem<STD::RealPtcl> &sph_system, unsigned int iron_grid_size,
            unsigned int silicate_grid_size) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag == 0) {
                sph_system[i].pres = sph_system[i].EoS->Pressure(sph_system[i].dens, sph_system[i].eng, silicate_grid_size);
                sph_system[i].snds = sph_system[i].EoS->SoundSpeed(sph_system[i].dens, sph_system[i].eng, silicate_grid_size);
            } else {
                sph_system[i].pres = sph_system[i].EoS->Pressure(sph_system[i].dens, sph_system[i].eng, iron_grid_size);
                sph_system[i].snds = sph_system[i].EoS->SoundSpeed(sph_system[i].dens, sph_system[i].eng, iron_grid_size);
            }
            if (sph_system[i].pres < 0) {
                sph_system[i].pres = 0;
            }
        }
    }

//    void CalcInternalEnergy(PS::ParticleSystem<STD::RealPtcl> &sph_system,
//                            const unsigned int aneos_grid_size, const unsigned int tillotson_grid_size) {
//#pragma omp parallel for
//        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//            // switch to tillotson for iron using id tag (it already is doing this?)
//            // need to use the iron table for interpolating against iron
//            sph_system[i].eng = sph_system[i].EoS->InternalEnergy(sph_system[i].dens, sph_system[i].ent,
//                                                                  aneos_grid_size);
//        }
//    }
//
//    void CalcEntropy(PS::ParticleSystem<STD::RealPtcl> &sph_system) {
//#pragma omp parallel for
//        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//            sph_system[i].ent = sph_system[i].EoS->Entropy(sph_system[i].dens, sph_system[i].eng);
//        }
//    }

    void SetConstantEntropy(PS::ParticleSystem<STD::RealPtcl> &sph_system, const double &mantle_entropy,
                            const double &core_entropy) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag == 0) {
                sph_system[i].ent = mantle_entropy;
            } else {
                sph_system[i].ent = core_entropy;
            }
        }
    }

    void CalcEntropy(PS::ParticleSystem<STD::RealPtcl> &sph_system, unsigned int iron_grid_size,
                      unsigned int silicate_grid_size) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag == 0) {
                sph_system[i].ent = sph_system[i].EoS->Entropy(sph_system[i].dens, sph_system[i].eng, silicate_grid_size);
            } else {
                sph_system[i].ent = sph_system[i].EoS->Entropy(sph_system[i].dens, sph_system[i].eng, iron_grid_size);
            }
        }
    }


//    void CalcEntropyAndInternalEnergy(PS::ParticleSystem<STD::RealPtcl> &sph_system, const double &entropy,
//                                      const unsigned int aneos_grid_size, const unsigned int tillotson_grid_size) {
//#pragma omp parallel for
//        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//            sph_system[i].eng = sph_system[i].EoS->InternalEnergy(sph_system[i].dens, entropy, aneos_grid_size);
//            sph_system[i].ent = sph_system[i].EoS->Entropy(sph_system[i].dens, sph_system[i].eng);
//        }
//    }

    void CalcAll(PS::ParticleSystem<STD::RealPtcl> &sph_system, unsigned int iron_grid_size,
    unsigned int silicate_grid_size) {
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            // switch to tillotson for iron using id tag
            // if 0, then mantle, if 1, then core
            // need to use the iron table for interpolating against iron
            // hack, fix later
            if (sph_system[i].tag == 0) {
                sph_system[i].eng = sph_system[i].EoS->InternalEnergy(sph_system[i].dens, sph_system[i].ent,
                                                                      silicate_grid_size);
                 sph_system[i].temp = sph_system[i].EoS->Temperature(sph_system[i].dens, sph_system[i].eng,
                                                                    silicate_grid_size);
            } else {
                sph_system[i].eng = sph_system[i].EoS->InternalEnergy(sph_system[i].dens, sph_system[i].ent,
                                                                      iron_grid_size);
                 sph_system[i].temp = sph_system[i].EoS->Temperature(sph_system[i].dens, sph_system[i].eng,
                                                                    iron_grid_size);
            }
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
                    const PS::F64vec dr = ith.pos - jth.pos;
                    const PS::F64 dr2 = dr * dr;
                    const PS::F64 dr_inv = 1.0 / sqrt(dr2 + ith.getEps2());
                    const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
                    grav[i].acc -= G * m_dr3_inv * dr;
                    grav[i].pot -= G * jth.mass * dr_inv;
                }
            }
        }
    };

    template<class TPtclJ>
    const double CalcGravityForce<TPtclJ>::G = 6.67e-11;
}

