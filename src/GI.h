#pragma once

#include "parse.h"
#include <unordered_set>

#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif

std::unordered_set<unsigned int> create_removal_list(const unsigned int lowest_index, const unsigned int highest_index,
                                                     const unsigned int number_of_removed_items) {
    std::unordered_set<unsigned int> removal_list;

    if (number_of_removed_items == 0)
        return removal_list;

    while (removal_list.size() < number_of_removed_items) {
        const unsigned int num = rand() % static_cast<unsigned int>(highest_index - lowest_index) + lowest_index;

        // This works even if num is already in the list, because the unordered_set filters out duplicates
        removal_list.insert(num);
    }
    return removal_list;
}


template<class Ptcl>
class GI_universal : public Problem<Ptcl> {
public:

    static void setEoS(PS::ParticleSystem<Ptcl> &sph_system) {
#pragma omp parallel for
        for (PS::U64 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag % 2 == 0) {
                sph_system[i].setPressure(&ADunite);
            } else {
                sph_system[i].setPressure(&Iron);
            }
        }
    }

    static double end_time;
    static double damping;
//    static std::string impactor_fname;
//    static std::string target_fname;

    static void setupIC(PS::ParticleSystem<Ptcl> &sph_system, system_t &sysinfo, PS::DomainInfo &dinfo,
                        ParameterFile &parameter_file) {
        const unsigned int mode = parameter_file.getValueOf("mode", 2);

        if (mode == 1) {
            end_time = parameter_file.getValueOf("end_time", 100000.0);
            static constexpr PS::F64 R = 6400.0e+3;
            static constexpr PS::F64 M = 6.0e+24;
            static constexpr double Grav = 6.67e-11;
            static constexpr double L_EM = 3.5e+34;
            damping = parameter_file.getValueOf("damping", 1.0);
//            impactor_fname = parameter_file.getValueOf("impactor_fname", "imp.dat");
//            target_fname = parameter_file.getValueOf("target_fname", "tar.dat");

            if (PS::Comm::getRank() != 0) return;
            std::vector<Ptcl> tar, imp;
            {
                std::ifstream fin("input/tar.dat");
                double time;
                std::size_t N;
                fin >> time;
                fin >> N;
                while (!fin.eof()) {
                    Ptcl ith;
                    fin >> ith.id >> ith.tag >> ith.mass >> ith.pos.x >> ith.pos.y >> ith.pos.z >> ith.vel.x
                        >> ith.vel.y
                        >> ith.vel.z >> ith.dens >> ith.eng >> ith.pres >> ith.pot >> ith.ent >> ith.temp;
                    tar.push_back(ith);
                }
                tar.pop_back();
            }
            {
                std::ifstream fin("input/imp.dat");
                double time;
                std::size_t N;
                fin >> time;
                fin >> N;
                while (!fin.eof()) {
                    Ptcl ith;
                    fin >> ith.id >> ith.tag >> ith.mass >> ith.pos.x >> ith.pos.y >> ith.pos.z >> ith.vel.x
                        >> ith.vel.y
                        >> ith.vel.z >> ith.dens >> ith.eng >> ith.pres >> ith.pot >> ith.ent >> ith.temp;
                    imp.push_back(ith);
                }
                imp.pop_back();
            }
            {
                sph_system.setNumberOfParticleLocal(tar.size() + imp.size());
                std::size_t cnt = 0;
                for (int i = 0; i < tar.size(); ++i) {
                    sph_system[cnt] = tar[i];
                    ++cnt;
                }
                PS::F64 num_particles_tar = tar.size();
                // target: id < 2; impactor: id >= 2
                for (int i = 0; i < imp.size(); ++i) {
                    //ad-hoc modification
                    imp[i].tag += 2; // impactor particles are tagged so that mantle and core particles take id values 2 and 3 (i.e. values 0 and 1 assigned to target)
                    imp[i].id += num_particles_tar;
                    sph_system[cnt] = imp[i];
                    ++cnt;
                }
            }


            //Shift positions
            PS::F64vec pos_tar = 0; // position of the target
            PS::F64vec pos_imp = 0; // position of the impactor
            PS::F64 mass_tar = 0; // mass of the target
            PS::F64 mass_imp = 0; // mass of the impactor
            PS::F64 mass_total = 0; // total mass (mass impactor + mass target)
            PS::F64 radi_tar = 0; // radius of the target
            PS::F64 radi_imp = 0; // radius of the impacor
            for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                if (sph_system[i].tag <= 1) {
                    //target
                    pos_tar += sph_system[i].mass * sph_system[i].pos;
                    mass_tar += sph_system[i].mass;
                } else {
                    //impactor
                    pos_imp += sph_system[i].mass * sph_system[i].pos;
                    mass_imp += sph_system[i].mass;
                }
            }
            //accumulate
            mass_total = mass_tar + mass_imp;
            pos_tar /= mass_tar;
            pos_imp /= mass_imp;
            for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                if (sph_system[i].tag <= 1) {
                    //target
                    radi_tar = std::max(radi_tar, sqrt((pos_tar - sph_system[i].pos) * (pos_tar - sph_system[i].pos)));
                } else {
                    //impactor
                    radi_imp = std::max(radi_imp, sqrt((pos_imp - sph_system[i].pos) * (pos_imp - sph_system[i].pos)));
                }
            }

//            const double x_init = radi_tar + radi_imp + parameter_file.getValueOf("delta_x", 0.0);
            double input = parameter_file.getValueOf("L_init_vs_L_em", 0.10);
            const double L_init = L_EM * input;
            PS::F64 v_imp = parameter_file.getValueOf("impVel", 0.0); // the impact velocity

            PS::F64 v_esc = sqrt(2.0 * Grav * (mass_tar + mass_imp) / (radi_tar + radi_imp));

            // the momentum balance must be M_imp  * M_tar/Mtotal * Vimp - M_tar * M_imp/Mtotal * Vimp = 0
            // where p_imp = M_imp  * M_tar/Mtotal * Vimp
            // and p_tar = M_tar * M_imp/Mtotal * Vimp
            // such that p_imp - p_tar = 0
            PS::F64 v_impactor = (mass_tar / mass_total) * v_imp;
            PS::F64 v_target = (mass_imp / mass_total) * v_imp;

            PS::F64 impAngle =
                    parameter_file.getValueOf("impact_angle", 0.0) / 180.0 * math::pi; //converting from degree to radian

//            const double v_inf = sqrt(std::max(v_imp * v_imp - v_esc * v_esc, 0.0));
            PS::F64 x_init = cos(impAngle) * (radi_imp + radi_tar);
            PS::F64 y_init = sin(impAngle) * (radi_imp + radi_tar);

//            std::cout << "v_init = " << v_init << std::endl;
            std::cout << "y_init / Rtar = " << y_init / radi_tar << std::endl;
//            std::cout << "v_imp  = " << v_imp << " = " << v_imp / v_esc << "v_esc" << std::endl;
            std::cout << "m_imp  = " << mass_imp / M << std::endl;
            //shift'em
            for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                // target particles
                sph_system[i].vel.x = 0.0;
                sph_system[i].vel.y = 0.0;
                sph_system[i].vel.z = 0.0;
                if (sph_system[i].tag >= 2) { // this currently sets parameters for impactor-tagged particles
                    sph_system[i].pos -= pos_imp;
                    sph_system[i].pos.x += x_init;
                    sph_system[i].pos.y += y_init;
//                    sph_system[i].vel.x -= v_init;
                    sph_system[i].vel.x -= v_impactor;
                } else { // this currently sets parameters for target-tagged particles
                    sph_system[i].vel.x += v_target;
                }
            }

            //
//            const double b = L_init / v_inf / mass_imp;
//            const double a = -Grav * mass_tar / (v_inf * v_inf);
//            const double rp = (sqrt(a * a + b * b) + a);
//            const double r_imp = 2.0 / (v_imp * v_imp / (Grav * mass_tar) + 1.0 / a);
            std::cout << "target radius: " << radi_tar << std::endl;
            std::cout << "impactor radius: " << radi_imp << std::endl;
            std::cout << "impactor initial x: " << x_init << std::endl;
            std::cout << "impactor initial y: " << y_init << std::endl;
//            std::cout << "a = " << a / radi_tar << " R_tar" << std::endl;
//            std::cout << "b = " << b / radi_tar << " R_tar" << std::endl;
//            std::cout << "r_imp = " << r_imp / radi_tar << " R_tar" << std::endl;
//            std::cout << "r_p   = " << rp / radi_tar << " R_tar" << std::endl;
            std::cout << "angle = " << impAngle << std::endl;

            std::cout << "setup..." << std::endl;






        } else if (mode == 2) {
            const double Corr = .98;//Correction Term

            /////////
            //place ptcls
            /////////
            std::vector<Ptcl> ptcl;
            std::vector<Ptcl> tar;//Target
            std::vector<Ptcl> imp;//Impactor
            /////////

            // Use parameters from input file, or defaults if none provided

            PS::F64 UnitMass = parameter_file.getValueOf("UnitMass", 6.0e+24);
            PS::F64 UnitRadi = parameter_file.getValueOf("UnitRadi", 6400e+3);
            PS::F64 coreFracRadi = parameter_file.getValueOf("coreFracRadi", 3500.0e+3 / 6400.0e+3);
            PS::F64 coreFracMass = parameter_file.getValueOf("coreFracMass", 0.3);
            PS::F64 imptarMassRatio = parameter_file.getValueOf("imptarMassRatio", 0.1);

            PS::F64 impVel = parameter_file.getValueOf("impVel", 0.);
            PS::F64 impAngle =
                    parameter_file.getValueOf("impact_angle", 0.) / 180.0 * math::pi; //converting from degree to radian


            end_time = parameter_file.getValueOf("end_time", 1.0e+4);
            damping = parameter_file.getValueOf("damping", 1.);
            PS::F64 Nptcl = parameter_file.getValueOf("total_number_of_particles", 100000);


            const PS::F64 Expand = 1.1;
            const PS::F64 tarMass = UnitMass;
            const PS::F64 tarRadi = UnitRadi;
            const PS::F64 tarCoreMass = tarMass * coreFracMass;
            const PS::F64 tarCoreRadi = tarRadi * coreFracRadi;
            const PS::F64 impMass = imptarMassRatio * tarMass;
            const PS::F64 impRadi = Expand * cbrt(impMass / tarMass) * UnitRadi;
            const PS::F64 impCoreMass = impMass * coreFracMass;
            const PS::F64 impCoreRadi = impRadi * coreFracRadi;


            const double offset = 5.0 * UnitRadi;
            /* the following line predicts the number of grid points in one direction
              The volume of a recutangular box whose radius is L^3
              The volume of a sphere whose radius is L/2 is 4 \pi/3 (L/2)^3
              dx is defined as the grid size in one direction
              Nptcl = (volume of a sphere)/(dx)^3 * L^3, where L = 2.0
                  dx = (4.0/3.0 * math::pi/Nptcl)^{1/3}
              The number of grid point is 2.0/dx
              we multiply by 1.1 so that enough particles are created to generate a sphere */
            const int gridpoint = int(2.0 / pow(4.0 / 3.0 * math::pi * 1.1 / Nptcl, 0.333));
            const PS::F64 dx = 2.0 / gridpoint;
            const PS::F64 Grav = 6.67e-11;

            //target
            int tarNptcl = 0;
            int tarNmntl = 0;
            int tarNcore = 0;


            //impactor
            double tarCoreShrinkFactor = 1.0;
            int impNmntl = 0;
            int impNcore = 0;
            int impNptcl = 0;

            const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();

            PS::S32 id = 0;

            for (PS::F64 x = -1.0; x <= 1.0; x += dx) {
                for (PS::F64 y = -1.0; y <= 1.0; y += dx) {
                    for (PS::F64 z = -1.0; z <= 1.0; z += dx) {
                        const PS::F64 r = sqrt(x * x + y * y + z * z) * UnitRadi;
                        if (r >= tarRadi || r <= tarCoreRadi) continue;
                        ++tarNmntl;
                    }
                }
            }

            while (tarCoreShrinkFactor *= 0.99) {
                tarNcore = 0;
                for (PS::F64 x = -1.0; x <= 1.0; x += dx) {
                    for (PS::F64 y = -1.0; y <= 1.0; y += dx) {
                        for (PS::F64 z = -1.0; z <= 1.0; z += dx) {
                            const PS::F64 r = tarCoreShrinkFactor * sqrt(x * x + y * y + z * z) * UnitRadi;
                            if (r >= Corr * tarCoreRadi) continue;
                            ++tarNcore;
                        }
                    }
                }
                if ((double) (tarNcore) / (double) (tarNcore + tarNmntl) > coreFracMass) break;
            }

            ///////////////////
            //Dummy end
            ///////////////////



            // checking if there are enough mantle particles
            if (tarNmntl < static_cast<int>(Nptcl * (1.0 - coreFracMass))) {
                std::cout
                        << "Too few mantle particles. Increase the grid size. The easiest fix is to increase the gridpoint in GI.h "
                        << std::endl;
                exit(0);
            }

            // checking if there are enough core particles
            if (tarNcore < static_cast<int>(Nptcl * coreFracMass)) {
                std::cout << "Too few core particles. Increase the grid size and/or change the core shrink factor."
                          << std::endl;
                exit(0);
            }

            //removing particles to reach the exact Nptcl

            int index = 0;

            std::unordered_set<unsigned int> removal_list;
            removal_list = create_removal_list(0, tarNmntl,
                                               tarNmntl - static_cast<int>(Nptcl * (1.0 - coreFracMass)));

            std::cout << "creating target" << std::endl;
            for (PS::F64 x = -1.0; x <= 1.0; x += dx) {
                for (PS::F64 y = -1.0; y <= 1.0; y += dx) {
                    for (PS::F64 z = -1.0; z <= 1.0; z += dx) {
                        const PS::F64 r = sqrt(x * x + y * y + z * z) * UnitRadi;
                        if (r >= tarRadi || r <= tarCoreRadi) continue;
                        Ptcl ith;
                        ith.pos.x = UnitRadi * x;
                        ith.pos.y = UnitRadi * y;
                        ith.pos.z = UnitRadi * z;
                        ith.dens = (tarMass - tarCoreMass) / (4.0 / 3.0 * math::pi * (tarRadi * tarRadi * tarRadi -
                                                                                      tarCoreRadi * tarCoreRadi *
                                                                                      tarCoreRadi));
                        ith.mass = tarMass;
                        ith.eng = 0.1 * Grav * tarMass / tarRadi;
                        ith.id = id++;
                        ith.setPressure(&ADunite);
                        ith.tag = 0;

                        if (removal_list.count(index)) {
                            id += -1;
                        } else if (ith.id / NptclIn1Node == PS::Comm::getRank()) {
                            tar.push_back(ith);
                        }
                        index += 1;
                    }
                }
            }

            std::cout << "# of mantle particles = " << id << std::endl;

            // making the core condition
            removal_list.clear();
            removal_list = create_removal_list(tarNmntl, tarNmntl + tarNcore,
                                               tarNcore - static_cast<int>(Nptcl * coreFracMass));

            index = tarNmntl;


            for (PS::F64 x = -1.0; x <= 1.0; x += dx) {
                for (PS::F64 y = -1.0; y <= 1.0; y += dx) {
                    for (PS::F64 z = -1.0; z <= 1.0; z += dx) {
                        const PS::F64 r = tarCoreShrinkFactor * sqrt(x * x + y * y + z * z) * UnitRadi;
                        if (r >= Corr * tarCoreRadi) continue;
                        Ptcl ith;
                        ith.pos.x = tarCoreShrinkFactor * UnitRadi * x;
                        ith.pos.y = tarCoreShrinkFactor * UnitRadi * y;
                        ith.pos.z = tarCoreShrinkFactor * UnitRadi * z;
                        ith.dens = tarCoreMass /
                                   (4.0 / 3.0 * math::pi * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr *
                                    Corr);
                        ith.mass = tarMass;
                        ith.eng = 0.1 * Grav * tarMass / tarRadi;
                        ith.id = id++;
                        ith.setPressure(&Iron);
                        ith.tag = 1;

                        if (removal_list.count(index)) {
                            id += -1;
                        } else {
                            if (ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
                        }
                        index += 1;

                    }
                }
            }

            std::cout << "# of total particles = " << id << std::endl;

            tarNmntl = static_cast<int>(Nptcl * (1.0 - coreFracMass));
            tarNcore = static_cast<int>(Nptcl * coreFracMass);

            for (PS::U32 i = 0; i < tar.size(); ++i) {
                tar[i].mass /= (PS::F64) (Nptcl);
            }
            for (PS::U32 i = 0; i < tar.size(); ++i) {
                ptcl.push_back(tar[i]);
            }

            tarNptcl = tarNcore + tarNmntl;
            impNptcl = impNcore + impNmntl;

            std::cout << "Target  :" << tarNptcl << std::endl;
            std::cout << "    radius           : " << tarRadi << std::endl;
            std::cout << "    total-to-core    : " << (double) (tarNcore) / (double) (tarNptcl) << std::endl;
            std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
            std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
            std::cout << "    core density     : "
                      << tarCoreMass /
                         (4.0 * math::pi / 3.0 * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr)
                      << std::endl;
            std::cout << "    mantle density   : " << (tarMass - tarCoreMass) / (4.0 * math::pi / 3.0 *
                                                                                 (tarRadi * tarRadi * tarRadi -
                                                                                  tarCoreRadi * tarCoreRadi *
                                                                                  tarCoreRadi))
                      << std::endl;
            std::cout << "    mean density     : " << tarMass / (4.0 * math::pi / 3.0 * tarRadi * tarRadi * tarRadi)
                      << std::endl;

            //assert(Nptcl == tarNptcl + impNptcl);

            std::cout << "Total number of particles:" << Nptcl << std::endl;

            const PS::S32 numPtclLocal = ptcl.size();
            sph_system.setNumberOfParticleLocal(numPtclLocal);
            for (PS::U32 i = 0; i < ptcl.size(); ++i) {
                sph_system[i] = ptcl[i];
            }
            //Fin.
            std::cout << "# of ptcls = " << ptcl.size() << std::endl;
            std::cout << "setup..." << std::endl;
        }
    }

    static void postTimestepProcess(PS::ParticleSystem<Ptcl>& sph_system, system_t& sys){
        if(1){
            //Shift Origin
            PS::F64vec com_loc;  //center of mass
            PS::F64vec mom_loc;  //moment
            PS::F64vec mom_ang_loc; // angular momentum
            PS::F64 mass_loc = 0;  // mass
            PS::F64 eng_loc = 0;  // energy
            for (PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                com_loc += sph_system[i].pos * sph_system[i].mass;
                mom_loc += sph_system[i].vel * sph_system[i].mass;
                mass_loc += sph_system[i].mass;
                eng_loc += sph_system[i].mass * (sph_system[i].eng + sph_system[i].vel * sph_system[i].vel + sph_system[i].pot);
            }
            PS::F64vec com = PS::Comm::getSum(com_loc);
            PS::F64vec mom = PS::Comm::getSum(mom_loc);
            PS::F64 mass = PS::Comm::getSum(mass_loc);
            com /= mass;
            PS::F64 eng = PS::Comm::getSum(eng_loc);
            PS::F64 mom2 = mom * mom;
            PS::F64 total_mom = sqrt(mom2);
            for (PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i) {
                PS::F64 x = sph_system[i].pos.x - com.x;
                PS::F64 y = sph_system[i].pos.y - com.y;
                PS::F64 z = sph_system[i].pos.z - com.x;
                mom_ang_loc.x += sph_system[i].mass * ((y * sph_system[i].vel.z) - (z * sph_system[i].vel.y));
                mom_ang_loc.y += sph_system[i].mass * ((x * sph_system[i].vel.z) - (z * sph_system[i].vel.x));
                mom_ang_loc.z += sph_system[i].mass * ((x * sph_system[i].vel.y) - (y * sph_system[i].vel.x));
            }
            PS::F64vec mom_ang = PS::Comm::getSum(mom_ang_loc);
            PS::F64 mom_ang2 = mom_ang * mom_ang;
            PS::F64 total_mom_ang = sqrt(mom_ang2);
            if (PS::Comm::getRank() == 0) {
                std::cout << "Linear Momentum: " << mom << std::endl;
                std::cout << "Total Linear Momentum: " << total_mom << std::endl;
                std::cout << "Angular Momentum: " << mom_ang << std::endl;
                std::cout << "Total Angular Momentum: " << total_mom_ang << std::endl;
                std::cout << "Energy: " << eng << std::endl;
            }
        }
    }
};