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
    static double end_time;
    static double damping;
    static constexpr PS::F64 R = 6400.0e+3;
    static constexpr PS::F64 M = 6.0e+24;
    static constexpr double Grav = 6.67e-11;
    static constexpr double L_EM = 3.5e+34;

    static void setupIC(PS::ParticleSystem<Ptcl> &sph_system, system_t &sysinfo, PS::DomainInfo &dinfo,
                        ParameterFile &parameter_file) {
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

        const unsigned int mode = parameter_file.getValueOf("mode", 2);
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

        switch (mode) {
            case 1:
                // This mode will enable to create a target and an imapctor from input/tar.dat and input/imp.dat
            {

                if (PS::Comm::getRank() != 0) return;
                {
                    std::ifstream fin("input/tar.dat");
                    double time;
                    std::size_t N;
                    fin >> time;
                    fin >> N;
                    while (!fin.eof()) {
                        Ptcl ith;
                        fin >> ith.id >> ith.tag >> ith.mass >> ith.pos.x >> ith.pos.y >> ith.pos.z >> ith.vel.x >> ith.vel.y
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
                        fin >> ith.id >> ith.tag >> ith.mass >> ith.pos.x >> ith.pos.y >> ith.pos.z >> ith.vel.x >> ith.vel.y
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
                    for (int i = 0; i < imp.size(); ++i) {
                        //ad-hoc modification
                        imp[i].tag += 2;
                        sph_system[cnt] = imp[i];
                        ++cnt;
                    }
                }


                //Shift positions
                PS::F64vec pos_tar = 0;
                PS::F64vec pos_imp = 0;
                PS::F64vec vel_imp = 0;
                PS::F64 mass_tar = 0;
                PS::F64 mass_imp = 0;
                PS::F64 radi_tar = 0;
                PS::F64 radi_imp = 0;
                for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                    if (sph_system[i].tag <= 1) {
                        //target
                        pos_tar += sph_system[i].mass * sph_system[i].pos;
                        mass_tar += sph_system[i].mass;
                    } else {
                        //impactor
                        pos_imp += sph_system[i].mass * sph_system[i].pos;
                        vel_imp += sph_system[i].mass * sph_system[i].vel;
                        mass_imp += sph_system[i].mass;
                    }
                }
                //accumulate
                pos_tar = PS::Comm::getSum(pos_tar);
                pos_imp = PS::Comm::getSum(pos_imp);
                vel_imp = PS::Comm::getSum(vel_imp);
                mass_tar = PS::Comm::getSum(mass_tar);
                mass_imp = PS::Comm::getSum(mass_imp);
                pos_tar /= mass_tar;
                pos_imp /= mass_imp;
                vel_imp /= mass_imp;
                for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                    if (sph_system[i].tag <= 1) {
                        //target
                        radi_tar = std::max(radi_tar,
                                            sqrt((pos_tar - sph_system[i].pos) * (pos_tar - sph_system[i].pos)));
                    } else {
                        //impactor
                        radi_imp = std::max(radi_imp,
                                            sqrt((pos_imp - sph_system[i].pos) * (pos_imp - sph_system[i].pos)));
                    }
                }
                radi_tar = PS::Comm::getMaxValue(radi_tar);
                radi_imp = PS::Comm::getMaxValue(radi_imp);

                const double v_esc = sqrt(2.0 * Grav * (mass_tar + mass_imp) / (radi_tar + radi_imp));
                const double x_init = 3.0 * radi_tar;
                double input = 0;
                input = parameter_file.getValueOf("L_init_vs_L_em", 0.10);
                const double L_init = L_EM * input;
                input = parameter_file.getValueOf("v_imp_vs_v_esc", 0.10);
                const double v_imp = v_esc * input;

                const double v_inf = sqrt(std::max(v_imp * v_imp - v_esc * v_esc, 0.0));
                double y_init = radi_tar;//Initial guess.
                double v_init;
                std::cout << "v_esc = " << v_esc << std::endl;
                for (int it = 0; it < 10; ++it) {
                    v_init = sqrt(v_inf * v_inf + 2.0 * Grav * mass_tar / sqrt(x_init * x_init + y_init * y_init));
                    y_init = L_init / (mass_imp * v_init);
                }

                std::cout << "v_init = " << v_init << std::endl;
                std::cout << "y_init / Rtar = " << y_init / radi_tar << std::endl;
                std::cout << "v_imp  = " << v_imp << " = " << v_imp / v_esc << "v_esc" << std::endl;
                std::cout << "m_imp  = " << mass_imp / M << std::endl;
                //shift'em
                for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
                    if (sph_system[i].tag <= 1) {
                    } else {
                        sph_system[i].pos -= pos_imp;
                        sph_system[i].vel -= vel_imp;
                        sph_system[i].pos.x += x_init;
                        sph_system[i].pos.y += y_init;
                        sph_system[i].vel.x -= v_init;
                    }
                }

                std::cout << "Total number of particles:" << Nptcl << std::endl;

                //
                const double b = L_init / v_inf / mass_imp;
                const double a = -Grav * mass_tar / (v_inf * v_inf);
                const double rp = (sqrt(a * a + b * b) + a);
                const double r_imp = 2.0 / (v_imp * v_imp / (Grav * mass_tar) + 1.0 / a);
                std::cout << "a = " << a / radi_tar << " R_tar" << std::endl;
                std::cout << "b = " << b / radi_tar << " R_tar" << std::endl;
                std::cout << "r_imp = " << r_imp / radi_tar << " R_tar" << std::endl;
                std::cout << "r_p   = " << rp / radi_tar << " R_tar" << std::endl;
                std::cout << "angle = " << asin(r_imp / (radi_tar + radi_imp)) / M_PI << "pi" << std::endl;

                std::cout << "setup..." << std::endl;
            }
            case 2:
                // This mode will create an initial condition
            {
                ///////////////////
                //Dummy put to determine # of ptcls
                ///////////////////


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
                            ith.mass = tarMass + impMass;
                            ith.eng = 0.1 * Grav * tarMass / tarRadi;
                            ith.id = id++;
                            ith.setPressure(&AGranite);
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
                            ith.mass = tarMass + impMass;
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
                assert(Nptcl == tarNptcl + impNptcl);

                std::cout << "Total number of particles:" << Nptcl << std::endl;

                const PS::S32 numPtclLocal = ptcl.size();
                sph_system.setNumberOfParticleLocal(numPtclLocal);
                for (PS::U32 i = 0; i < ptcl.size(); ++i) {
                    sph_system[i] = ptcl[i];
                }


                break;
            }
        }

        std::cout << "Target  :" << tarNptcl << std::endl;
        std::cout << "    radius           : " << tarRadi << std::endl;
        std::cout << "    total-to-core    : " << (double) (tarNcore) / (double) (tarNptcl) << std::endl;
        std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
        std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
        std::cout << "    core density     : "
                  << tarCoreMass / (4.0 * math::pi / 3.0 * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr)
                  << std::endl;
        std::cout << "    mantle density   : " << (tarMass - tarCoreMass) / (4.0 * math::pi / 3.0 *
                                                                             (tarRadi * tarRadi * tarRadi -
                                                                              tarCoreRadi * tarCoreRadi * tarCoreRadi))
                  << std::endl;
        std::cout << "    mean density     : " << tarMass / (4.0 * math::pi / 3.0 * tarRadi * tarRadi * tarRadi)
                  << std::endl;

        if (mode == 1) {
            std::cout << "Impactor:" << impNptcl << std::endl;
            std::cout << "    radius           : " << impRadi << std::endl;
            std::cout << "    total-to-core    : " << (double) (impNcore) / (double) (impNptcl) << std::endl;
            std::cout << "    # of core ptcls  : " << impNcore << std::endl;
            std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
            std::cout << "    core density     : " << impCoreMass /
                                                      (4.0 * math::pi / 3.0 * impCoreRadi * impCoreRadi * impCoreRadi *
                                                       Corr * Corr * Corr) << std::endl;
            std::cout << "    mantle density   : " << (impMass - impCoreMass) / (4.0 * math::pi / 3.0 *
                                                                                 (impRadi * impRadi * impRadi -
                                                                                  impCoreRadi * impCoreRadi *
                                                                                  impCoreRadi)) << std::endl;
            std::cout << "    mean density     : " << impMass / (4.0 * math::pi / 3.0 * impRadi * impRadi * impRadi)
                      << std::endl;
            std::cout << "    velocity         : " << impVel << std::endl;
            std::cout << "Tar-to-Imp mass ratio: " << (double) (impNmntl) / (double) (tarNmntl) << std::endl;
        }

        //Fin.
        std::cout << "# of ptcls = " << ptcl.size() << std::endl;
        std::cout << "setup..." << std::endl;
    }

    static void setEoS(PS::ParticleSystem<Ptcl> &sph_system) {
        for (PS::U64 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            // TODO: Modify the lines below for all particles that need new EoS
            if (sph_system[i].tag % 2 == 0) {
                sph_system[i].setPressure(&AGranite);
            } else {
                sph_system[i].setPressure(&Iron);
            }
        }
    }

    static void addExternalForce(PS::ParticleSystem<Ptcl> &sph_system, system_t &sysinfo) {
        if (sysinfo.time >= 5000) return;
        std::cout << "Add Ext. Force!!!" << std::endl;
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            sph_system[i].acc += -sph_system[i].vel * 0.05 / sph_system[i].dt;
        }
    }

//    static void postTimestepProcess(PS::ParticleSystem<Ptcl> &sph_system, system_t &sys) {
//        //SANITY Check
//        for (PS::U64 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//            sph_system[i].eng = std::max(sph_system[i].eng, 1.0e+4);
//            sph_system[i].dens = std::max(sph_system[i].dens, 100.0);
//        }
//        if (sys.step % 100 == 0 || 1) {
//            //Shift Origin
//            PS::F64vec com_loc = 0;//center of mass
//            PS::F64vec mom_loc = 0;//moment
//            PS::F64 mass_loc = 0;//
//            PS::F64 eng_loc = 0;//
//            for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//                com_loc += sph_system[i].pos * sph_system[i].mass;
//                mom_loc += sph_system[i].vel * sph_system[i].mass;
//                mass_loc += sph_system[i].mass;
//                eng_loc += sph_system[i].mass *
//                           (sph_system[i].eng + sph_system[i].vel * sph_system[i].vel + sph_system[i].pot);
//            }
//            PS::F64vec com = PS::Comm::getSum(com_loc);
//            PS::F64vec mom = PS::Comm::getSum(mom_loc);
//            PS::F64 mass = PS::Comm::getSum(mass_loc);
//            PS::F64 eng = PS::Comm::getSum(eng_loc);
//            std::cout << "Mom: " << mom << std::endl;
//            std::cout << "Eng: " << eng << std::endl;
//        }
//#if 0
//        //Shift Origin
//        PS::F64vec com_loc = 0;//center of mass of target core
//        PS::F64vec mom_loc = 0;//moment of target core
//        PS::F64 mass_loc = 0;//mass of target core
//        for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
//            if(sph_system[i].tag != 0) continue;
//            com_loc += sph_system[i].pos * sph_system[i].mass;
//            mom_loc += sph_system[i].vel * sph_system[i].mass;
//            mass_loc += sph_system[i].mass;
//        }
//        PS::F64vec com = PS::Comm::getSum(com_loc);
//        PS::F64vec mom = PS::Comm::getSum(mom_loc);
//        PS::F64 mass = PS::Comm::getSum(mass_loc);
//        com /= mass;
//        mom /= mass;
//        for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
//            sph_system[i].pos -= com;
//            sph_system[i].vel -= mom;
//        }
//#endif
//#if 1
//        std::size_t Nptcl = sph_system.getNumberOfParticleLocal();
//        for (PS::S32 i = 0; i < Nptcl; ++i) {
//            if (sqrt(sph_system[i].pos * sph_system[i].pos) / R > 150.0) {
//                //bounded particles should not be killed.
//                if (0.5 * sph_system[i].vel * sph_system[i].vel + sph_system[i].pot < 0) continue;
//                std::cout << "KILL" << std::endl;
//                sph_system[i] = sph_system[--Nptcl];
//                --i;
//            }
//        }
//        sph_system.setNumberOfParticleLocal(Nptcl);
//#endif
//    }

};


template<class Ptcl>
class GI : public Problem<Ptcl> {
public:
    static double end_time;
    static double damping;

    static void setupIC(PS::ParticleSystem<Ptcl> &sph_system, system_t &sysinfo, PS::DomainInfo &dinfo,
                        ParameterFile &parameter_file) {
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

        const unsigned int mode = parameter_file.getValueOf("mode", 2);
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

        switch (mode) {
            case 1:
                // This mode will enable to create a target and an imapctor from input/tar.dat and input/imp.dat
            {
                std::cout << "creating target from tar.dat" << std::endl;
                FILE *tarFile;
                tarFile = fopen("input/tar.txt", "r");
                FileHeader tarheader;
                int nptcltar;
                nptcltar = tarheader.readAscii(tarFile);
                std::cout << "num tar ptcl: " << nptcltar << std::endl;
                for (int i = 0; i < nptcltar; i++) {
                    Ptcl ith;
                    ith.readAscii(tarFile);
                    if (ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
                }

                for (PS::U32 i = 0; i < tar.size(); ++i) {
                    ptcl.push_back(tar[i]);
                }

                for (PS::U32 i = 0; i < tar.size(); ++i) {
                    if (tar[i].tag == 0) {
                        tarNmntl += 1;
                    } else {
                        tarNcore += 1;
                    }
                }

                tarNptcl = tarNmntl + tarNcore;

                std::cout << "creating impactor from imp.dat" << std::endl;
                FILE *impFile;
                impFile = fopen("input/imp.dat", "r");
                FileHeader impheader;
                int nptclimp;
                nptclimp = impheader.readAscii(impFile);
                std::cout << "num imp ptcl: " << nptclimp << std::endl;
                for (int i = 0; i < nptclimp; i++) {
                    Ptcl ith;
                    ith.readAscii(impFile);
                    ith.vel.x += (-1) * cos(impAngle) * impVel;
                    ith.vel.y += (-1) * sin(impAngle) * impVel;
                    ith.pos.x += (impRadi + tarRadi) * cos(impAngle);
                    ith.pos.y += (impRadi + tarRadi) * sin(impAngle);

                    if (ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
                }

                for (PS::U32 i = 0; i < imp.size(); ++i) {
                    ptcl.push_back(imp[i]);
                    if (imp[i].tag == 0) {
                        impNmntl += 1;
                    } else {
                        impNcore += 1;
                    }
                }

                impNptcl = impNmntl + impNcore;

                Nptcl = tarNptcl + impNptcl;

                break;
            }
            case 2:
                // This mode will create an initial condition
            {
                ///////////////////
                //Dummy put to determine # of ptcls
                ///////////////////


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
                            ith.mass = tarMass + impMass;
                            ith.eng = 0.1 * Grav * tarMass / tarRadi;
                            ith.id = id++;
                            ith.setPressure(&AGranite);
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
                            ith.mass = tarMass + impMass;
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
                break;
            }
        }


        tarNptcl = tarNcore + tarNmntl;
        impNptcl = impNcore + impNmntl;

        std::cout << "Target  :" << tarNptcl << std::endl;
        std::cout << "    radius           : " << tarRadi << std::endl;
        std::cout << "    total-to-core    : " << (double) (tarNcore) / (double) (tarNptcl) << std::endl;
        std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
        std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
        std::cout << "    core density     : "
                  << tarCoreMass / (4.0 * math::pi / 3.0 * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr)
                  << std::endl;
        std::cout << "    mantle density   : " << (tarMass - tarCoreMass) / (4.0 * math::pi / 3.0 *
                                                                             (tarRadi * tarRadi * tarRadi -
                                                                              tarCoreRadi * tarCoreRadi * tarCoreRadi))
                  << std::endl;
        std::cout << "    mean density     : " << tarMass / (4.0 * math::pi / 3.0 * tarRadi * tarRadi * tarRadi)
                  << std::endl;

        if (mode == 1) {
            std::cout << "Impactor:" << impNptcl << std::endl;
            std::cout << "    radius           : " << impRadi << std::endl;
            std::cout << "    total-to-core    : " << (double) (impNcore) / (double) (impNptcl) << std::endl;
            std::cout << "    # of core ptcls  : " << impNcore << std::endl;
            std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
            std::cout << "    core density     : " << impCoreMass /
                                                      (4.0 * math::pi / 3.0 * impCoreRadi * impCoreRadi * impCoreRadi *
                                                       Corr * Corr * Corr) << std::endl;
            std::cout << "    mantle density   : " << (impMass - impCoreMass) / (4.0 * math::pi / 3.0 *
                                                                                 (impRadi * impRadi * impRadi -
                                                                                  impCoreRadi * impCoreRadi *
                                                                                  impCoreRadi)) << std::endl;
            std::cout << "    mean density     : " << impMass / (4.0 * math::pi / 3.0 * impRadi * impRadi * impRadi)
                      << std::endl;
            std::cout << "    velocity         : " << impVel << std::endl;
            std::cout << "Tar-to-Imp mass ratio: " << (double) (impNmntl) / (double) (tarNmntl) << std::endl;
        }

        assert(Nptcl == tarNptcl + impNptcl);

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

    static void setEoS(PS::ParticleSystem<Ptcl> &sph_system) {
        for (PS::U64 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            // TODO: Modify the lines below for all particles that need new EoS
            if (sph_system[i].tag % 2 == 0) {
                sph_system[i].setPressure(&AGranite);
            } else {
                sph_system[i].setPressure(&Iron);
            }
        }
    }

    static void addExternalForce(PS::ParticleSystem<Ptcl> &sph_system, system_t &sysinfo) {
        if (sysinfo.time >= 5000) return;
        std::cout << "Add Ext. Force!!!" << std::endl;
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            sph_system[i].acc += -sph_system[i].vel * 0.05 / sph_system[i].dt;
        }
    }
};

template<class Ptcl>
class GI_imp : public Problem<Ptcl> {
public:
    static constexpr double end_time = 100000.0;
    static constexpr PS::F64 R = 6400.0e+3;
    static constexpr PS::F64 M = 6.0e+24;
    static constexpr double Grav = 6.67e-11;
    static constexpr double L_EM = 3.5e+34;
    static constexpr double damping = 1.0;

    static void setupIC(PS::ParticleSystem<Ptcl> &sph_system, system_t &sysinfo, PS::DomainInfo &dinfo,
                        ParameterFile &parameter_file) {
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
                fin >> ith.id >> ith.tag >> ith.mass >> ith.pos.x >> ith.pos.y >> ith.pos.z >> ith.vel.x >> ith.vel.y
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
                fin >> ith.id >> ith.tag >> ith.mass >> ith.pos.x >> ith.pos.y >> ith.pos.z >> ith.vel.x >> ith.vel.y
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
            for (int i = 0; i < imp.size(); ++i) {
                //ad-hoc modification
                imp[i].tag += 2;
                sph_system[cnt] = imp[i];
                ++cnt;
            }
        }


        //Shift positions
        PS::F64vec pos_tar = 0;
        PS::F64vec pos_imp = 0;
        PS::F64vec vel_imp = 0;
        PS::F64 mass_tar = 0;
        PS::F64 mass_imp = 0;
        PS::F64 radi_tar = 0;
        PS::F64 radi_imp = 0;
        for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag <= 1) {
                //target
                pos_tar += sph_system[i].mass * sph_system[i].pos;
                mass_tar += sph_system[i].mass;
            } else {
                //impactor
                pos_imp += sph_system[i].mass * sph_system[i].pos;
                vel_imp += sph_system[i].mass * sph_system[i].vel;
                mass_imp += sph_system[i].mass;
            }
        }
        //accumulate
        pos_tar = PS::Comm::getSum(pos_tar);
        pos_imp = PS::Comm::getSum(pos_imp);
        vel_imp = PS::Comm::getSum(vel_imp);
        mass_tar = PS::Comm::getSum(mass_tar);
        mass_imp = PS::Comm::getSum(mass_imp);
        pos_tar /= mass_tar;
        pos_imp /= mass_imp;
        vel_imp /= mass_imp;
        for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag <= 1) {
                //target
                radi_tar = std::max(radi_tar, sqrt((pos_tar - sph_system[i].pos) * (pos_tar - sph_system[i].pos)));
            } else {
                //impactor
                radi_imp = std::max(radi_imp, sqrt((pos_imp - sph_system[i].pos) * (pos_imp - sph_system[i].pos)));
            }
        }
        radi_tar = PS::Comm::getMaxValue(radi_tar);
        radi_imp = PS::Comm::getMaxValue(radi_imp);

        const double v_esc = sqrt(2.0 * Grav * (mass_tar + mass_imp) / (radi_tar + radi_imp));
        const double x_init = 3.0 * radi_tar;
        double input = 0;
        std::cout << "Input L_init / L_EM" << std::endl;
        std::cin >> input;
        const double L_init = L_EM * input;
        std::cout << "Input v_imp / v_esc" << std::endl;
        std::cin >> input;
        const double v_imp = v_esc * input;

        const double v_inf = sqrt(std::max(v_imp * v_imp - v_esc * v_esc, 0.0));
        double y_init = radi_tar;//Initial guess.
        double v_init;
        std::cout << "v_esc = " << v_esc << std::endl;
        for (int it = 0; it < 10; ++it) {
            v_init = sqrt(v_inf * v_inf + 2.0 * Grav * mass_tar / sqrt(x_init * x_init + y_init * y_init));
            y_init = L_init / (mass_imp * v_init);
        }

        std::cout << "v_init = " << v_init << std::endl;
        std::cout << "y_init / Rtar = " << y_init / radi_tar << std::endl;
        std::cout << "v_imp  = " << v_imp << " = " << v_imp / v_esc << "v_esc" << std::endl;
        std::cout << "m_imp  = " << mass_imp / M << std::endl;
        //shift'em
        for (PS::U32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag <= 1) {
            } else {
                sph_system[i].pos -= pos_imp;
                sph_system[i].vel -= vel_imp;
                sph_system[i].pos.x += x_init;
                sph_system[i].pos.y += y_init;
                sph_system[i].vel.x -= v_init;
            }
        }
        //
        const double b = L_init / v_inf / mass_imp;
        const double a = -Grav * mass_tar / (v_inf * v_inf);
        const double rp = (sqrt(a * a + b * b) + a);
        const double r_imp = 2.0 / (v_imp * v_imp / (Grav * mass_tar) + 1.0 / a);
        std::cout << "a = " << a / radi_tar << " R_tar" << std::endl;
        std::cout << "b = " << b / radi_tar << " R_tar" << std::endl;
        std::cout << "r_imp = " << r_imp / radi_tar << " R_tar" << std::endl;
        std::cout << "r_p   = " << rp / radi_tar << " R_tar" << std::endl;
        std::cout << "angle = " << asin(r_imp / (radi_tar + radi_imp)) / M_PI << "pi" << std::endl;

        std::cout << "setup..." << std::endl;
    }

//    static void postTimestepProcess(PS::ParticleSystem<Ptcl> &sph_system, system_t &sys) {
//        //SANITY Check
//        for (PS::U64 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//            sph_system[i].eng = std::max(sph_system[i].eng, 1.0e+4);
//            sph_system[i].dens = std::max(sph_system[i].dens, 100.0);
//        }
//        if (sys.step % 100 == 0 || 1) {
//            //Shift Origin
//            PS::F64vec com_loc = 0;//center of mass
//            PS::F64vec mom_loc = 0;//moment
//            PS::F64 mass_loc = 0;//
//            PS::F64 eng_loc = 0;//
//            for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//                com_loc += sph_system[i].pos * sph_system[i].mass;
//                mom_loc += sph_system[i].vel * sph_system[i].mass;
//                mass_loc += sph_system[i].mass;
//                eng_loc += sph_system[i].mass *
//                           (sph_system[i].eng + sph_system[i].vel * sph_system[i].vel + sph_system[i].pot);
//            }
//            PS::F64vec com = PS::Comm::getSum(com_loc);
//            PS::F64vec mom = PS::Comm::getSum(mom_loc);
//            PS::F64 mass = PS::Comm::getSum(mass_loc);
//            PS::F64 eng = PS::Comm::getSum(eng_loc);
//            std::cout << "Mom: " << mom << std::endl;
//            std::cout << "Eng: " << eng << std::endl;
//        }
//#if 0
//        //Shift Origin
//        PS::F64vec com_loc = 0;//center of mass of target core
//        PS::F64vec mom_loc = 0;//moment of target core
//        PS::F64 mass_loc = 0;//mass of target core
//        for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
//            if(sph_system[i].tag != 0) continue;
//            com_loc += sph_system[i].pos * sph_system[i].mass;
//            mom_loc += sph_system[i].vel * sph_system[i].mass;
//            mass_loc += sph_system[i].mass;
//        }
//        PS::F64vec com = PS::Comm::getSum(com_loc);
//        PS::F64vec mom = PS::Comm::getSum(mom_loc);
//        PS::F64 mass = PS::Comm::getSum(mass_loc);
//        com /= mass;
//        mom /= mass;
//        for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
//            sph_system[i].pos -= com;
//            sph_system[i].vel -= mom;
//        }
//#endif
//#if 1
//        std::size_t Nptcl = sph_system.getNumberOfParticleLocal();
//        for (PS::S32 i = 0; i < Nptcl; ++i) {
//            if (sqrt(sph_system[i].pos * sph_system[i].pos) / R > 150.0) {
//                //bounded particles should not be killed.
//                if (0.5 * sph_system[i].vel * sph_system[i].vel + sph_system[i].pot < 0) continue;
//                std::cout << "KILL" << std::endl;
//                sph_system[i] = sph_system[--Nptcl];
//                --i;
//            }
//        }
//        sph_system.setNumberOfParticleLocal(Nptcl);
//#endif
//    }

    static void setEoS(PS::ParticleSystem<Ptcl> &sph_system) {
#pragma omp parallel for
        for (PS::U64 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].tag % 2 == 0) {
                sph_system[i].setPressure(&Granite);
            } else {
                sph_system[i].setPressure(&Iron);
            }
        }
    }
};



