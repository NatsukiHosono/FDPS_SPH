#include <parse.h>

#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif
template <class Ptcl> class GI : public Problem<Ptcl>{
    public:
	static double end_time;
        static double damping;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo,
                        const std::string &input_file){
		const double Corr = .98;//Correction Term
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		std::vector<Ptcl> tar;//Target
		std::vector<Ptcl> imp;//Impactor
		/////////

		// Use parameters from input file, or defaults if none provided
		// TODO: Currently the input file has to be in the same directory as the executable
		//       Change this into a command-line parameter.
        std::string input_directory("input/");
        ParameterFile parameter_file(input_directory + input_file);
        std::cout << "reading from input/" << input_file << std::endl;
	PS::F64 Nptcl  = parameter_file.getValueOf("Nptcl", 100000);
		PS::F64 UnitMass = parameter_file.getValueOf("UnitMass", 6.0e+24);
		PS::F64 UnitRadi = parameter_file.getValueOf("UnitRadi", 6400e+3);
		PS::F64 coreFracRadi = parameter_file.getValueOf("coreFracRadi", 3500.0e+3 / 6400.0e+3);
		PS::F64 coreFracMass = parameter_file.getValueOf("coreFracMass", 0.3);
		PS::F64 imptarMassRatio = parameter_file.getValueOf("imptarMassRatio", 0.1);
        int mode = parameter_file.getValueOf("mode", 1 );
        PS::F64 impVel = parameter_file.getValueOf("impVel",0.);
        end_time = parameter_file.getValueOf("end_time",1.0e+4);
        damping = parameter_file.getValueOf("damping",1.);
        
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
		// the following line predicts the number of grid points in one direction
		const int  gridpoint = int(2.0/pow(4.18*1.1/Nptcl,0.333));
		const PS::F64 dx = 2.0/gridpoint;  //gridpoint;
		const PS::F64 Grav = 6.67e-11;


		int tarNmntl;
		int tarNcore;
		int tarNptcl;
		int impNptcl;
		
		int impNcore;
		int impNmntl;
		

		double tarCoreShrinkFactor = 1.0;
		double impCoreShrinkFactor = 1.0;


		
		//tarNptcl = tarNcore + tarNmntl;
		//impNptcl = impNcore + impNmntl;
		
		//if (mode==0){
		//	Nptcl    = tarNptcl + impNptcl;
		//}

		/*
		std::cout << "Target  :" << tarNptcl << std::endl;
		std::cout << "    radius           : " << tarRadi << std::endl;
		std::cout << "    total-to-core    : " << (double)(tarNcore) / (double)(tarNptcl) << std::endl;
		std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
		std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
		std::cout << "    core density     : " << tarCoreMass / (4.0 * math::pi / 3.0 * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr) << std::endl;
		std::cout << "    mantle density   : " << (tarMass - tarCoreMass) / (4.0 * math::pi / 3.0 * (tarRadi * tarRadi * tarRadi - tarCoreRadi * tarCoreRadi * tarCoreRadi)) << std::endl;
		std::cout << "    mean density     : " << tarMass / (4.0 * math::pi / 3.0 * tarRadi * tarRadi * tarRadi) << std::endl;
		std::cout << "Impactor:" << impNptcl << std::endl;
		std::cout << "    radius           : " << impRadi << std::endl;
		std::cout << "    total-to-core    : " << (double)(impNcore) / (double)(impNptcl) << std::endl;
		std::cout << "    # of core ptcls  : " << impNcore << std::endl;
		std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
		std::cout << "    core density     : " << impCoreMass / (4.0 * math::pi / 3.0 * impCoreRadi * impCoreRadi * impCoreRadi * Corr * Corr * Corr) << std::endl;
		std::cout << "    mantle density   : " << (impMass - impCoreMass) / (4.0 * math::pi / 3.0 * (impRadi * impRadi * impRadi - impCoreRadi * impCoreRadi * impCoreRadi)) << std::endl;
		std::cout << "    mean density     : " << impMass / (4.0 * math::pi / 3.0 * impRadi * impRadi * impRadi) << std::endl;
        std::cout << "    velocity         : " << impVel << std::endl;
		std::cout << "Total:" << Nptcl << std::endl;
		std::cout << "Tar-to-Imp mass ratio: " << (double)(impNmntl) / (double)(tarNmntl) << std::endl;

		*/
		
		const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();
		///////////////////
		//Real put
		///////////////////
		PS::S32 id = 0;
        
        switch (mode){
            case 1:
                std::cout << "creating target from tar.dat" << std::endl;
                FILE * tarFile;
                tarFile = fopen("input/tar.dat","r");
                FileHeader tarheader;
                tarNptcl = tarheader.readAscii(tarFile);
		
                std::cout << "Number of particles of the target: " << tarNptcl << std::endl;
		
                for(int i=0; i < tarNptcl; i++){
                    Ptcl ith;
                    ith.readAscii(tarFile);
                    if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
                }
		
                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
                    ptcl.push_back(tar[i]);
                }


		//std::cout << tar[0].mass << std::endl;
		//exit(0);
		
                std::cout << "creating impactor from imp.dat" << std::endl;
                FILE * impFile;
                impFile = fopen("input/imp.dat","r");
                FileHeader impheader;
                int impNptcl;
                impNptcl = impheader.readAscii(impFile);
                std::cout << "Number of particles of the impactor: " <<  impNptcl << std::endl;
                for(int i=0; i<impNptcl; i++){
                    Ptcl ith;
                    ith.readAscii(impFile);
                    ith.vel.x += (-1) * impVel;
                    if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
                }
                for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
                    imp[i].mass /= (PS::F64)(Nptcl);
                }
                for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
                    ptcl.push_back(imp[i]);
                }

		Nptcl = tarNptcl + impNptcl;
		
                break;
            case 2:
                //Put Tar.


	        tarNmntl = int(Nptcl * (1.0-coreFracMass));
		tarNcore = int(Nptcl * coreFracMass);

		//Assuming there is no impactor here

	      
                std::cout << "creating target" << std::endl;
                for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
                    for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
                        for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
                            const PS::F64 r = sqrt(x*x + y*y + z*z) * UnitRadi;
                            if(r >= tarRadi || r <= tarCoreRadi) continue;
                            Ptcl ith;
                            ith.pos.x = UnitRadi * x;
                            ith.pos.y = UnitRadi * y;
                            ith.pos.z = UnitRadi * z;
			    // I am a little confused by this following line. Shouldn't SPH compute density?
			    // Is this just a guess and SPH will actually compute the density from the next step? - Miki
                            ith.dens = (tarMass - tarCoreMass) / (4.0 / 3.0 * math::pi * (tarRadi * tarRadi * tarRadi - tarCoreRadi * tarCoreRadi * tarCoreRadi));
                            ith.mass = tarMass; // + impMass; 
                            ith.eng  = 0.1 * Grav * tarMass / tarRadi;
                            ith.id   = id++;
                            // TODO: Modify this line for all particles that need new EoS
                            ith.setPressure(&AGranite);
                            ith.tag = 0;
                            if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
                        }
                    }
                }

		if (tar.size() < int(Nptcl*(1.0-coreFracMass))){
		  std::cout << "The estimated number of grid points is too small. You may want to modify the parameter gridpoint in GI.h " << std::endl;
		  exit(0);
		}

		// excess of mantle particles
		int Nexcess;
		Nexcess = tar.size() - int(Nptcl*(1.0-coreFracMass));

		
		// erasing excess of mantle particles
		for(PS::U32 i = 0 ; i < Nexcess ; ++ i){
		  tar.erase(tar.begin()+ rand () % (tar.size() + 1));
		}

		
		std::cout << "# of mantle particles = " <<  tar.size() << std::endl ;

		
                for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
                    for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
                        for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
                            const PS::F64 r = tarCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
                            if(r >= Corr * tarCoreRadi) continue;
                            Ptcl ith;
                            ith.pos.x = tarCoreShrinkFactor * UnitRadi * x;
                            ith.pos.y = tarCoreShrinkFactor * UnitRadi * y;
                            ith.pos.z = tarCoreShrinkFactor * UnitRadi * z;
                            ith.dens = tarCoreMass / (4.0 / 3.0 * math::pi * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr);
                            ith.mass = tarMass; // + impMass;
                            ith.eng  = 0.1 * Grav * tarMass / tarRadi;
                            ith.id   = id++;
                            // TODO: Modify this line for all particles that need new EoS
                            ith.setPressure(&Iron);
                            ith.tag = 1;
                            if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
                        }
                    }
		}

		// excess of core particles
		Nexcess = tar.size() - Nptcl;
		if (tar.size() < int(Nptcl)){
		  std::cout << "The number of core particle is not enough. You may want to modify the parameter gridpoint in GI.h " << std::endl;
		  exit(0);
		}
		

		// erasing excess of core particles
		for(PS::U32 i = 0 ; i < Nexcess ; ++ i){
		  tar.erase(tar.begin()+ rand () % (tar.size() - int(Nptcl*(1.0-coreFracMass)) + 1) +  int(Nptcl*(1.0-coreFracMass)));
		  //tar.erase(tar.begin()+ rand () % (tar.size() + 1));
		}


		std::cout << "# of total SPH particles " <<tar.size() << std::endl ;
		exit(0);
		
	       

                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
                    tar[i].mass /= (PS::F64)(Nptcl);
                }
                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
                    ptcl.push_back(tar[i]);
                }
                break;
            case 3:
	      // Miki did not touch this -- this will probably break.
                //imp
                std::cout << "creating impactor" << std::endl;
                for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
                    for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
                        for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
                            const PS::F64 r = Expand * sqrt(x*x + y*y + z*z) * UnitRadi;
                            if(r >= impRadi || r <= impCoreRadi) continue;
                            Ptcl ith;
                            ith.pos.x = Expand * UnitRadi * x + offset;
                            ith.pos.y = Expand * UnitRadi * y;
                            ith.pos.z = Expand * UnitRadi * z;
                            ith.dens = (impMass - impCoreMass) / (4.0 / 3.0 * math::pi * (impRadi * impRadi * impRadi - impCoreRadi * impCoreRadi * impCoreRadi));
                            ith.mass = tarMass + impMass;
                            ith.eng  = 0.1 * Grav * tarMass / tarRadi;
                            ith.id   = id++;
                            // TODO: Modify this line for all particles that need new EoS
                            ith.setPressure(&AGranite);
                            ith.tag = 2;
                            if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
                        }
                    }
                }
                for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
                    for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
                        for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
                            const PS::F64 r = Expand * impCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
                            if(r >= impCoreRadi) continue;
                            Ptcl ith;
                            ith.pos.x = Expand * impCoreShrinkFactor * UnitRadi * x + offset;
                            ith.pos.y = Expand * impCoreShrinkFactor * UnitRadi * y;
                            ith.pos.z = Expand * impCoreShrinkFactor * UnitRadi * z;
                            ith.dens = impCoreMass / (4.0 / 3.0 * math::pi * impCoreRadi * impCoreRadi * impCoreRadi * Corr * Corr * Corr);
                            ith.mass = tarMass + impMass;
                            ith.eng  = 0.1 * Grav * tarMass / tarRadi;
                            ith.id   = id++;
                            // TODO: Modify this line for all particles that need new EoS
                            ith.setPressure(&Iron);
                            ith.tag = 3;
                            if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
                        }
                    }
                }
                for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
                    imp[i].mass /= (PS::F64)(Nptcl);
                }
                for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
                    ptcl.push_back(imp[i]);
                }
                break;
        }

		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			sph_system[i] = ptcl[i];
		}
		//Fin.
		std::cout << "# of ptcls = " << ptcl.size() << std::endl;
		std::cout << "setup..." << std::endl;
	}

	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			// TODO: Modify the lines below for all particles that need new EoS
			if(sph_system[i].tag % 2 == 0){
				sph_system[i].setPressure(&AGranite);
			}else{
				sph_system[i].setPressure(&Iron);
			}
		}
	}

	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo){
		if(sysinfo.time >= 5000) return;
		std::cout << "Add Ext. Force!!!" << std::endl;
		#pragma omp parallel for
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].acc += - sph_system[i].vel * 0.05 / sph_system[i].dt;
		}
	}
};
