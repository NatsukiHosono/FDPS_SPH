#include <parse.h>
#include <unordered_set>

#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif

std::unordered_set <unsigned int>  create_removal_list (const unsigned int lowest_index, const unsigned int highest_index,const unsigned int number_of_removed_items){
      std::unordered_set<unsigned int> removal_list;

      if (number_of_removed_items == 0)
	return removal_list;
      
      while (removal_list.size() < number_of_removed_items){
	const unsigned int num = rand ()  % static_cast<unsigned int>(highest_index - lowest_index) + lowest_index;
	//std::cout << num << std::endl;
	
	// This works even if num is already in the list, because the unordered_set filters out duplicates
	removal_list.insert(num);
      }
      return  removal_list;
}


template <class Ptcl> class GI : public Problem<Ptcl>{
    public:
	static double end_time;
    static double damping;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo,
                        ParameterFile &parameter_file){
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

		const unsigned int mode = parameter_file.getValueOf("mode", 2 );
		PS::F64 impVel = parameter_file.getValueOf("impVel",0.);
		PS::F64 impAngle = parameter_file.getValueOf("impact_angle",0.) /180.0 * math::pi; //converting from degree to radian

		
		end_time = parameter_file.getValueOf("end_time",1.0e+4);
		damping = parameter_file.getValueOf("damping",1.);
		PS::F64 Nptcl  = parameter_file.getValueOf("total_number_of_particles", 100000);

        
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
		  we multiply 1.1 so that the enough particles are created to generate a sphere */
		const int  gridpoint = int(2.0/pow(4.0/3.0 * math::pi * 1.1/Nptcl,0.333));
		const PS::F64 dx =  2.0/gridpoint;	
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
        
        switch (mode){
            case 1:
	      {		
                std::cout << "creating target from tar.dat" << std::endl;
                FILE * tarFile;
                tarFile = fopen("input/tar.dat","r");
                FileHeader tarheader;
                int nptcltar;
                nptcltar = tarheader.readAscii(tarFile);
                std::cout << "num tar ptcl: " << nptcltar << std::endl;
                for(int i=0; i<nptcltar; i++){
                    Ptcl ith;
                    ith.readAscii(tarFile);		    
                    if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
                }

                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
                    ptcl.push_back(tar[i]);
                }

                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		  if (tar[i].tag==0){
		      tarNmntl += 1;
		  }else{
		    tarNcore += 1;
		  }
                }

		tarNptcl = tarNmntl + tarNcore;
                
                std::cout << "creating impactor from imp.dat" << std::endl;
                FILE * impFile;
                impFile = fopen("input/imp.dat","r");
                FileHeader impheader;
                int nptclimp;
                nptclimp = impheader.readAscii(impFile);
                std::cout << "num imp ptcl: " << nptclimp << std::endl;
                for(int i=0; i<nptclimp; i++){
                    Ptcl ith;
                    ith.readAscii(impFile);

		    // This needs to be updated -- the code should determine the impactor and target sizes

                    ith.vel.x += (-1) * cos(impAngle) * impVel;
		    ith.vel.y += (-1) * sin(impAngle) * impVel;
		    ith.pos.x +=  (impRadi + tarRadi) * cos(impAngle) * 1.2 ; 
		    ith.pos.y +=  (impRadi + tarRadi) * sin(impAngle) * 1.2 ;

		    
                    if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
                }

                for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
                    ptcl.push_back(imp[i]);
                }

		for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		  if (imp[i].tag==0){
		      impNmntl += 1;
		  }else{
		    impNcore += 1;
		  }
                }

		impNptcl = impNmntl + impNcore;

		Nptcl = tarNptcl + impNptcl;
		
                break;
	      }
            case 2:
	      {
		///////////////////
		//Dummy put to determine # of ptcls
		///////////////////
		

		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		  for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
		    for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
		      const PS::F64 r = sqrt(x*x + y*y + z*z) * UnitRadi;
		      if(r >= tarRadi || r <= tarCoreRadi) continue;
		      ++ tarNmntl;
		    }
		  }
		}

		while(tarCoreShrinkFactor *= 0.99){
		  tarNcore = 0;
		  for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		    for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
		      for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
			const PS::F64 r = tarCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
			if(r >= Corr * tarCoreRadi) continue;
			++ tarNcore;
		      }
		    }
		  }
		  if((double)(tarNcore) / (double)(tarNcore + tarNmntl) > coreFracMass) break;
		}

		///////////////////
		//Dummy end
		///////////////////


		
		// checking if there are enough mantle particles
	      if (tarNmntl <  static_cast<int>(Nptcl * (1.0-coreFracMass))){
		std::cout << "Too few mantle particles. Increase the grid size. The easiest fix is to increase the gridpoint in GI.h " << std::endl;
		exit(0);
	      }

	      // checking if there are enough core particles
	      if (tarNcore <  static_cast<int>(Nptcl * coreFracMass)){
		std::cout << "Too few core particles. Increase the grid size and/or change the core shrink factor." << std::endl;
		exit(0);
	      }

	      //removing particles to reach the exact Nptcl

	      int index = 0;

	      std::unordered_set<unsigned int> removal_list;
	      removal_list = create_removal_list (0, tarNmntl, tarNmntl - static_cast<int>(Nptcl * (1.0-coreFracMass)));

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
                            ith.dens = (tarMass - tarCoreMass) / (4.0 / 3.0 * math::pi * (tarRadi * tarRadi * tarRadi - tarCoreRadi * tarCoreRadi * tarCoreRadi));
                            ith.mass = tarMass + impMass;
                            ith.eng  = 0.1 * Grav * tarMass / tarRadi;
                            ith.id   = id++;
                            ith.setPressure(&AGranite);
                            ith.tag = 0;

			    if (removal_list.count(index)){
			      id += -1;
			    }else{
			      if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
			    }
			    index += 1;
                        }
                    }
                }

		std::cout << "# of mantle particles = " <<  id  << std::endl;

		// making the core condition
		removal_list.clear();
		removal_list = create_removal_list (tarNmntl, tarNmntl + tarNcore, tarNcore - static_cast<int>(Nptcl * coreFracMass));
		
		index = tarNmntl;
	       
		
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
                            ith.mass = tarMass + impMass;
                            ith.eng  = 0.1 * Grav * tarMass / tarRadi;
                            ith.id   = id++;
                            ith.setPressure(&Iron);
                            ith.tag = 1;

			    if (removal_list.count(index)){
			      id += -1;
			    }else{
			      if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
			    }
			    index += 1;
			    
                        }
                    }
                }

		std::cout << "# of total particles = " <<  id  << std::endl;
		
		tarNmntl = static_cast<int>(Nptcl * (1.0 - coreFracMass));
		tarNcore = static_cast<int>(Nptcl * coreFracMass);
		
                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
                    tar[i].mass /= (PS::F64)(Nptcl);
                }
                for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
                    ptcl.push_back(tar[i]);
                }
                break;
	      }
        }

	


	
	        tarNptcl = tarNcore + tarNmntl;
		impNptcl = impNcore + impNmntl;

	        std::cout << "Target  :" << tarNptcl << std::endl;
		std::cout << "    radius           : " << tarRadi << std::endl;
		std::cout << "    total-to-core    : " << (double)(tarNcore) / (double)(tarNptcl) << std::endl;
		std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
		std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
		std::cout << "    core density     : " << tarCoreMass / (4.0 * math::pi / 3.0 * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr) << std::endl;
		std::cout << "    mantle density   : " << (tarMass - tarCoreMass) / (4.0 * math::pi / 3.0 * (tarRadi * tarRadi * tarRadi - tarCoreRadi * tarCoreRadi * tarCoreRadi)) << std::endl;
		std::cout << "    mean density     : " << tarMass / (4.0 * math::pi / 3.0 * tarRadi * tarRadi * tarRadi) << std::endl;

		if (mode==1){
		  std::cout << "Impactor:" << impNptcl << std::endl;
		  std::cout << "    radius           : " << impRadi << std::endl;
		  std::cout << "    total-to-core    : " << (double)(impNcore) / (double)(impNptcl) << std::endl;
		  std::cout << "    # of core ptcls  : " << impNcore << std::endl;
		  std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
		  std::cout << "    core density     : " << impCoreMass / (4.0 * math::pi / 3.0 * impCoreRadi * impCoreRadi * impCoreRadi * Corr * Corr * Corr) << std::endl;
		  std::cout << "    mantle density   : " << (impMass - impCoreMass) / (4.0 * math::pi / 3.0 * (impRadi * impRadi * impRadi - impCoreRadi * impCoreRadi * impCoreRadi)) << std::endl;
		  std::cout << "    mean density     : " << impMass / (4.0 * math::pi / 3.0 * impRadi * impRadi * impRadi) << std::endl;
		  std::cout << "    velocity         : " << impVel << std::endl;
		  std::cout << "Tar-to-Imp mass ratio: " << (double)(impNmntl) / (double)(tarNmntl) << std::endl;
		}

		assert(Nptcl == tarNptcl + impNptcl);
		
		std::cout << "Total number of particles:" << Nptcl << std::endl;
	
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
