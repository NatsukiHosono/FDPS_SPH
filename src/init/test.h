
/*
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	/////////
	const PS::F64 Mass = 1.0;
	const PS::F64 Radi = 1.0;
	const PS::F64 Grav = 1.0;
	*end_time = 0.77;
	const PS::F64 dx = 1.0 / 8.0;
	PS::S32 id = 0;
	for(PS::F64 x = -Radi ; x <= Radi ; x += dx){
		for(PS::F64 y = -Radi ; y <= Radi ; y += dx){
			for(PS::F64 z = -Radi ; z <= Radi ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r > Radi) continue;
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.pos *= sqrt(r);

				ith.dens = (r <= 1.0e-10) ? 1.0 / dx : Mass / (2.0 * math::pi * Radi * Radi * r);
				ith.mass = Mass;
				ith.eng  = 0.05 * Grav * Mass / Radi;
				ith.id   = id++;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = ptcl[i].mass / (PS::F64)(ptcl.size());
	}

	if(PS::Comm::getRank() == 0){
		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			sph_system[i] = ptcl[i];
		}
	}else{
		sph_system.setNumberOfParticleLocal(0);
	}
	//Fin.
	std::cout << "setup..." << std::endl;
}
*/

template <class Ptcl> class GI : public Problem<Ptcl>{
	public:
	static const double END_TIME;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){

		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		/////////
		const PS::F64 Mass = 1.0;
		const PS::F64 Radi = 1.0;
		const PS::F64 Grav = 1.0;
		const PS::F64 dx = 1.0 / 8.0;
		PS::S32 id = 0;
		for(PS::F64 x = -Radi ; x <= Radi ; x += dx){
			for(PS::F64 y = -Radi ; y <= Radi ; y += dx){
				for(PS::F64 z = -Radi ; z <= Radi ; z += dx){
					const PS::F64 r = sqrt(x*x + y*y + z*z);
					if(r > Radi) continue;
					Ptcl ith;
					ith.pos.x = x;
					ith.pos.y = y;
					ith.pos.z = z;
					ith.pos *= sqrt(r);

					ith.dens = (r <= 1.0e-10) ? 1.0 / dx : Mass / (2.0 * math::pi * Radi * Radi * r);
					ith.mass = Mass;
					ith.eng  = 0.05 * Grav * Mass / Radi;
					ith.id   = id++;
					ptcl.push_back(ith);
				}
			}
		}
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].mass = ptcl[i].mass / (PS::F64)(ptcl.size());
		}

		if(PS::Comm::getRank() == 0){
			const PS::S32 numPtclLocal = ptcl.size();
			sph_system.setNumberOfParticleLocal(numPtclLocal);
			for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
				sph_system[i] = ptcl[i];
			}
		}else{
			sph_system.setNumberOfParticleLocal(0);
		}
		//Fin.
		std::cout << "setup..." << std::endl;

	}

	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure(&Monoatomic);
		}
	}

	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo){
	}
};

template <class Ptcl>
const double GI<Ptcl>::END_TIME = 0.77;
