#pragma once

template <class ThisPtcl> void OutputBinary(PS::ParticleSystem<ThisPtcl>& sph_system, const system_t& sysinfo, char* out_dir){
	//Binary
	char filename[256];
	std::ofstream fout;
    sprintf(filename, "results/%s/%05d_%05d.bin", out_dir, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
	fout.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	if(!fout){
		std::cout << "cannot write restart file." << std::endl;
		exit(1);
	}
	fout.write(reinterpret_cast<const char * const>(&sysinfo), sizeof(system_t));
	for(std::size_t i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		ThisPtcl ith = sph_system[i];
		fout.write((char*)&ith, sizeof(ThisPtcl));
	}
	fout.close();
}

template <class ThisPtcl> void OutputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& sph_system, system_t& sysinfo, const PS::F64 end_time, char* out_dir){
	if(sysinfo.time > sysinfo.output_time){
		FileHeader header;
		header.time = sysinfo.time;
		header.Nbody = sph_system.getNumberOfParticleLocal();
		char filename[256];
        sprintf(filename, "results/%s/%05d", out_dir, sysinfo.output_id);
		sph_system.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
		if(PS::Comm::getRank() == 0){
            std::cout << "//================================" << std::endl;
            std::cout << "output " << filename << "." << std::endl;
			std::cout << "//================================" << std::endl;
		}
		sysinfo.output_time += end_time / PARAM::NUMBER_OF_SNAPSHOTS;
		++ sysinfo.output_id;
	}
}

template <class ThisPtcl> void InputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& sph_system, system_t& sysinfo){
	FileHeader header;
	char filename[256];
    sprintf(filename, "results/%05d", sysinfo.step);
	sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
	sysinfo.time = header.time;
	std::cout << header.time << std::endl;
}

template <class ThisPtcl> void InputBinary(PS::ParticleSystem<ThisPtcl>& sph_system, system_t& sysinfo){
	char filename[256];
	//sprintf(filename, "result/%05d_%05d_%05d.bin", sysinfo->step, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
	sprintf(filename, "./%05d_%05d.bin", PS::Comm::getNumberOfProc(), PS::Comm::getRank());
	std::ifstream fin(filename, std::ios::in | std::ios::binary);
	if(!fin){
		std::cout << "cannot open restart file." << std::endl;
		exit(1);
	}
	std::vector<ThisPtcl> ptcl;
	fin.read((char*)&sysinfo, sizeof(system_t));
	//sysinfo.step ++;
	while(1){
		ThisPtcl ith;
		fin.read((char*)&ith, sizeof(ThisPtcl));
		if(fin.eof() == true) break;
		ptcl.push_back(ith);
	}
	fin.close();
	sph_system.setNumberOfParticleLocal(ptcl.size());
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
		sph_system[i] = ptcl[i];
	}
	std::cout << "read binary" << std::endl;
}

