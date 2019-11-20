#pragma once

#include <sys/stat.h>
#include <dirent.h>
#include <cassert>

template <class ThisPtcl> void OutputBinary(PS::ParticleSystem<ThisPtcl>& sph_system, const system_t& sysinfo, std::string &out_dir){
	//Binary
	char filename[256];
	std::ofstream fout;
    sprintf(filename, "%s/%05d_%05d.bin", out_dir.c_str(), PS::Comm::getNumberOfProc(), PS::Comm::getRank());
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

template <class ThisPtcl> void OutputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& sph_system, system_t& sysinfo, const PS::F64 output_interval, std::string &out_dir){
	if(sysinfo.time > sysinfo.output_time){
		FileHeader header;
		header.time = sysinfo.time;
		header.Nbody = sph_system.getNumberOfParticleLocal();
		char filename[20];
        sprintf(filename, "results.%05lld", sysinfo.output_id);
		std::string full_filename = out_dir + filename;

		sph_system.writeParticleAscii(full_filename.c_str(), "%s_%05d_%05d.dat", header);
		if(PS::Comm::getRank() == 0){
            std::cout << "//================================" << std::endl;
            std::cout << "output " << full_filename << "." << std::endl;
			std::cout << "//================================" << std::endl;
		}
		sysinfo.output_time += output_interval;
		++ sysinfo.output_id;
	}
}

template <class ThisPtcl> void InputFileWithTimeInterval(PS::ParticleSystem<ThisPtcl>& sph_system, system_t& sysinfo){
	FileHeader header;
	char filename[256];
    sprintf(filename, "results/%05lld", sysinfo.step);
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

void createOutputDirectory(const std::string &directory_name){
    // check if the output directory exists, if not create it.
	if (DIR *output_directory = opendir(directory_name.c_str())){
		int error = closedir(output_directory);
		assert(error == 0);
	}
	else{
		int error = mkdir(directory_name.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
		assert(error == 0);
	}
}


