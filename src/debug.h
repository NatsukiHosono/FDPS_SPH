//
// Created by Scott Hull on 8/20/20.
//
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <experimental/filesystem>


#ifndef FDPS_SPH_DEBUG_H
#define FDPS_SPH_DEBUG_H

#endif //FDPS_SPH_DEBUG_H

class EnergyOutFile {
public:
    std::ofstream EnergyDebugFile;

    void write_timestep(PS::ParticleSystem<STD::RealPtcl> &sph_system, system_t &sysinfo) {
        EnergyDebugFile.open("energy_debug/energy_debug" + std::to_string(sysinfo.output_id) + ".txt");
#pragma omp parallel for
        for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            if (sph_system[i].eng < 0) {
                std::string p_id = std::to_string(sph_system[i].id);
                std::string p_tag = std::to_string(sph_system[i].tag);
                std::string p_eng = std::to_string(sph_system[i].eng);
                EnergyDebugFile << p_id << "," << p_tag << "," << p_eng << "\n";
            }
        }
        EnergyDebugFile.close();
    };

    EnergyOutFile() {
//        std::experimental::filesystem::remove_all("energy_debug");
//        std::experimental::filesystem::create_directory("energy_debug");
    };
};
