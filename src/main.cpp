#include <particle_simulator.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "param.h"
#include "mathfunc.h"
#include "kernel.h"
#include "EoS.h"
#include "class.h"
#include "GI.h"
#include "force.h"
#include "io.h"
#include "integral.h"


template<class Ptcl> double GI_universal<Ptcl>::end_time;
template<class Ptcl> double GI_universal<Ptcl>::damping;


int main(int argc, char *argv[]) {
    namespace PTCL = STD;
    typedef GI_universal<PTCL::RealPtcl> PROBLEM;  // note: this must be changed to GI_imp for impact simulations...come up with a fix soon for mode dependence
    //////////////////
    //Create vars.
    //////////////////
    PS::Initialize(argc, argv);
    PS::ParticleSystem<PTCL::RealPtcl> sph_system;
    sph_system.initialize();
    PS::DomainInfo dinfo;
    dinfo.initialize();
    system_t sysinfo;

    sph_system.setAverageTargetNumberOfSampleParticlePerProcess(200);
    //////////////////
    //Setup Initial
    //////////////////
    bool newSim = true;
    std::string input_file("input.txt");

    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input_file = std::string(argv[i + 1]);
        }
        if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--resume") == 0) {
            sysinfo.step = atoi(argv[i + 1]);
            newSim = false;
        }
    }
    ParameterFile parameter_file(input_file);
    std::cout << "Reading parameters from " << input_file << std::endl;
    std::string output_directory = parameter_file.getValueOf("output_directory", std::string("results/"));
    if (output_directory.back() != '/')
        output_directory.back() += '/';
    createOutputDirectory(output_directory);
    if (newSim) {
        PROBLEM::setupIC(sph_system, sysinfo, dinfo, parameter_file);
        PROBLEM::setEoS(sph_system);
        PTCL::CalcPressure(sph_system);
    } else {
        InputFileWithTimeInterval<PTCL::RealPtcl>(sph_system, sysinfo);
        PROBLEM::setEoS(sph_system);
    }
    PS::F64 output_interval = parameter_file.getValueOf("output_interval", 50.);

#pragma omp parallel for
    for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
        sph_system[i].initialize();
    }
    OutputFileWithTimeInterval(sph_system, sysinfo, output_interval, output_directory);

    //Dom. info
    dinfo.decomposeDomainAll(sph_system);
    sph_system.exchangeParticle(dinfo);
    //plant tree
    PS::TreeForForceShort<PTCL::RESULT::Dens, PTCL::EPI::Dens, PTCL::EPJ::Dens>::Gather dens_tree;
    PS::TreeForForceShort<PTCL::RESULT::Drvt, PTCL::EPI::Drvt, PTCL::EPJ::Drvt>::Gather drvt_tree;
    PS::TreeForForceShort<PTCL::RESULT::Hydro, PTCL::EPI::Hydro, PTCL::EPJ::Hydro>::Symmetry hydr_tree;
#ifdef SELF_GRAVITY
    PS::TreeForForceLong<PTCL::RESULT::Grav, PTCL::EPI::Grav, PTCL::EPJ::Grav>::Monopole grav_tree;
#endif

    dens_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, 128);
    drvt_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, 128);
    hydr_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, 128);
#ifdef SELF_GRAVITY
    grav_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 8, 256);
#endif
    for (short int loop = 0; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP; ++loop) {
        dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo);
    }

    PTCL::CalcPressure(sph_system);
    drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo);
    hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo);
#ifdef SELF_GRAVITY
    grav_tree.calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(),
                                       PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
#endif
    sysinfo.dt = getTimeStepGlobal<PTCL::RealPtcl>(sph_system);
    PROBLEM::addExternalForce(sph_system, sysinfo);
    OutputFileWithTimeInterval(sph_system, sysinfo, output_interval, output_directory);

    if (PS::Comm::getRank() == 0) {
        std::cout << "//================================" << std::endl;
        std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt
                  << std::endl;
        std::cout << "step = " << sysinfo.step << std::endl;
        std::cout << "//================================" << std::endl;
    }

    const unsigned int mode = parameter_file.getValueOf("mode", 1); // get modelling mode from input file
    const double initial_mantle_entropy = parameter_file.getValueOf("mantle_entropy",
                                                             3.12432609E+03); // initial constant entropy value
    const double initial_core_entropy = parameter_file.getValueOf("core_entropy",
                                                             3.12432609E+03); // initial constant entropy value
    const unsigned int aneos_grid_size = parameter_file.getValueOf("aneos_grid_size",
                                                                   120); // get grid size of ANEOS input file
    const unsigned int tillotson_grid_size = parameter_file.getValueOf("tillotson_grid_size",
                                                                       120); // get grid size of Tillotson input file

    while (sysinfo.time < PROBLEM::end_time) {
#pragma omp parallel for
        for (int i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            sph_system[i].initialKick(sysinfo.dt);
            sph_system[i].fullDrift(sysinfo.dt);
            sph_system[i].predict(sysinfo.dt);
        }
        sysinfo.time += sysinfo.dt;
        sph_system.adjustPositionIntoRootDomain(dinfo);
        dinfo.decomposeDomainAll(sph_system);
        sph_system.exchangeParticle(dinfo);
        PROBLEM::setEoS(sph_system);

        for (short int loop = 0; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP; ++loop) {
            dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo);
        }
        PTCL::CalcPressure(sph_system);
        drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo);
        hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo);
#ifdef SELF_GRAVITY
        grav_tree.calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(),
                                           PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
#endif
        PROBLEM::addExternalForce(sph_system, sysinfo);
#pragma omp parallel for
        for (int i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
            sph_system[i].finalKick(sysinfo.dt);
            sph_system[i].dampMotion(PROBLEM::damping);
        }
        //    Calculate initial internal energy for mode 1 initial target/impactor creation
        if (mode == 2) {
            PTCL::SetConstantEntropy(sph_system, initial_mantle_entropy, initial_core_entropy);
            PTCL::CalcAll(sph_system, aneos_grid_size, tillotson_grid_size);
        }
        PROBLEM::postTimestepProcess(sph_system, sysinfo);
        sysinfo.dt = getTimeStepGlobal<PTCL::RealPtcl>(sph_system);
        OutputFileWithTimeInterval<PTCL::RealPtcl>(sph_system, sysinfo, output_interval, output_directory);
        ++sysinfo.step;
        if (PS::Comm::getRank() == 0) {
            std::cout << "//================================" << std::endl;
            std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = "
                      << sysinfo.dt << std::endl;
            std::cout << "step = " << sysinfo.step << std::endl;
            std::cout << "//================================" << std::endl;
        }
        if (sysinfo.step % 30 == 0) {
            OutputBinary(sph_system, sysinfo, output_directory);
        }
    }

    PS::Finalize();
    return 0;
}

