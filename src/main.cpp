#include <particle_simulator.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>

#include "param.h"
#include "mathfunc.h"
#include "kernel.h"
#include "EoS.h"
#include "class.h"
#include "init/GI.h"
#include "force.h"
#include "io.h"
#include "integral.h"

int main(int argc, char* argv[]){
	namespace PTCL = STD;

	// There should be a capitalization rule for names, you currently use ALLCAPS
	// for namespaces, but here you also use it for GI which is a class type. I would
	// recommend using 'CamelCase' for everything except variable names, which you already
	// name as 'names_with_underscores', but that is a personal preference and you are welcome
	// to use your own.

	// Also there needs to be a way to include tests into the program. For example you could add a
	// problem description called 'test' instead of 'GI' and if the program is started with a '--test'
	// flag this test setup is executed. The results could then be compared to a reference test result.

	typedef GI<PTCL::RealPtcl> PROBLEM;
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
	if(argc == 1){
		PROBLEM::setupIC(sph_system, sysinfo, dinfo);
		PROBLEM::setEoS(sph_system);
		PTCL::CalcPressure(sph_system);
	}else{
		sysinfo.step = atoi(argv[1]);
		InputFileWithTimeInterval<PTCL::RealPtcl>(sph_system, sysinfo);
		PROBLEM::setEoS(sph_system);
	}

	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].initialize();
	}
	OutputFileWithTimeInterval(sph_system, sysinfo, PROBLEM::END_TIME);

	//Dom. info
	dinfo.decomposeDomainAll(sph_system);
	sph_system.exchangeParticle(dinfo);
	//plant tree
	PS::TreeForForceShort<PTCL::RESULT::Dens , PTCL::EPI::Dens , PTCL::EPJ::Dens >::Gather   dens_tree;
	PS::TreeForForceShort<PTCL::RESULT::Drvt , PTCL::EPI::Drvt , PTCL::EPJ::Drvt >::Gather   drvt_tree;
	PS::TreeForForceShort<PTCL::RESULT::Hydro, PTCL::EPI::Hydro, PTCL::EPJ::Hydro>::Symmetry hydr_tree;

	// It would be advantageous if you could settle on a single style for input parameters, currently
	// the code mixes variables, templates, and preprocessor macros. Ultimately it would be nice to
	// have text files outside the source code that specify the model, but settling on a single header
	// file would also be fine for the moment. If you want to save memory by not allocating a tree if this
	// parameter is not set, consider making it an owning pointer (std::unique_ptr) like the following:

	std::unique_ptr<PS::TreeForForceLong <PTCL::RESULT::Grav , PTCL::EPI::Grav , PTCL::EPJ::Grav >::Monopole> grav_tree;
	if (PROBLEM::self_gravity)
	  grav_tree.reset(new PS::TreeForForceLong <PTCL::RESULT::Grav , PTCL::EPI::Grav , PTCL::EPJ::Grav >::Monopole());

	// This way only the memory for a single pointer is allocated if the tree is not needed, and you can
	// check for that by using if (grav_tree)

	dens_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, 128);
	drvt_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, 128);
	hydr_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, 128);
        if (PROBLEM::self_gravity)
                grav_tree->initialize(sph_system.getNumberOfParticleLocal(), 0.5, 8, 256);

        for(short int loop = 0 ; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP ; ++ loop){
		dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo);
	}

	PTCL::CalcPressure(sph_system);
	drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo);
	hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo);
        if (PROBLEM::self_gravity)
          grav_tree->calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(), PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

        sysinfo.dt = getTimeStepGlobal<PTCL::RealPtcl>(sph_system);
	PROBLEM::addExternalForce(sph_system, sysinfo);
	OutputFileWithTimeInterval(sph_system, sysinfo, PROBLEM::END_TIME);

	if(PS::Comm::getRank() == 0){
		std::cout << "//================================" << std::endl;
		std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
		std::cout << "step = " << sysinfo.step << std::endl;
		std::cout << "//================================" << std::endl;
	}

	while(sysinfo.time < PROBLEM::END_TIME){
		#pragma omp parallel for
		for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].initialKick(sysinfo.dt);
			sph_system[i].fullDrift(sysinfo.dt);
			sph_system[i].predict(sysinfo.dt);
		}
		sysinfo.time += sysinfo.dt;
		sph_system.adjustPositionIntoRootDomain(dinfo);
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		PROBLEM::setEoS(sph_system);

		for(short int loop = 0 ; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP ; ++ loop){
			dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo);
		}
		PTCL::CalcPressure(sph_system);
		drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo);
		hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo);
	        if (PROBLEM::self_gravity)
	          grav_tree->calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(), PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

	        PROBLEM::addExternalForce(sph_system, sysinfo);
		#pragma omp parallel for
		for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].finalKick(sysinfo.dt);
		}
		PROBLEM::postTimestepProcess(sph_system, sysinfo);
		sysinfo.dt = getTimeStepGlobal<PTCL::RealPtcl>(sph_system);
		OutputFileWithTimeInterval<PTCL::RealPtcl>(sph_system, sysinfo, PROBLEM::END_TIME);
		++ sysinfo.step;
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
			std::cout << "step = " << sysinfo.step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		if(sysinfo.step % 30 == 0){
			OutputBinary(sph_system, sysinfo);
		}
	}

	PS::Finalize();
	return 0;
}

