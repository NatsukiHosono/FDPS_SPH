#pragma once
#include "class_device.hpp"

//const int N_THREAD_GPU = 2688;
const int N_THREAD_GPU = 1024;
const int N_WALK_LIMIT = 1000;
const int NI_LIMIT     = 1000 * N_WALK_LIMIT;
const int NJ_LIMIT     = 10000 * N_WALK_LIMIT;

PS::S32 DrvtDispatchKernel(const PS::S32, const PS::S32, const STD::EPI::Drvt** , const PS::S32*, const STD::EPJ::Drvt** , const PS::S32*);
PS::S32 HydrDispatchKernel(const PS::S32, const PS::S32, const STD::EPI::Hydro**, const PS::S32*, const STD::EPJ::Hydro**, const PS::S32*);
PS::S32 GravDispatchKernel(const PS::S32, const PS::S32, const STD::EPI::Grav** , const PS::S32*, const STD::EPJ::Grav** , const PS::S32*, const PS::SPJMonopole**, const PS::S32*);

PS::S32 DrvtRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, STD::RESULT::Drvt** );
PS::S32 HydrRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, STD::RESULT::Hydro**);
PS::S32 GravRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, STD::RESULT::Grav** );

namespace DENS{
	PS::S32 DispatchKernel(const PS::S32, const PS::S32, const STD::EPI::Dens**, const PS::S32*, const STD::EPJ::Dens**, const PS::S32*);
	PS::S32 RetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, STD::RESULT::Dens**);
	PS::S32 DispatchKernelDynamicAllocate(const PS::S32, const PS::S32, const STD::EPI::Dens**, const PS::S32*, const STD::EPJ::Dens**, const PS::S32*);
	PS::S32 RetrieveKernelDynamicAllocate(const PS::S32, const PS::S32, const PS::S32*, STD::RESULT::Dens**);
};

