#include <particle_simulator.hpp>
#include "EoS.h"
#include "param.h"
#include "mathfunc.h"
#include "kernel.h"
#include "class.h"
#include "kernel.cuh"
#include "class_device.hpp"

struct kernel_t_device{
	__device__ static real pow8(const real x){
		const real x2 = x  * x;
		const real x4 = x2 * x2;
		return x4 * x4;
	}
	__device__ static real plus(const real x){
		return (x > 0.0) ? x : 0.0;
	}
	__device__ static real pow7(const real x){
		real x2 = x * x;
		real x4 = x2 * x2;
		return x4 * x2 * x;
	}
	//W
	__device__ real W(const real r, const real h) const{
		const real H = supportRadius() * h;
		const real s = r / H;
		real r_value;
		//r_value = (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))) * pow8(plus(1.0 - s));
		r_value = __fmaf_rn(__fmaf_rn(__fmaf_rn(s, 32.0, 25.0), s, 8.0), s, 1.0) * pow8(plus(1.0 - s));
		r_value *= (1365./64.) / (H * H * H * M_PI);
		return r_value;
	}
	//gradW
	__device__ real gradW(const real r, const real h) const{
		const real H = supportRadius() * h;
		const real s = r / H;
		real r_value;
		r_value = pow7(plus(1.0 - s)) * (plus(1.0 - s) * (8.0 + s * (50.0 + s * (96.0))) - 8.0 * (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))));
		r_value *= (1365./64.) / (H * H * H * M_PI);
		return r_value / (H + 0.01 * h);
	}
	__device__ static real supportRadius(){
		return 3.5;
	}
};

static struct{
	int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	Drvt::EpiDev *epi_d, *epi_h;
	Drvt::EpjDev *epj_d, *epj_h;
	Drvt::ForceDev *res_d, *res_h;
}drvt_host;

static struct{
	int *ni_displc_d, *nj_displc_d, *ni_displc_h, *nj_displc_h;
	Hydr::EpiDev *epi_d, *epi_h;
	Hydr::EpjDev *epj_d, *epj_h;
	Hydr::ForceDev *res_d, *res_h;
}hydr_host;

__global__ void deviceCalcDensity(const Dens::EpiDev *epi, const int *ni_displc, const Dens::EpjDev *epj, const int *nj_displc, Dens::ForceDev *dens){
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	kernel_t_device kernel;
	const Dens::EpiDev& ith = epi[id];
	const int j_head = nj_displc[ith.id_walk];
	const int j_tail = nj_displc[ith.id_walk + 1];
	real dens_buf = 0.0;
	for(int j = j_head ; j < j_tail ; ++ j){
		const Dens::EpjDev& jth = epj[j];
		const double3 dr = make_double3(jth.rx - ith.rx, jth.ry - ith.ry, jth.rz - ith.rz);
		//const real r = sqrtf(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
		const real r = __fsqrt_rn(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
		dens_buf += jth.mass * kernel.W(r, ith.smth);
	}
	dens[id].dens = max(dens_buf, 5.0);
}

__global__ void deviceCalcDerivative(const Drvt::EpiDev *epi, const int *ni_displc, const Drvt::EpjDev *epj, const int *nj_displc, Drvt::ForceDev *force){
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	kernel_t_device kernel;
	const Drvt::EpiDev& ith = epi[id];
	const int j_head = nj_displc[ith.id_walk];
	const int j_tail = nj_displc[ith.id_walk + 1];
	double4 force_buf = make_double4(0.0, 0.0, 0.0, 0.0);
	for(int j = j_head ; j < j_tail ; ++ j){
		const Drvt::EpjDev& jth = epj[j];
		double3 dr = make_double3(jth.rx - ith.rx, jth.ry - ith.ry, jth.rz - ith.rz);
		double3 dv = make_double3(jth.vx - ith.vx, jth.vy - ith.vy, jth.vz - ith.vz);
		const real r = sqrtf(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + 1.0e-4;
		const real rinv = 1.0 / r;
		const real drdv = dr.x * dv.x + dr.y * dv.y + dr.z * dv.z;

		const real jth_mass_ith_abs_gradW = jth.mass * kernel.gradW(r, ith.smth) * rinv;
		force_buf.w += - drdv * jth_mass_ith_abs_gradW;
		force_buf.x += - (dr.y * dv.z - dr.z * dv.y) * jth_mass_ith_abs_gradW;
		force_buf.y += - (dr.z * dv.x - dr.x * dv.z) * jth_mass_ith_abs_gradW;
		force_buf.z += - (dr.x * dv.y - dr.y * dv.x) * jth_mass_ith_abs_gradW;
	}
	force[id].rot_vx = force_buf.x;
	force[id].rot_vy = force_buf.y;
	force[id].rot_vz = force_buf.z;
	force[id].div_v  = force_buf.w;
}

__global__ void deviceCalcHydroForce(const Hydr::EpiDev *epi, const int *ni_displc, const Hydr::EpjDev *epj, const int *nj_displc, Hydr::ForceDev *force){
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	kernel_t_device kernel;
	const Hydr::EpiDev& ith = epi[id];
	const int j_head = nj_displc[ith.id_walk];
	const int j_tail = nj_displc[ith.id_walk + 1];

	double v_sig_max = 0.0;
	double4 force_buf = make_double4(0.0, 0.0, 0.0, 0.0);
	const real ith_pres_over_dens2 = ith.pres / (ith.dens * ith.dens);
	for(int j = j_head ; j < j_tail ; ++ j){
		const Hydr::EpjDev& jth = epj[j];
		const double3 dr = make_double3(jth.rx - ith.rx, jth.ry - ith.ry, jth.rz - ith.rz);
		const double3 dv = make_double3(jth.vx - ith.vx, jth.vy - ith.vy, jth.vz - ith.vz);
		/*
		const real r  = sqrtf(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z + 1.0e-8f);
		const real rinv = __fdiv_rn(1.0f, r);
		*/
		const real r2 = __fmaf_rn(dr.x, dr.x, __fmaf_rn(dr.y, dr.y, __fmaf_rn(dr.z, dr.z, 1.0e-4)));
		const real rinv = __frsqrt_rn(r2);
		const real r = r2 * rinv;
		const real drdv = dr.x * dv.x + dr.y * dv.y + dr.z * dv.z;
		//AV
		const real w_ij = (drdv < 0.0) ? drdv * rinv : 0.0;
		const real v_sig = ith.snds + jth.snds - 3.0 * w_ij;
		v_sig_max = (v_sig_max < v_sig) ? v_sig : v_sig_max;
		const real AV = - 0.5 * v_sig * w_ij / (0.5 * (ith.dens + jth.dens));
		/*
		const real ith_abs_gradW = kernel.gradW(r, ith.smth);
		const real jth_abs_gradW = kernel.gradW(r, jth.smth);
		const real abs_gradW_rinv = 0.5f * (ith_abs_gradW + jth_abs_gradW) * rinv;
		*/
		const real abs_gradW_rinv = 0.5 * kernel.gradW(r, 0.5 * (ith.smth + jth.smth)) * rinv;
		const double3 gradW = make_double3(abs_gradW_rinv * dr.x, abs_gradW_rinv * dr.y, abs_gradW_rinv * dr.z);
		const real acc = jth.mass * (ith_pres_over_dens2 + jth.pres / (jth.dens * jth.dens) + AV);
		force_buf.x += acc * gradW.x;
		force_buf.y += acc * gradW.y;
		force_buf.z += acc * gradW.z;
		force_buf.w += jth.mass * (ith_pres_over_dens2 + 0.5 * AV) * (dv.x * gradW.x + dv.y * gradW.y + dv.z * gradW.z);
	}
	force[id].ax = force_buf.x;
	force[id].ay = force_buf.y;
	force[id].az = force_buf.z;
	force[id].eng_dot = force_buf.w;
	force[id].dt = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
}

namespace DENS{
	static struct{
		Dens::ForceDev *res;
		int *ni_displc, *nj_displc;
		Dens::EpiDev *epi;
		Dens::EpjDev *epj;
	}host, dev;

	int DispatchKernel(const PS::S32 tag, const int n_walk, const STD::EPI::Dens* epi[], const int* n_epi, const STD::EPJ::Dens* epj[], const int* n_epj){
		static bool isFirst = true;
		if(isFirst == true){
			std::cout << "Alloc Cuda Vars.." << std::endl;
			(cudaMalloc    ((void**)&dev.ni_displc, (N_WALK_LIMIT + 1) * sizeof(int)));
			(cudaMalloc    ((void**)&dev.nj_displc, (N_WALK_LIMIT + 1) * sizeof(int)));
			(cudaMallocHost((void**)&host.ni_displc, (N_WALK_LIMIT + 1) * sizeof(int)));
			(cudaMallocHost((void**)&host.nj_displc, (N_WALK_LIMIT + 1) * sizeof(int)));
			(cudaMalloc    ((void**)&dev.epi, NI_LIMIT * sizeof(Dens::EpiDev)));
			(cudaMalloc    ((void**)&dev.epj, NJ_LIMIT * sizeof(Dens::EpjDev)));
			(cudaMalloc    ((void**)&dev.res, NI_LIMIT * sizeof(Dens::ForceDev)));
			(cudaMallocHost((void**)&host.epi, NI_LIMIT * sizeof(Dens::EpiDev)));
			(cudaMallocHost((void**)&host.epj, NJ_LIMIT * sizeof(Dens::EpjDev)));
			(cudaMallocHost((void**)&host.res, NI_LIMIT * sizeof(Dens::ForceDev)));
			isFirst = false;
		}
		host.ni_displc[0] = host.nj_displc[0] = 0;
		for(std::size_t i = 0; i < n_walk ; ++ i){
			host.ni_displc[i+1] = host.ni_displc[i] + n_epi[i];
			host.nj_displc[i+1] = host.nj_displc[i] + n_epj[i];
		}
		const PS::S32 ni_total = host.ni_displc[n_walk];
		if(ni_total >= NI_LIMIT){
			std::cout << ni_total << " >= " << NI_LIMIT << std::endl;
			assert(ni_total < NI_LIMIT);
		}
		const int ni_total_reg = host.ni_displc[n_walk] + ((ni_total % N_THREAD_GPU != 1) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
		//make data for the device on the host
		int cnt = 0;
		int cnt_j = 0;
		for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
			for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
				host.epi[cnt].rx = epi[walk][i].pos.x;
				host.epi[cnt].ry = epi[walk][i].pos.y;
				host.epi[cnt].rz = epi[walk][i].pos.z;
				host.epi[cnt].mass = epi[walk][i].mass;
				host.epi[cnt].smth = epi[walk][i].smth;
				host.epi[cnt].id_walk = walk;
				++ cnt;
			}
			for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
				host.epj[cnt_j].rx = epj[walk][j].pos.x;
				host.epj[cnt_j].ry = epj[walk][j].pos.y;
				host.epj[cnt_j].rz = epj[walk][j].pos.z;
				host.epj[cnt_j].mass = epj[walk][j].mass;
				host.epj[cnt_j].smth = epj[walk][j].smth;
				++ cnt_j;
			}
		}

		(cudaMemcpy(dev.epi, host.epi, ni_total_reg * sizeof(Dens::EpiDev), cudaMemcpyHostToDevice));
		(cudaMemcpy(dev.epj, host.epj, cnt_j * sizeof(Dens::EpjDev), cudaMemcpyHostToDevice));
		(cudaMemcpy(dev.ni_displc, host.ni_displc, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
		(cudaMemcpy(dev.nj_displc, host.nj_displc, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

		const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
		dim3 size_grid(n_grid, 1, 1);
		dim3 size_thread(N_THREAD_GPU, 1, 1);
		deviceCalcDensity<<<size_grid, size_thread>>> (dev.epi, dev.ni_displc, dev.epj, dev.nj_displc, dev.res);
		return 0;
	}

	int RetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, STD::RESULT::Dens* force[]){
		int ni_tot = 0;
		for(int i = 0 ; i < n_walk ; ++ i){
			ni_tot += ni[i];
		}
		(cudaMemcpy(host.res, dev.res, ni_tot * sizeof(Dens::ForceDev), cudaMemcpyDeviceToHost));
		int cnt = 0;
		for(int walk = 0 ; walk < n_walk ; ++ walk){
			for(int i = 0 ; i < ni[walk] ; ++ i){
				force[walk][i].dens = max(host.res[cnt].dens, 5.0);
				++ cnt;
			}
		}
		return 0;
	}
};

int DrvtDispatchKernel(const PS::S32 tag, const int n_walk, const STD::EPI::Drvt* epi[], const int* n_epi, const STD::EPJ::Drvt* epj[], const int* n_epj){
	static bool isFirst = true;
	if(isFirst == true){
		std::cout << "Alloc Cuda Vars.." << std::endl;
		(cudaMalloc    ((void**)&drvt_host.ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&drvt_host.nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&drvt_host.ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&drvt_host.nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&drvt_host.epi_d, NI_LIMIT * sizeof(Drvt::EpiDev)));
		(cudaMalloc    ((void**)&drvt_host.epj_d, NJ_LIMIT * sizeof(Drvt::EpjDev)));
		(cudaMalloc    ((void**)&drvt_host.res_d, NI_LIMIT * sizeof(Drvt::ForceDev)));
		(cudaMallocHost((void**)&drvt_host.epi_h, NI_LIMIT * sizeof(Drvt::EpiDev)));
		(cudaMallocHost((void**)&drvt_host.epj_h, NJ_LIMIT * sizeof(Drvt::EpjDev)));
		(cudaMallocHost((void**)&drvt_host.res_h, NI_LIMIT * sizeof(Drvt::ForceDev)));
		isFirst = false;
	}
	drvt_host.ni_displc_h[0] = drvt_host.nj_displc_h[0] = 0;
	for(std::size_t i = 0; i < n_walk ; ++ i){
		drvt_host.ni_displc_h[i+1] = drvt_host.ni_displc_h[i] + n_epi[i];
		drvt_host.nj_displc_h[i+1] = drvt_host.nj_displc_h[i] + n_epj[i];
	}
	const PS::S32 ni_total = drvt_host.ni_displc_h[n_walk];
	if(ni_total >= NI_LIMIT){
		std::cout << ni_total << " >= " << NI_LIMIT << std::endl;
		assert(ni_total < NI_LIMIT);
	}
	if(drvt_host.nj_displc_h[n_walk] >= NJ_LIMIT){
		std::cout << drvt_host.nj_displc_h[n_walk] << " >= " << NJ_LIMIT << std::endl;
		assert(drvt_host.nj_displc_h[n_walk] < NJ_LIMIT);
	}

	const int ni_total_reg = drvt_host.ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 1) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
	//make data for device on host
	int cnt = 0;
	int cnt_j = 0;
	for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
		for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
			drvt_host.epi_h[cnt].rx = epi[walk][i].pos.x;
			drvt_host.epi_h[cnt].ry = epi[walk][i].pos.y;
			drvt_host.epi_h[cnt].rz = epi[walk][i].pos.z;
			drvt_host.epi_h[cnt].vx = epi[walk][i].vel.x;
			drvt_host.epi_h[cnt].vy = epi[walk][i].vel.y;
			drvt_host.epi_h[cnt].vz = epi[walk][i].vel.z;
			drvt_host.epi_h[cnt].dens = epi[walk][i].dens;
			drvt_host.epi_h[cnt].smth = epi[walk][i].smth;
			drvt_host.epi_h[cnt].id_walk = walk;
			++ cnt;
		}
		for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
			drvt_host.epj_h[cnt_j].rx = epj[walk][j].pos.x;
			drvt_host.epj_h[cnt_j].ry = epj[walk][j].pos.y;
			drvt_host.epj_h[cnt_j].rz = epj[walk][j].pos.z;
			drvt_host.epj_h[cnt_j].vx = epj[walk][j].vel.x;
			drvt_host.epj_h[cnt_j].vy = epj[walk][j].vel.y;
			drvt_host.epj_h[cnt_j].vz = epj[walk][j].vel.z;
			drvt_host.epj_h[cnt_j].mass = epj[walk][j].mass;
			drvt_host.epj_h[cnt_j].smth = epj[walk][j].smth;
			++ cnt_j;
		}
	}

	(cudaMemcpy(drvt_host.epi_d, drvt_host.epi_h, ni_total_reg * sizeof(Drvt::EpiDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(drvt_host.epj_d, drvt_host.epj_h, cnt_j * sizeof(Drvt::EpjDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(drvt_host.ni_displc_d, drvt_host.ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
	(cudaMemcpy(drvt_host.nj_displc_d, drvt_host.nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

	const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
	dim3 size_grid(n_grid, 1, 1);
	dim3 size_thread(N_THREAD_GPU, 1, 1);
	deviceCalcDerivative<<<size_grid, size_thread>>> (drvt_host.epi_d, drvt_host.ni_displc_d, drvt_host.epj_d, drvt_host.nj_displc_d, drvt_host.res_d);
	return 0;
}

int DrvtRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, STD::RESULT::Drvt* force[]){
	int ni_tot = 0;
	for(int i = 0 ; i < n_walk ; ++ i){
		ni_tot += ni[i];
	}
	(cudaMemcpy(drvt_host.res_h, drvt_host.res_d, ni_tot * sizeof(Drvt::ForceDev), cudaMemcpyDeviceToHost));
	int cnt = 0;
	for(int walk = 0 ; walk < n_walk ; ++ walk){
		for(int i = 0 ; i < ni[walk] ; ++ i){
			force[walk][i].div_v   = drvt_host.res_h[cnt].div_v;
			force[walk][i].rot_v.x = drvt_host.res_h[cnt].rot_vx;
			force[walk][i].rot_v.y = drvt_host.res_h[cnt].rot_vy;
			force[walk][i].rot_v.z = drvt_host.res_h[cnt].rot_vz;
			++ cnt;
		}
	}
	return 0;
}

int HydrDispatchKernel(const PS::S32 tag, const int n_walk, const STD::EPI::Hydro** epi, const int* n_epi, const STD::EPJ::Hydro** epj, const int* n_epj){
	static bool isFirst = true;
	if(isFirst == true){
		std::cout << "Alloc Cuda Vars.." << std::endl;
		(cudaMalloc    ((void**)&hydr_host.ni_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&hydr_host.nj_displc_d, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&hydr_host.ni_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMallocHost((void**)&hydr_host.nj_displc_h, (N_WALK_LIMIT + 1) * sizeof(int)));
		(cudaMalloc    ((void**)&hydr_host.epi_d, NI_LIMIT * sizeof(Hydr::EpiDev)));
		(cudaMalloc    ((void**)&hydr_host.epj_d, NJ_LIMIT * sizeof(Hydr::EpjDev)));
		(cudaMalloc    ((void**)&hydr_host.res_d, NI_LIMIT * sizeof(Hydr::ForceDev)));
		(cudaMallocHost((void**)&hydr_host.epi_h, NI_LIMIT * sizeof(Hydr::EpiDev)));
		(cudaMallocHost((void**)&hydr_host.epj_h, NJ_LIMIT * sizeof(Hydr::EpjDev)));
		(cudaMallocHost((void**)&hydr_host.res_h, NI_LIMIT * sizeof(Hydr::ForceDev)));
		isFirst = false;
	}
	hydr_host.ni_displc_h[0] = hydr_host.nj_displc_h[0] = 0;
	for(std::size_t i = 0; i < n_walk ; ++ i){
		hydr_host.ni_displc_h[i+1] = hydr_host.ni_displc_h[i] + n_epi[i];
		hydr_host.nj_displc_h[i+1] = hydr_host.nj_displc_h[i] + n_epj[i];
	}
	const PS::S32 ni_total = hydr_host.ni_displc_h[n_walk];
	const int ni_total_reg = hydr_host.ni_displc_h[n_walk] + ((ni_total % N_THREAD_GPU != 0) ? (N_THREAD_GPU - (ni_total % N_THREAD_GPU)) : 0);
	//make data for device on host
	int cnt = 0;
	int cnt_j = 0;
	for(std::size_t walk = 0 ; walk < n_walk ; ++ walk){
		for(std::size_t i = 0 ; i < n_epi[walk] ; ++ i){
			hydr_host.epi_h[cnt].rx      = epi[walk][i].pos.x;
			hydr_host.epi_h[cnt].ry      = epi[walk][i].pos.y;
			hydr_host.epi_h[cnt].rz      = epi[walk][i].pos.z;
			hydr_host.epi_h[cnt].vx      = epi[walk][i].vel.x;
			hydr_host.epi_h[cnt].vy      = epi[walk][i].vel.y;
			hydr_host.epi_h[cnt].vz      = epi[walk][i].vel.z;
			hydr_host.epi_h[cnt].dens    = epi[walk][i].dens;
			hydr_host.epi_h[cnt].pres    = epi[walk][i].pres;
			hydr_host.epi_h[cnt].snds    = epi[walk][i].snds;
			hydr_host.epi_h[cnt].smth    = epi[walk][i].smth;
			hydr_host.epi_h[cnt].Bal     = epi[walk][i].Bal;
			hydr_host.epi_h[cnt].id_walk = walk;
			hydr_host.epi_h[cnt].grad_smth = epi[walk][i].grad_smth;
			++ cnt;
		}
		for(std::size_t j = 0 ; j < n_epj[walk] ; ++ j){
			hydr_host.epj_h[cnt_j].rx   = epj[walk][j].pos.x;
			hydr_host.epj_h[cnt_j].ry   = epj[walk][j].pos.y;
			hydr_host.epj_h[cnt_j].rz   = epj[walk][j].pos.z;
			hydr_host.epj_h[cnt_j].vx   = epj[walk][j].vel.x;
			hydr_host.epj_h[cnt_j].vy   = epj[walk][j].vel.y;
			hydr_host.epj_h[cnt_j].vz   = epj[walk][j].vel.z;
			hydr_host.epj_h[cnt_j].dens = epj[walk][j].dens;
			hydr_host.epj_h[cnt_j].pres = epj[walk][j].pres;
			hydr_host.epj_h[cnt_j].snds = epj[walk][j].snds;
			hydr_host.epj_h[cnt_j].mass = epj[walk][j].mass;
			hydr_host.epj_h[cnt_j].smth = epj[walk][j].smth;
			hydr_host.epj_h[cnt_j].Bal  = epj[walk][j].Bal;
			hydr_host.epj_h[cnt_j].grad_smth = epj[walk][j].grad_smth;
			++ cnt_j;
		}
	}

	(cudaMemcpy(hydr_host.epi_d, hydr_host.epi_h, ni_total_reg * sizeof(Hydr::EpiDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(hydr_host.epj_d, hydr_host.epj_h, cnt_j * sizeof(Hydr::EpjDev), cudaMemcpyHostToDevice));
	(cudaMemcpy(hydr_host.ni_displc_d, hydr_host.ni_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));
	(cudaMemcpy(hydr_host.nj_displc_d, hydr_host.nj_displc_h, (n_walk + 1) * sizeof(int), cudaMemcpyHostToDevice));

	const int n_grid = ni_total_reg / N_THREAD_GPU + ((ni_total_reg % N_THREAD_GPU == 0) ? 0 : 1);
	dim3 size_grid(n_grid, 1, 1);
	dim3 size_thread(N_THREAD_GPU, 1, 1);
	deviceCalcHydroForce<<<size_grid, size_thread>>> (hydr_host.epi_d, hydr_host.ni_displc_d, hydr_host.epj_d, hydr_host.nj_displc_d, hydr_host.res_d);
	return 0;
}

int HydrRetrieveKernel(const PS::S32 tag, const PS::S32 n_walk, const PS::S32* ni, STD::RESULT::Hydro** force){
	int ni_tot = 0;
	for(int i = 0 ; i < n_walk ; ++ i){
		ni_tot += ni[i];
	}
	(cudaMemcpy(hydr_host.res_h, hydr_host.res_d, ni_tot * sizeof(Hydr::ForceDev), cudaMemcpyDeviceToHost));
	int cnt = 0;
	for(int walk = 0 ; walk < n_walk ; ++ walk){
		for(int i = 0 ; i < ni[walk] ; ++ i){
			force[walk][i].acc.x = hydr_host.res_h[cnt].ax;
			force[walk][i].acc.y = hydr_host.res_h[cnt].ay;
			force[walk][i].acc.z = hydr_host.res_h[cnt].az;
			force[walk][i].eng_dot = hydr_host.res_h[cnt].eng_dot;
			force[walk][i].dt = hydr_host.res_h[cnt].dt;
			++ cnt;
		}
	}
	return 0;
}


