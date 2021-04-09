#pragma once

//Wendland C6
struct WendlandC6{
	WendlandC6(){}
	//W
	PS::F64 W(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))) * math::pow8(math::plus(1.0 - s));
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		r_value *= (78./7.) / (H * H * math::pi);
		#else
		r_value *= (1365./64.) / (H * H * H * math::pi);
        #endif
		return r_value;
	}
	//gradW
	PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = math::pow7(math::plus(1.0 - s)) * (math::plus(1.0 - s) * (8.0 + s * (50.0 + s * (96.0))) - 8.0 * (1.0 + s * (8.0 + s * (25.0 + s * (32.0)))));
		#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		r_value *= (78./7.) / (H * H * math::pi);
		#else
		r_value *= (1365./64.) / (H * H * H * math::pi);
		#endif
		return dr * r_value / (sqrt(dr * dr) * H  + 1.0e-6 * h);
	}
	static PS::F64 supportRadius(){
		return 2.5;
	}
	PS::F64 intWs2(const PS::F64vec dr, const PS::F64 h) const {
        // int_0^s_ij W(s) s^2 ds
        const PS::F64 H = supportRadius() * h;
        const PS::F64 s = sqrt(dr * dr) / H;
        const PS::F64 coeff = (1.0 / (64.0 * math::pi * math::pow3(H)));
        const PS::F64 var1 = ((624.0 * math::pow2(s)) - (4851.0 * s) + 16016.0);
        const PS::F64 var2 = s * var1 - 28665.0;
        const PS::F64 var3 = (5.0 * s) * var2 + 144144.0;
        const PS::F64 var4 = s * var3 - 70070.0;
        const PS::F64 var5 = var4 * math::pow6(s);
        PS::F64 r_value;
        r_value = math::pow3(s) * ((12870.0 * math::pow4(s)) - (3003.0 * math::pow2(s)) + var5 + 455.0);
        #ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        r_value *= (78./7.) / (H * H * math::pi);
        #else
        r_value *= coeff;
        #endif
        return r_value;
	}
};

typedef WendlandC6 kernel_t;

