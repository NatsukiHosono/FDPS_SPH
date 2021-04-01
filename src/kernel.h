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
	PS::F64 intW(const PS::F64vec dr, const PS::F64 h) const {
        const PS::F64 H = supportRadius() * h;
        const PS::F64 s = sqrt(dr * dr) / H;
        PS::F64 r_value;
        r_value = ((65.0 / 12.0) * math::pow12(s)) - ((519.0 / 11.0) * math::pow11(s)) +
                ((906.0 / 5.0) * math::pow10(s)) - ((1204 / 3.0) * math::pow9(s)) + ((2247.0 / 4.0) * math::pow8(s)) -
                (510.0 * math::pow7(s)) + (294.0 * math::pow6(s)) - ((492.0 / 5.0) * math::pow5(s)) +
                ((57.0 / 4.0) * math::pow4(s)) + (math::pow3(s) / 3);
        r_value *= (1365.0 / 16.0);
        return r_value;
	}
};

typedef WendlandC6 kernel_t;

