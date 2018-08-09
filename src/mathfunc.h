#pragma once

namespace math{
        // You probably created these functions, because you are concerned about speed in
        // performance critical parts of your code? In this case you should replace the call
        // to atan by just stating pi = 3.14159265358979323846. atan is actually an expensive
        // function to call
	const PS::F64 pi = atan(1.0) * 4.0;

	// There is std::numeric_limits<PS::F64>::signaling_NaN() for this purpose, or std::numeric_limits<PS::F64>::quiet_NaN()
	// if you do not want to throw an exception when this value is accessed
	const PS::F64 NaN = std::numeric_limits<PS::F64>::signaling_NaN();

	// This is available as std::numeric_limits<PS::F64>::max()
	const PS::F64 VERY_LARGE_VALUE = std::numeric_limits<PS::F64>::max();
	template <typename type> inline type plus(const type arg){
		return (arg > 0) ? arg : 0;
	}
	template <typename type> inline type sign(const type arg){
		return (arg > 0) ? 1.0 : - 1.0;
	}
	template <typename type> inline type pow2(const type arg){
		return arg * arg;
	}
	template <typename type> inline type pow3(const type arg){
		return arg * arg * arg;
	}
	template <typename type> inline type pow4(const type arg){
		const type arg2 = arg * arg;
		return arg2 * arg2;
	}
	template <typename type> inline type pow5(const type arg){
		const type arg2 = arg * arg;
		return arg2 * arg2 * arg;
	}
	template <typename type> inline type pow6(const type arg){
		const type arg3 = arg * arg * arg;
		return arg3 * arg3;
	}
	template <typename type> inline type pow7(const type arg){
		const type arg3 = arg * arg * arg;
		return arg3 * arg3 * arg;
	}
	template <typename type> inline type pow8(const type arg){
		const type arg2 = arg * arg;
		const type arg4 = arg2 * arg2;
		return arg4 * arg4;
	}
}

