// Avisynth v2.5.  Copyright 2002 Ben Rudiak-Gould et al.
// http://www.avisynth.org

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA, or visit
// http://www.gnu.org/copyleft/gpl.html .
//
// Linking Avisynth statically or dynamically with other modules is making a
// combined work based on Avisynth.  Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Avisynth give you
// permission to link Avisynth with independent modules that communicate with
// Avisynth solely through the interfaces defined in avisynth.h, regardless of the license
// terms of these independent modules, and to copy and distribute the
// resulting combined work under terms of your choice, provided that
// every copy of the combined work is accompanied by a complete copy of
// the source code of Avisynth (the version of Avisynth used to produce the
// combined work), being distributed under the terms of the GNU General
// Public License plus this exception.  An independent module is a module
// which is not derived from or based on Avisynth, such as 3rd-party filters,
// import and export plugins, or graphical user interfaces.

// VS 2013
#if _MSC_VER >= 1800

// VS 2017 v15.3
#if _MSC_VER >= 1911
#define C17_ENABLE
#define C17_MATH_ENABLE
#endif

#ifdef C17_MATH_ENABLE
#include <cmath>
#endif
#include <immintrin.h>
#include "avs/minmax.h"
#include "JincResizeMT.h"

#define myalignedfree(ptr) if (ptr!=nullptr) { _aligned_free(ptr); ptr=nullptr;}

static bool init_coeff_table(EWAPixelCoeff *out, int quantize_x, int quantize_y,
	int filter_size, int dst_width, int dst_height, int mod_align)
{
	out->filter_size = filter_size;
	if (mod_align > 0)
		out->coeff_stride = (filter_size + (mod_align - 1)) & ~(mod_align - 1);
	else
		out->coeff_stride = filter_size;

	// This will be reserved to exact size in coff generating procedure
	out->factor = nullptr;

	// Allocate metadata
	out->meta = new EWAPixelCoeffMeta[static_cast<int64_t>(dst_width) * dst_height];
	if (out->meta == nullptr) return(false);

	// Alocate factor map
	out->factor_map = new int[static_cast<int64_t>(quantize_x) * quantize_y];

	// Zeroed memory
	if (out->factor_map == nullptr) return(false);

	memset(out->factor_map, 0, static_cast<int64_t>(quantize_x) * quantize_y * sizeof(int));
	memset(out->meta, 0, static_cast<int64_t>(dst_width) * dst_height * sizeof(EWAPixelCoeffMeta));

	return(true);
}

static AVS_FORCEINLINE unsigned portable_clz(size_t x)
{
	unsigned long index;
	return (_BitScanReverse(&index, static_cast<unsigned long>(x))) ? (31 - index) : 32;
}

// Taylor series coefficients of 2*BesselJ1(pi*x)/(pi*x) as (x^2) -> 0
static const double jinc_taylor_series[31] =
{
	1.0,
	-1.23370055013616982735431137,
	0.507339015802096027273126733,
	-0.104317403816764804365258186,
	0.0128696438477519721233840271,
	-0.00105848577966854543020422691,
	6.21835470803998638484476598e-05,
	-2.73985272294670461142756204e-06,
	9.38932725442064547796003405e-08,
	-2.57413737759717407304931036e-09,
	5.77402672521402031756429343e-11,
	-1.07930605263598241754572977e-12,
	1.70710316782347356046974552e-14,
	-2.31434518382749184406648762e-16,
	2.71924659665997312120515390e-18,
	-2.79561335187943028518083529e-20,
	2.53599244866299622352138464e-22,
	-2.04487273140961494085786452e-24,
	1.47529860450204338866792475e-26,
	-9.57935105257523453155043307e-29,
	5.62764317309979254140393917e-31,
	-3.00555258814860366342363867e-33,
	1.46559362903641161989338221e-35,
	-6.55110024064596600335624426e-38,
	2.69403199029404093412381643e-40,
	-1.02265499954159964097119923e-42,
	3.59444454568084324694180635e-45,
	-1.17313973900539982313119019e-47,
	3.56478606255557746426034301e-50,
	-1.01100655781438313239513538e-52,
	2.68232117541264485328658605e-55
};

static const double jinc_zeros[16] =
{
	1.2196698912665045,
	2.2331305943815286,
	3.2383154841662362,
	4.2410628637960699,
	5.2427643768701817,
	6.2439216898644877,
	7.2447598687199570,
	8.2453949139520427,
	9.2458926849494673,
	10.246293348754916,
	11.246622794877883,
	12.246898461138105,
	13.247132522181061,
	14.247333735806849,
	15.247508563037300,
	16.247661874700962
};

#define EPS_JINC_PI 1e-6

#ifndef M_PI // GCC seems to have it
static const double M_PI = 3.14159265358979323846;
#endif

static const double JINC_ZERO_SQR = jinc_zeros[0] * jinc_zeros[0];

//  Modified from boost package math/tools/`rational.hpp`
//
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
static double evaluate_rational(const double *num, const double *denom, double z, int count)
{
	double s1, s2;

	if (std::fabs(z) <= 1.0)
	{
		s1 = num[count - 1];
		s2 = denom[count - 1];
		for (int i = count - 2; i >= 0; --i)
		{
			s1 *= z;
			s2 *= z;
			s1 += num[i];
			s2 += denom[i];
		}
	}
	else
	{
		z = 1.0 / z;
		s1 = num[0];
		s2 = denom[0];
		for (int i = 1; i < count; ++i)
		{
			s1 *= z;
			s2 *= z;
			s1 += num[i];
			s2 += denom[i];
		}
	}

	return(s1 / s2);
}

//  Modified from boost package `BesselJ1.hpp`
//
//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
static double jinc_sqr_boost_l(double x2)
{
	const double bPC[7] =
	{
		-4.4357578167941278571e+06,
		-9.9422465050776411957e+06,
		-6.6033732483649391093e+06,
		-1.5235293511811373833e+06,
		-1.0982405543459346727e+05,
		-1.6116166443246101165e+03,
		0.0
	};
	const double bQC[7] =
	{
		-4.4357578167941278568e+06,
		-9.9341243899345856590e+06,
		-6.5853394797230870728e+06,
		-1.5118095066341608816e+06,
		-1.0726385991103820119e+05,
		-1.4550094401904961825e+03,
		1.0
	};
	const double bPS[7] =
	{
		3.3220913409857223519e+04,
		8.5145160675335701966e+04,
		6.6178836581270835179e+04,
		1.8494262873223866797e+04,
		1.7063754290207680021e+03,
		3.5265133846636032186e+01,
		0.0
	};
	const double bQS[7] =
	{
		7.0871281941028743574e+05,
		1.8194580422439972989e+06,
		1.4194606696037208929e+06,
		4.0029443582266975117e+05,
		3.7890229745772202641e+04,
		8.6383677696049909675e+02,
		1.0
	};

	const double y2 = M_PI * M_PI*x2;
	const double xp = sqrt(y2);
	const double y2p = 64.0 / y2;
	const double sx = sin(xp);
	const double cx = cos(xp);

	return ((sqrt(xp / M_PI)*2.0 / y2)*(evaluate_rational(bPC, bQC, y2p, 7)*(sx - cx) + (8.0 / xp)*evaluate_rational(bPS, bQS, y2p, 7)*(sx + cx)));
}

#ifndef C17_MATH_ENABLE
#define MAX_TERMS 50
static double bessel_j1(double x)
{
    const double EPS = 1e-15;
	double TabD[MAX_TERMS];

    double term = x*0.5; // first value m=0, x/2

	TabD[0] = term;

    double x2 = (x*x)*0.25;

    unsigned long m=1;
	while ((std::fabs(term) >= EPS) && (m<MAX_TERMS))
	{
		term *= -x2/(m*(m+1)); // recurrence
		TabD[m++] = term;
	}

	double sum = TabD[m-1];

	for (int i=((int)m)-2; i>=0; i--)
		sum += TabD[i];
	
    return sum;
}
#endif

// jinc(sqrt(x2))
static double jinc_sqr(double x2)
{
	if (x2 < 1.49)        // the 1-tap radius
	{
		double res = 0.0;

		for (int j = 16; j > 0; --j)
			res = res * x2 + jinc_taylor_series[j - 1];

		return(res);
	}
	else if (x2 < 4.97)   // the 2-tap radius
	{
		double res = 0.0;

		for (int j = 21; j > 0; --j)
			res = res * x2 + jinc_taylor_series[j - 1];

		return(res);
	}
	else if (x2 < 10.49)  // the 3-tap radius
	{
		double res = 0.0;

		for (int j = 26; j > 0; --j)
			res = res * x2 + jinc_taylor_series[j - 1];

		return(res);
	}
	else if (x2 < 17.99)  // the 4-tap radius
	{
		double res = 0.0;

		for (int j = 31; j > 0; --j)
			res = res * x2 + jinc_taylor_series[j - 1];

		return(res);
	}
	else if (x2 < 52.57)  // the 5~7-tap radius
	{
		const double x = M_PI * sqrt(x2);
#ifdef C17_MATH_ENABLE
		return(2.0*std::cyl_bessel_j(1, x) / x);
#else
		return(2.0*bessel_j1(x) / x);
#endif
	}
	else if (x2 < 68.07)  // the 8-tap radius // Modify from pull request #4
	{
		return(jinc_sqr_boost_l(x2));
	}
	else                  // the 9~16-tap radius
	{
		const double x = M_PI * sqrt(x2);

#ifdef C17_MATH_ENABLE
		return(2.0*std::cyl_bessel_j(1, x) / x);
#else
		return(2.0*bessel_j1(x) / x);
#endif
	}
}

static double sample_sqr(double(*filter)(double), double x2, double blur2, double radius2)
{
	if (blur2 > 0.0)
		x2 /= blur2;

	if (x2 < radius2)
		return(filter(x2));

	return(0.0);
}

static double GetFactor2D(double dx, double dy, double radius, double blur, WEIGHTING_TYPE wt)
{
	const double arg2 = dx * dx + dy * dy;
	double arg = sqrt(arg2);

	if (arg > radius) arg = radius; // limit for safety ?

	const double blur2 = blur * blur;
	const double radius2 = radius * radius;

	switch (wt)
	{
	case SP_WT_NONE:
		return(sample_sqr(jinc_sqr, arg2, blur2, radius2));
		break;
	case SP_WT_JINC:
		return(sample_sqr(jinc_sqr, arg2, blur2, radius2)*sample_sqr(jinc_sqr, JINC_ZERO_SQR*(arg2 / radius2), 1.0, radius2));
		break;
	case SP_WT_TRD2:
		if (arg < (radius / 2))
			return(sample_sqr(jinc_sqr, arg2, blur2, radius2));
		else
			return(sample_sqr(jinc_sqr, arg2, blur2, radius2)*((2.0 - (2.0*(arg / radius)))));
		break;
	default: return(0.0); break;
	}
}

// Here is our simple jinc_pi(x), not 2.0x because of auto-normalizing of kernel for convolution in resampling program generator, but M_PI scaled argument
inline static double jinc_pi(double arg)
{
	arg *= M_PI;

	if (std::fabs(arg) < EPS_JINC_PI)
	{
		const double a = -1.0 / 16.0, b = 1.0 / 384.0;

		arg *= arg;
		return((arg*b + a)*arg + 0.5);
	}
	else
	{
#ifdef C17_ENABLE
		return(std::cyl_bessel_j(1, arg) / arg);
#else
		return(bessel_j1(arg) / arg);
#endif
	}
}

// jinc function value from x,y origin to dx, dy coordinate from origin 
inline static double jpdi(double dx, double dy, int x, int y)
{
	double dist = sqrt((x - dx)*(x - dx) + (y - dy)*(y - dy));
	return jinc_pi(dist);
}

/*
2D kernel of sum of jincs of max size 5x5 with trimmed out corner samples (XX), so 21 jincs in sum total
kernel samples placement in 2D full numbering (x,y)
where k(+0,+0) = 1.0 - center sample of kernel
XX       k(-1,+2) k(+0,+2) k(+1,+2) XX
k(-2,+1) k(-1,+1) k(+0,+1) k(+1,+1) k(+2,+1)
k(-2,+0) k(-1,+0) k(+0,+0) k(+1,+0) k(+2,+0)
k(-2,-1) k(-1,-1) k(+0,-1) k(+1,-1) k(+2,-1)
XX       k(-1,-2) k(+0,-2) k(+1,-2) XX

copy of k10, k20, k11, k21 by symmethry:
XX       k21      k20      k21      XX
k21      k11      k10      k11      k21
k20      k10      1.0      k10      k20
k21      k11      k10      k11      k21
XX       k21      k20      k21      XX
*/
static double GetFactor2D_JINCSUM_21(double dx, double dy, float k10, float k20, float k11, float k21, double radius_sq)
{
	const double dist_sq = dx * dx + dy * dy;

	if (dist_sq > radius_sq) return(0.0); // make kernel round, may be option to be square in the future

	auto sum = 0.0;
	// 1st row
	sum += jpdi(dx, dy, -1, +2)*k21 + jpdi(dx, dy, +0, +2)*k20 + jpdi(dx, dy, +1, +2)*k21;
	// 2nd row
	sum += jpdi(dx, dy, -2, +1)*k21 + jpdi(dx, dy, -1, +1)*k11 + jpdi(dx, dy, +0, +1)*k10 + jpdi(dx, dy, +1, +1)*k11 + jpdi(dx, dy, +2, +1)*k21;
	// 3rd row
	sum += jpdi(dx, dy, -2, +0)*k20 + jpdi(dx, dy, -1, +0)*k10 + jpdi(dx, dy, +0, +0)*1.0 + jpdi(dx, dy, +1, +0)*k10 + jpdi(dx, dy, +2, +0)*k20;
	// 4th row
	sum += jpdi(dx, dy, -2, -1)*k21 + jpdi(dx, dy, -1, -1)*k11 + jpdi(dx, dy, +0, -1)*k10 + jpdi(dx, dy, +1, -1)*k11 + jpdi(dx, dy, +2, -1)*k21;
	// 5th row
	sum += jpdi(dx, dy, -1, -2)*k21 + jpdi(dx, dy, +0, -2)*k20 + jpdi(dx, dy, +1, -2)*k21;

	return sum;
}

/* Coefficient table generation for fp16*/
#if defined(CLANG)
__attribute__((__target__("avx2,f16c")))
#endif
bool generate_coeff_table_fp16_c(const JincMT_generate_coeff_params& params)
{
	JincMT_Lut *func = params.func;
	EWAPixelCoeff* out_fp16 = params.out_fp16;
	int quantize_x = params.quantize_x;
	int quantize_y = params.quantize_y;
	int samples = params.samples;
	int src_width = params.src_width;
	int src_height = params.src_height;
	int dst_width = params.dst_width;
	int dst_height = params.dst_height;
	double radius = params.radius;
	int mod_align = params.mod_align;

	const float k10 = params.k10;
	const float k20 = params.k20;
	const float k11 = params.k11;
	const float k21 = params.k21;
	SP_KERNEL_TYPE kernel_type = params.kernel_type;

	const double filter_step_x = min(static_cast<double>(dst_width) / params.crop_width, 1.0);
	const double filter_step_y = min(static_cast<double>(dst_height) / params.crop_height, 1.0);

	const float filter_support_x = static_cast<float>(radius / filter_step_x);
	const float filter_support_y = static_cast<float>(radius / filter_step_y);

	const float filter_support = max(filter_support_x, filter_support_y);
	const int filter_size = max(static_cast<int>(ceil(filter_support_x * 2.0)), static_cast<int>(ceil(filter_support_y * 2.0)));

	const float start_x = static_cast<float>(params.crop_left + (params.crop_width / dst_width - 1.0) / 2.0);

	const float x_step = static_cast<float>(params.crop_width / dst_width);
	const float y_step = static_cast<float>(params.crop_height / dst_height);

	float xpos = start_x;
	float ypos = static_cast<float>(params.crop_top + (params.crop_height - dst_height) / (dst_height * static_cast<int64_t>(2)));

	// Initialize EWAPixelCoeff data structure
	if (!init_coeff_table(out_fp16, quantize_x, quantize_y, filter_size, dst_width, dst_height, mod_align)) return(false);
	out_fp16->coeff_stride /= 2; // for 16bit samples

	size_t tmp_array_capacity = out_fp16->coeff_stride * 2 * filter_size; // in number of float32
	float* tmp_array = static_cast<float*>(_aligned_malloc(tmp_array_capacity * sizeof(float), 64));
	if (tmp_array == nullptr) return(false);

	size_t tmp_array_fp16_capacity = params.initial_capacity / 2;
	float* tmp_array_fp16 = static_cast<float*>(_aligned_malloc(tmp_array_fp16_capacity * sizeof(float), 64));
	if (tmp_array_fp16 == nullptr) return(false);
	size_t tmp_array_fp16_size = 0;
	int tmp_array_fp16_top = 0;
	unsigned base_clz_fp16 = portable_clz(tmp_array_fp16_capacity);
	const double initial_growth_factor = params.initial_factor;

	const double radius2 = radius * radius;

	// Use to advance the coeff pointer
	const int coeff_per_pixel_fp16 = out_fp16->coeff_stride * filter_size;

	for (int y = 0; y < dst_height; ++y)
	{
		for (int x = 0; x < dst_width; ++x)
		{
			bool is_border = false;

			EWAPixelCoeffMeta* meta_fp16 = &out_fp16->meta[y * dst_width + x];

			// Here, the window_*** variable specified a begin/size/end
			// of EWA window to process.
			int window_end_x = static_cast<int>(xpos + filter_support);
			int window_end_y = static_cast<int>(ypos + filter_support);

			if (window_end_x >= src_width)
			{
				window_end_x = src_width - 1;
				is_border = true;
			}
			if (window_end_y >= src_height)
			{
				window_end_y = src_height - 1;
				is_border = true;
			}

			int window_begin_x = window_end_x - filter_size + 1;
			int window_begin_y = window_end_y - filter_size + 1;

			if (window_begin_x < 0)
			{
				window_begin_x = 0;
				is_border = true;
			}
			if (window_begin_y < 0)
			{
				window_begin_y = 0;
				is_border = true;
			}

			meta_fp16->start_x = window_begin_x;
			meta_fp16->start_y = window_begin_y;

			// Quantize xpos and ypos
			const int quantized_x_int = static_cast<int>(xpos * quantize_x);
			const int quantized_y_int = static_cast<int>(ypos * quantize_y);
			const int quantized_x_value = quantized_x_int % quantize_x;
			const int quantized_y_value = quantized_y_int % quantize_y;
			const float quantized_xpos = static_cast<float>(quantized_x_int) / quantize_x;
			const float quantized_ypos = static_cast<float>(quantized_y_int) / quantize_y;

			if (!is_border && out_fp16->factor_map[quantized_y_value * quantize_x + quantized_x_value] != 0)
			{
				// Not border pixel and already have coefficient calculated at this quantized position
				meta_fp16->coeff_meta = out_fp16->factor_map[quantized_y_value * quantize_x + quantized_x_value] - 1;
			}
			else
			{
				// then need computation
				float divider = 0.0f;

				// This is the location of current target pixel in source pixel
				// Quantized
				//const float current_x = clamp(is_border ? xpos : quantized_xpos, 0.f, src_width - 1.f);
				//const float current_y = clamp(is_border ? ypos : quantized_ypos, 0.f, src_height - 1.f);

				if (!is_border)
				{
					// Change window position to quantized position
					window_begin_x = static_cast<int>(quantized_xpos + filter_support) - filter_size + 1;
					window_begin_y = static_cast<int>(quantized_ypos + filter_support) - filter_size + 1;
				}

				// Windowing positon
				int window_x = window_begin_x;
				int window_y = window_begin_y;

				// First loop calcuate coeff
				int curr_factor_ptr = 0;
				memset(tmp_array, 0, tmp_array_capacity * sizeof(float)); // clean with zeroes

				// fp16 array
				const size_t new_size_fp16 = tmp_array_fp16_size + coeff_per_pixel_fp16;
				if (new_size_fp16 > tmp_array_fp16_capacity)
				{
					size_t new_capacity_fp16 = (size_t)(tmp_array_fp16_capacity * (1.0 + (initial_growth_factor - 1.0)
						* (1.0 - static_cast<double>(max(0, static_cast<int>(base_clz_fp16 - portable_clz(tmp_array_fp16_capacity)))) / 32.0)));
					if (new_capacity_fp16 < new_size_fp16)
						new_capacity_fp16 = new_size_fp16;
					float* new_tmp_fp16 = static_cast<float*>(_aligned_malloc(new_capacity_fp16 * sizeof(float), 64));
					if (new_tmp_fp16 == nullptr)
					{
						myalignedfree(tmp_array_fp16);
						return(false);
					}
					memcpy(new_tmp_fp16, tmp_array_fp16, tmp_array_fp16_size * sizeof(float));
					myalignedfree(tmp_array_fp16);
					tmp_array_fp16 = new_tmp_fp16;
					tmp_array_fp16_capacity = new_capacity_fp16;
				}
				memset(tmp_array_fp16 + tmp_array_fp16_size, 0, coeff_per_pixel_fp16 * sizeof(float));
				int curr_factor_ptr_fp16 = tmp_array_fp16_top;
				tmp_array_fp16_size = new_size_fp16;

				for (int ly = 0; ly < filter_size; ++ly)
				{
					for (int lx = 0; lx < filter_size; ++lx)
					{
						// Euclidean distance to sampling pixel
						const double dx = (clamp(is_border ? xpos : quantized_xpos, 0.0f, static_cast<float>(src_width - 1)) - window_x) * filter_step_x;
						const double dy = (clamp(is_border ? ypos : quantized_ypos, 0.0f, static_cast<float>(src_height - 1)) - window_y) * filter_step_y;

						float factor;

						switch (kernel_type)
						{
						case SP_JINCSINGLE:
							if (params.bUseLUTkernel)
							{
								//int index = static_cast<int>(llround((samples-1)*(dx*dx+dy*dy)/radius2 + DOUBLE_ROUND_MAGIC_NUMBER));
								const int index = static_cast<int>(llround((samples - 1) * (dx * dx + dy * dy) / radius2));
								factor = func->GetFactor(index);
							}
							else
								factor = (float)GetFactor2D(dx, dy, radius, params.blur, params.weighting_type);
							break;
						case SP_JINCSUM:
							factor = (float)GetFactor2D_JINCSUM_21(dx, dy, k10, k20, k11, k21, radius2);
							break;
						default: factor = 0.0; break;
						}

						tmp_array[curr_factor_ptr + static_cast<int64_t>(lx)] = factor;
						divider += factor;

						++window_x;
					}

					curr_factor_ptr += (out_fp16->coeff_stride * 2); // in size of float32

					window_x = window_begin_x;
					++window_y;
				}

				// Second loop to divide the coeff
				curr_factor_ptr = 0;
				for (int ly = 0; ly < filter_size; ++ly)
				{
					for (int lx = 0; lx < filter_size; ++lx)
					{
						tmp_array[curr_factor_ptr + static_cast<int64_t>(lx)] /= divider;
					}

					curr_factor_ptr += (out_fp16->coeff_stride * 2); // in size of float32
				}

				// Save factor to table
				if (!is_border)
				{
					out_fp16->factor_map[quantized_y_value * quantize_x + quantized_x_value] = tmp_array_fp16_top + 1;
				}

				// convert copy to fp16
				curr_factor_ptr = 0;
				curr_factor_ptr_fp16 = tmp_array_fp16_top;
				for (int cy = 0; cy < filter_size; ++cy)
				{
					for (int cx = 0; cx < filter_size; cx += 8)
					{
						const __m256 coeff = _mm256_load_ps(tmp_array + curr_factor_ptr + cx);
						//const __m128i coeff_fp16 = _mm256_cvtps_ph(coeff, _MM_FROUND_NO_EXC);
						const __m128i coeff_fp16 = _mm256_cvtps_ph(coeff, 0);
						_mm_store_si128((__m128i*)(tmp_array_fp16 + curr_factor_ptr_fp16 + (cx / 2)), coeff_fp16);
					}

					//curr_factor_ptr += out->coeff_stride;
					curr_factor_ptr += (out_fp16->coeff_stride) * 2; // in size of float32
					curr_factor_ptr_fp16 += out_fp16->coeff_stride;
				}

				meta_fp16->coeff_meta = tmp_array_fp16_top;
				tmp_array_fp16_top += coeff_per_pixel_fp16;

			}

			xpos += x_step;
		}

		ypos += y_step;
		xpos = start_x;
	}

	// Copy from tmp_array to real array
	out_fp16->factor = tmp_array_fp16;

	// free fp32 array
	myalignedfree(tmp_array);

	return(true);
}

template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx2,fma,f16c")))
#endif
void resize_plane_avx2_1x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx = 0;

	const T *src = reinterpret_cast<const T*>(MT_DataGF->src[idx]);
	T *JincMT_RESTRICT dst = reinterpret_cast<T*>(MT_DataGF->dst[idx]);

	const ptrdiff_t src_pitch = (ptrdiff_t)MT_DataGF->src_pitch[idx] / sizeof(T);
	const ptrdiff_t dst_pitch = (ptrdiff_t)MT_DataGF->dst_pitch[idx] / sizeof(T);
	const ptrdiff_t src_pitchx2 = src_pitch << 1;

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min = static_cast<T>(Val_Min[idx]);
	const T val_max = static_cast<T>(Val_Max[idx]);
	const __m256 min_val = _mm256_set1_ps(Val_Min[idx]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

    for (int y = Y_Min; y < Y_Max; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;

        for (int x = 0; x < dst_width; ++x)
        {
            const T *src_ptr = src + (meta->start_y * src_pitch + meta->start_x);
            const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
            __m256 result = _mm256_setzero_ps();
			__m256 result2 = _mm256_setzero_ps();

            if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
            {
                for (int ly = 0; ly < filter_size_mod2; ly += 2)
                {
                    for (int lx = 0; lx < filter_size; lx += 8)
                    {
                        const __m256 src_ps = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx + src_pitch)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							const __m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result = _mm256_fmadd_ps(src_ps, coeff, result);
						result2 = _mm256_fmadd_ps(src_ps2, coeff2, result2);
                    }

                    coeff_ptr += coeff_stridex2;
                    src_ptr += src_pitchx2;
                }

				result = _mm256_add_ps(result, result2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result = _mm256_fmadd_ps(src_ps, coeff, result);
					}
				}

                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(result), _mm256_extractf128_ps(result, 1));
                hsum = _mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()), _mm_setzero_si128()));

				dst[x] = clamp(final_res,val_min,val_max);
            }
            else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
            {
                for (int ly = 0; ly < filter_size_mod2; ly += 2)
                {
                    for (int lx = 0; lx < filter_size; lx += 8)
                    {
						const __m256 src_ps = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx + src_pitch)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result = _mm256_fmadd_ps(src_ps, coeff, result);
						result2 = _mm256_fmadd_ps(src_ps2, coeff2, result2);
                    }

					coeff_ptr += coeff_stridex2;
                    src_ptr += src_pitchx2;
                }

				result = _mm256_add_ps(result, result2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result = _mm256_fmadd_ps(src_ps, coeff, result);
					}
				}

                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(result), _mm256_extractf128_ps(result, 1));
                hsum = _mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()));

				dst[x] = clamp(final_res,val_min,val_max);
            }
            else
            {
                for (int ly = 0; ly < filter_size_mod2; ly += 2)
                {
                    for (int lx = 0; lx < filter_size; lx += 8)
                    {
						const __m256 src_ps = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr + lx)), min_val);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr + lx + src_pitch)), min_val);
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result = _mm256_fmadd_ps(src_ps, coeff, result);
						result2 = _mm256_fmadd_ps(src_ps2, coeff2, result2);
					}

                    coeff_ptr += coeff_stridex2;
                    src_ptr += src_pitchx2;
                }

				result = _mm256_add_ps(result, result2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr + lx)), min_val);
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result = _mm256_fmadd_ps(src_ps, coeff, result);
					}
				}

                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(result), _mm256_extractf128_ps(result, 1));

				dst[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum)));
            }
			meta++;
        } // for (x)
		meta_y += dst_width;
        dst += dst_pitch;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx2,f16c")))
#endif
void resize_plane_avx2_2x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = (PlaneYMode) ? 0 : 1;
	const uint8_t idx2 = (PlaneYMode) ? 3 : 2;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch1x2 = src_pitch1 << 1;
	const ptrdiff_t src_pitch2x2 = src_pitch2 << 1;

	const int Y_Min = (PlaneYMode) ? MT_DataGF->dst_Y_h_min : MT_DataGF->dst_UV_h_min;
	const int Y_Max = (PlaneYMode) ? MT_DataGF->dst_Y_h_max : MT_DataGF->dst_UV_h_max;
	const int dst_width = (PlaneYMode) ? MT_DataGF->dst_Y_w : MT_DataGF->dst_UV_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const __m256 min_val1 = _mm256_set1_ps(Val_Min[idx1]);
	const __m256 min_val2 = _mm256_set1_ps(Val_Min[idx2]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m256 result1 = _mm256_setzero_ps();
			__m256 result2 = _mm256_setzero_ps();
			__m256 result1_2 = _mm256_setzero_ps();
			__m256 result2_2 = _mm256_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps1_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps2_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps1_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps2_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
			}
			else
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps1_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx + src_pitch1)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val1);
						const __m256 src_ps2_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx + src_pitch2)), min_val2);
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val2);
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx2,f16c")))
#endif
void resize_plane_avx2_3x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = 0;
	const uint8_t idx2 = 1;
	const uint8_t idx3 = 2;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	const T *src3 = reinterpret_cast<const T*>(MT_DataGF->src[idx3]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);
	T *JincMT_RESTRICT dst3 = reinterpret_cast<T*>(MT_DataGF->dst[idx3]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch3 = (ptrdiff_t)MT_DataGF->src_pitch[idx3] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch3 = (ptrdiff_t)MT_DataGF->dst_pitch[idx3] / sizeof(T);
	const ptrdiff_t src_pitch1x2 = src_pitch1 << 1;
	const ptrdiff_t src_pitch2x2 = src_pitch2 << 1;
	const ptrdiff_t src_pitch3x2 = src_pitch3 << 1;

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_min3 = static_cast<T>(Val_Min[idx3]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const T val_max3 = static_cast<T>(Val_Max[idx3]);
	const __m256 min_val1 = _mm256_set1_ps(Val_Min[idx1]);
	const __m256 min_val2 = _mm256_set1_ps(Val_Min[idx2]);
	const __m256 min_val3 = _mm256_set1_ps(Val_Min[idx3]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m256 result1 = _mm256_setzero_ps();
			__m256 result2 = _mm256_setzero_ps();
			__m256 result3 = _mm256_setzero_ps();
			__m256 result1_2 = _mm256_setzero_ps();
			__m256 result2_2 = _mm256_setzero_ps();
			__m256 result3_2 = _mm256_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps1_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps2_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps3_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx + src_pitch2)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm256_fmadd_ps(src_ps3_2, coeff2, result3_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);
				result3 = _mm256_add_ps(result3, result3_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps1_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps2_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps3_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx + src_pitch2)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm256_fmadd_ps(src_ps3_2, coeff2, result3_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);
				result3 = _mm256_add_ps(result3, result3_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
			}
			else
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps1_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx + src_pitch1)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val1);
						const __m256 src_ps2_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx + src_pitch2)), min_val2);
						const __m256 src_ps3 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx)), min_val3);
						const __m256 src_ps3_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx + src_pitch3)), min_val3);
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm256_fmadd_ps(src_ps3_2, coeff2, result3_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);
				result3 = _mm256_add_ps(result3, result3_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val2);
						const __m256 src_ps3 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx)), min_val3);
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
				dst3[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx2,f16c")))
#endif
void resize_plane_avx2_4x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = 0;
	const uint8_t idx2 = 1;
	const uint8_t idx3 = 2;
	const uint8_t idx4 = 3;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	const T *src3 = reinterpret_cast<const T*>(MT_DataGF->src[idx3]);
	const T *src4 = reinterpret_cast<const T*>(MT_DataGF->src[idx4]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);
	T *JincMT_RESTRICT dst3 = reinterpret_cast<T*>(MT_DataGF->dst[idx3]);
	T *JincMT_RESTRICT dst4 = reinterpret_cast<T*>(MT_DataGF->dst[idx4]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch3 = (ptrdiff_t)MT_DataGF->src_pitch[idx3] / sizeof(T);
	const ptrdiff_t src_pitch4 = (ptrdiff_t)MT_DataGF->src_pitch[idx4] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch3 = (ptrdiff_t)MT_DataGF->dst_pitch[idx3] / sizeof(T);
	const ptrdiff_t dst_pitch4 = (ptrdiff_t)MT_DataGF->dst_pitch[idx4] / sizeof(T);

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_min3 = static_cast<T>(Val_Min[idx3]);
	const T val_min4 = static_cast<T>(Val_Min[idx4]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const T val_max3 = static_cast<T>(Val_Max[idx3]);
	const T val_max4 = static_cast<T>(Val_Max[idx4]);
	const __m256 min_val1 = _mm256_set1_ps(Val_Min[idx1]);
	const __m256 min_val2 = _mm256_set1_ps(Val_Min[idx2]);
	const __m256 min_val3 = _mm256_set1_ps(Val_Min[idx3]);
	const __m256 min_val4 = _mm256_set1_ps(Val_Min[idx4]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const T *src_ptr4 = src4 + (meta->start_y * src_pitch4 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m256 result1 = _mm256_setzero_ps();
			__m256 result2 = _mm256_setzero_ps();
			__m256 result3 = _mm256_setzero_ps();
			__m256 result4 = _mm256_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size; ++ly)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps4 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
					}

					coeff_ptr += coeff_stride;
					src_ptr1 += src_pitch1;
					src_ptr2 += src_pitch2;
					src_ptr3 += src_pitch3;
					src_ptr4 += src_pitch4;
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				hsum4 = _mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res4 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum4), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
				dst4[x] = clamp(final_res4, val_min4, val_max4);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size; ++ly)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps4 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
					}

					coeff_ptr += coeff_stride;
					src_ptr1 += src_pitch1;
					src_ptr2 += src_pitch2;
					src_ptr3 += src_pitch3;
					src_ptr4 += src_pitch4;
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				hsum4 = _mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()));
				const T final_res4 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum4), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
				dst4[x] = clamp(final_res4, val_min4, val_max4);
			}
			else
			{
				for (int ly = 0; ly < filter_size; ++ly)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val2);
						const __m256 src_ps3 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx)), min_val3);
						const __m256 src_ps4 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr4 + lx)), min_val4);
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
					}

					coeff_ptr += coeff_stride;
					src_ptr1 += src_pitch1;
					src_ptr2 += src_pitch2;
					src_ptr3 += src_pitch3;
					src_ptr4 += src_pitch4;
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
				dst3[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3)));
				dst4[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
		dst4 += dst_pitch4;
	} // for (y)
}


template void resize_plane_avx2_1x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_1x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_1x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_1x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_1x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_1x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_2x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_2x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_2x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_2x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_2x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_2x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_3x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_3x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_3x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_3x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_3x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_3x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_4x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_4x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

#endif