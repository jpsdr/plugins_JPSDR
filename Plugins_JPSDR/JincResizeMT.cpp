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

#include <cmath>
#include "avs/minmax.h"

// VS 2013
#if _MSC_VER >= 1800
#define AVX2_BUILD_POSSIBLE
#endif

// VS 2017 v15.3
#if _MSC_VER >= 1911
#define AVX512_BUILD_POSSIBLE
#define C17_ENABLE
#endif

#include "JincResizeMT.h"
#include "resize_plane_sse41.h"

#ifdef AVX2_BUILD_POSSIBLE
#include "resize_plane_avx2.h"
#endif

#ifdef AVX512_BUILD_POSSIBLE
#include "resize_plane_avx512.h"
#endif

#define myfree(ptr) if (ptr!=nullptr) { free(ptr); ptr=nullptr;}
#define myalignedfree(ptr) if (ptr!=nullptr) { _aligned_free(ptr); ptr=nullptr;}
#define mydeleteT(ptr) if (ptr!=nullptr) { delete[] ptr; ptr=nullptr;}
#define mydelete(ptr) if (ptr!=nullptr) { delete ptr; ptr=nullptr;}

extern ThreadPoolInterface *poolInterface;

static bool is_paramstring_empty_or_auto(const char *param)
{
	if (param == nullptr) return true;
	return (_stricmp(param, "auto") == 0); // true is match
}

static bool getChromaLocation(const char *chromaloc_name, IScriptEnvironment *env, ChromaLocation_e &_ChromaLocation)
{
	ChromaLocation_e index = AVS_CHROMA_UNUSED;

	if (_stricmp(chromaloc_name, "left") == 0) index = AVS_CHROMA_LEFT;
	if (_stricmp(chromaloc_name, "center") == 0) index = AVS_CHROMA_CENTER;
	if ((_stricmp(chromaloc_name, "top_left") == 0) || (_stricmp(chromaloc_name, "topleft") == 0))
		index = AVS_CHROMA_TOP_LEFT;
	if (_stricmp(chromaloc_name, "top") == 0) index = AVS_CHROMA_TOP; // not used in Avisynth
	if ((_stricmp(chromaloc_name, "bottom_left") == 0) || (_stricmp(chromaloc_name, "bottomleft") == 0))
		index = AVS_CHROMA_BOTTOM_LEFT; // not used in Avisynth
	if (_stricmp(chromaloc_name, "bottom") == 0) index = AVS_CHROMA_BOTTOM; // not used in Avisynth
	if (_stricmp(chromaloc_name, "dv") == 0) index = AVS_CHROMA_DV; // Special to Avisynth
	// compatibility
	if (_stricmp(chromaloc_name, "mpeg1") == 0) index = AVS_CHROMA_CENTER;
	if (_stricmp(chromaloc_name, "mpeg2") == 0) index = AVS_CHROMA_LEFT;
	if (_stricmp(chromaloc_name, "jpeg") == 0) index = AVS_CHROMA_CENTER;

	if (index != AVS_CHROMA_UNUSED)
	{
		_ChromaLocation = index;
		return true;
	}

	env->ThrowError("JincResizeMT: Unknown chroma placement");
	// empty
	return false;
}

static void chromaloc_parse_merge_with_props(const VideoInfo &vi, const char *chromaloc_name, const AVSMap *props, ChromaLocation_e &_ChromaLocation, ChromaLocation_e _ChromaLocation_Default, IScriptEnvironment *env)
{
	if (props != nullptr)
	{
		if (vi.Is420() || vi.Is422() || vi.IsYV411())
		{ // yes, YV411 can also have valid _ChromaLocation, if 'left'-ish one is given
			if (env->propNumElements(props, "_ChromaLocation") > 0)
				_ChromaLocation_Default = (ChromaLocation_e)env->propGetIntSaturated(props, "_ChromaLocation", 0, nullptr);
		}
		else
		{
			// Theoretically RGB and not subsampled formats must not have chroma location
			if (env->propNumElements(props, "_ChromaLocation") > 0)
			{
				// Uncommented for a while, just ignore when there is any
				// env->ThrowError("Error: _ChromaLocation property found at a non-subsampled source.");
			}
		}
	}

	if (is_paramstring_empty_or_auto(chromaloc_name) || !getChromaLocation(chromaloc_name, env, _ChromaLocation))
		_ChromaLocation = _ChromaLocation_Default;
}

static AVS_FORCEINLINE unsigned portable_clz(size_t x)
{
    unsigned long index;
    return (_BitScanReverse(&index, static_cast<unsigned long>(x))) ? (31 - index) : 32;
}


#ifndef M_PI // GCC seems to have it
static double M_PI = 3.14159265358979323846;
#endif

// Taylor series coefficients of 2*BesselJ1(pi*x)/(pi*x) as (x^2) -> 0
static double jinc_taylor_series[31] =
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

static double jinc_zeros[16] =
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

double JINC_ZERO_SQR = jinc_zeros[0] * jinc_zeros[0];

//  Modified from boost package math/tools/`rational.hpp`
//
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
static double evaluate_rational(const double* num, const double* denom, double z, int count)
{
    double s1, s2;
    if (z <= 1.0)
    {
        s1 = num[count - 1];
        s2 = denom[count - 1];
        for (auto i = count - 2; i >= 0; --i)
        {
            s1 *= z;
            s2 *= z;
            s1 += num[i];
            s2 += denom[i];
        }
    }
    else
    {
        z = 1.0f / z;
        s1 = num[0];
        s2 = denom[0];
        for (auto i = 1; i < count; ++i)
        {
            s1 *= z;
            s2 *= z;
            s1 += num[i];
            s2 += denom[i];
        }
    }

    return s1 / s2;
}

//  Modified from boost package `BesselJ1.hpp`
//
//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
static double jinc_sqr_boost_l(double x2)
{
    double bPC[7] =
    {
        -4.4357578167941278571e+06,
        -9.9422465050776411957e+06,
        -6.6033732483649391093e+06,
        -1.5235293511811373833e+06,
        -1.0982405543459346727e+05,
        -1.6116166443246101165e+03,
        0.0
    };
    double bQC[7] =
    {
        -4.4357578167941278568e+06,
        -9.9341243899345856590e+06,
        -6.5853394797230870728e+06,
        -1.5118095066341608816e+06,
        -1.0726385991103820119e+05,
        -1.4550094401904961825e+03,
        1.0
    };
    double bPS[7] =
    {
        3.3220913409857223519e+04,
        8.5145160675335701966e+04,
        6.6178836581270835179e+04,
        1.8494262873223866797e+04,
        1.7063754290207680021e+03,
        3.5265133846636032186e+01,
        0.0
    };
    double bQS[7] =
    {
        7.0871281941028743574e+05,
        1.8194580422439972989e+06,
        1.4194606696037208929e+06,
        4.0029443582266975117e+05,
        3.7890229745772202641e+04,
        8.6383677696049909675e+02,
        1.0
    };

    const auto y2 = M_PI * M_PI * x2;
    const auto xp = sqrt(y2);
    const auto y2p = 64.0 / y2;
    const auto sx = sin(xp);
    const auto cx = cos(xp);

    return ((sqrt(xp/M_PI)*2.0/y2)*(evaluate_rational(bPC,bQC,y2p,7)*(sx-cx)+(8.0/xp)*evaluate_rational(bPS,bQS,y2p,7)*(sx+cx)));
}

#ifndef C17_ENABLE
#define MAX_TERMS 50
static double bessel_j1(double x)
{
    const double EPS = 1e-15;
	double TabD[MAX_TERMS];

    double term = x / 2.0; // first value m=0

	TabD[0] = term;

    double x2 = (x * x) / 4.0;

    unsigned long m=1;
	while ((std::fabs(term) >= EPS) && (m<MAX_TERMS))
	{
		term *= -x2 / (m * (m + 1)); // recurrence
		TabD[m++] = term;
	}

	double sum = 0.0;

	for (int i=m-1; i>=0; i--)
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
        for (auto j = 16; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 4.97)   // the 2-tap radius
    {
        double res = 0.0;
        for (auto j = 21; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 10.49)  // the 3-tap radius
    {
        double res = 0.0;
        for (auto j = 26; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 17.99)  // the 4-tap radius
    {
        double res = 0.0;
        for (auto j = 31; j > 0; --j)
            res = res * x2 + jinc_taylor_series[j - 1];
        return res;
    }
    else if (x2 < 52.57)  // the 5~7-tap radius
    {
        const auto x = M_PI * sqrt(x2);
#ifdef C17_ENABLE
        return 2.0 * std::cyl_bessel_j(1, x) / x;
#else
		return 2.0 * bessel_j1(x)/x;
#endif
    }
    else if (x2 < 68.07)  // the 8-tap radius // Modify from pull request #4
    {
        return jinc_sqr_boost_l(x2);
    }
    else                  // the 9~16-tap radius
    {
        const auto x = M_PI * sqrt(x2);
#ifdef C17_ENABLE
        return 2.0 * std::cyl_bessel_j(1, x) / x;
#else
		return 2.0 * bessel_j1(x)/x;
#endif
    }
}


static double sample_sqr(double (*filter)(double), double x2, double blur2, double radius2)
{
    if (blur2 > 0.0)
        x2 /= blur2;

    if (x2 < radius2)
        return filter(x2);

    return 0.0;
}


Lut::Lut() : lut_size(LUT_SIZE_VALUE)
{
    lut = new double[lut_size];
}

Lut::~Lut()
{
	mydelete(lut);
}

bool Lut::InitLut(int lutsize, double radius, double blur, WEIGHTING_TYPE wt)
{
	if ((lut == nullptr) || (lutsize > lut_size)) return (false);

	const auto radius2 = radius * radius;
    const auto blur2 = blur * blur;

    for (auto i = 0; i < lutsize; ++i)
    {
        const auto t2 = i / (lutsize - 1.0);
		
		switch(wt)
		{
			case SP_WT_NONE :
				lut[i] = sample_sqr(jinc_sqr, radius2 * t2, blur2, radius2);
				break;
			case SP_WT_JINC :
				lut[i] = sample_sqr(jinc_sqr, radius2 * t2, blur2, radius2) * sample_sqr(jinc_sqr, JINC_ZERO_SQR * t2, 1.0, radius2);
				break;
			case SP_WT_TRD2 :
				if (i < (lutsize / 2))
					lut[i] = sample_sqr(jinc_sqr, radius2 * t2, blur2, radius2);
				else
					lut[i] = sample_sqr(jinc_sqr, radius2 * t2, blur2, radius2) * ((2.0 - (2.0 * t2 )));
				break;
			default : return(false); break;
		}
    }

	return(true);
}

float Lut::GetFactor(int index)
{
    if ((index >= lut_size) || (lut == nullptr))
        return 0.0f;
    return static_cast<float>(lut[index]);
}


static double GetFactor2D(double dx, double dy, double radius, double blur, WEIGHTING_TYPE wt)
{
	const auto arg2 = dx*dx + dy*dy;
	auto arg = sqrt(arg2);

	if (arg > radius) arg = radius; // limit for safety ?

	const auto blur2 = blur * blur;
	const auto radius2 = radius * radius;

	switch(wt)
	{
		case SP_WT_NONE :
			return(sample_sqr(jinc_sqr,arg2,blur2,radius2));
			break;
		case SP_WT_JINC :
			return(sample_sqr(jinc_sqr,arg2,blur2,radius2)*sample_sqr(jinc_sqr,JINC_ZERO_SQR*(arg2/radius2),1.0,radius2));
			break;
		case SP_WT_TRD2 :
			if (arg<(radius/2))
				return(sample_sqr(jinc_sqr,arg2,blur2,radius2));
			else
				return(sample_sqr(jinc_sqr,arg2,blur2,radius2)*((2.0-(2.0*(arg/radius))))); 
			break;
		default : return(0.0); break;
	}
}


// Here is our simple jinc_pi(x), not 2.0x because of auto-normalizing of kernel for convolution in resampling program generator, but M_PI scaled argument
inline static double jinc_pi(double arg)
{
	const auto x = M_PI * arg;
#ifdef C17_ENABLE
	return std::cyl_bessel_j(1, x) / x;
#else
	return bessel_j1(x) / x;
#endif
}

// jinc function value from x,y origin to dx, dy coordinate from origin 
inline static double jpdi(double dx, double dy, int x, int y)
{
	double dist = sqrt((x-dx)*(x-dx)+(y-dy)*(y-dy));
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
	double dist_sq = dx*dx+dy*dy;

	if (dist_sq > radius_sq) return(0.0); // make kernel round, may be option to be square in the future

	auto sum = 0.0;
// 1st row
	sum += jpdi(dx,dy,-1,+2)*k21+jpdi(dx,dy,+0,+2)*k20+jpdi(dx,dy,+1,+2)*k21;
// 2nd row
	sum += jpdi(dx,dy,-2,+1)*k21+jpdi(dx,dy,-1,+1)*k11+jpdi(dx,dy,+0,+1)*k10+jpdi(dx,dy,+1,+1)*k11+jpdi(dx,dy,+2,+1)*k21;
// 3rd row
	sum += jpdi(dx,dy,-2,+0)*k20+jpdi(dx,dy,-1,+0)*k10+jpdi(dx,dy,+0,+0)*1.0+jpdi(dx,dy,+1,+0)*k10+jpdi(dx,dy,+2,+0)*k20;
// 4th row
	sum += jpdi(dx,dy,-2,-1)*k21+jpdi(dx,dy,-1,-1)*k11+jpdi(dx,dy,+0,-1)*k10+jpdi(dx,dy,+1,-1)*k11+jpdi(dx,dy,+2,-1)*k21;
// 5th row
	sum += jpdi(dx,dy,-1,-2)*k21+jpdi(dx,dy,+0,-2)*k20+jpdi(dx,dy,+1,-2)*k21;

	return sum;
}

//static const double DOUBLE_ROUND_MAGIC_NUMBER = 6755399441055744.0;

static bool init_coeff_table(EWAPixelCoeff *out, int quantize_x, int quantize_y,
    int filter_size, int dst_width, int dst_height, int mod_align)
{
    out->filter_size = filter_size;
	if (mod_align > 0)
		out->coeff_stride = (filter_size + (mod_align-1)) & ~(mod_align-1);
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

static void delete_coeff_table(EWAPixelCoeff *out)
{
	myalignedfree(out->factor);
	mydeleteT(out->factor_map);
    mydeleteT(out->meta);
}

struct generate_coeff_params
{
    Lut *func;
    EWAPixelCoeff *out;
    int quantize_x;
    int quantize_y;
    int samples;
    int src_width;
    int src_height;
    int dst_width;
    int dst_height;
    double radius;
    double crop_left;
    double crop_top;
    double crop_width;
    double crop_height;
    int initial_capacity;
    double initial_factor;
	int mod_align;
	bool bUseLUTkernel;
	double blur;
	WEIGHTING_TYPE weighting_type;
	SP_KERNEL_TYPE kernel_type;
	float k10;
	float k20;
	float k11;
	float k21;
};

#ifndef C17_ENABLE
#define llround(x) (x<0.0) ? (long long)floor(x - 0.5) : (long long)floor(x + 0.5)
#define lrintf(x) (long)floor(x + 0.5)
#endif

/* Coefficient table generation */
static bool generate_coeff_table_c(const generate_coeff_params &params)
{
    Lut *func = params.func;
    EWAPixelCoeff *out = params.out;
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
    if (!init_coeff_table(out, quantize_x, quantize_y, filter_size, dst_width, dst_height, mod_align)) return(false);

    size_t tmp_array_capacity = params.initial_capacity;
    float* tmp_array = static_cast<float*>(_aligned_malloc(tmp_array_capacity * sizeof(float), 64));
    if (tmp_array==nullptr) return(false);
    size_t tmp_array_size = 0;
    int tmp_array_top = 0;
    unsigned base_clz = portable_clz(tmp_array_capacity);
    const double initial_growth_factor = params.initial_factor;
    const double radius2 = radius * radius;

    // Use to advance the coeff pointer
    const int coeff_per_pixel = out->coeff_stride * filter_size;

    for (int y = 0; y < dst_height; ++y)
    {
        for (int x = 0; x < dst_width; ++x)
        {
            bool is_border = false;

            EWAPixelCoeffMeta* meta = &out->meta[y * dst_width + x];

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

            meta->start_x = window_begin_x;
            meta->start_y = window_begin_y;

            // Quantize xpos and ypos
            const int quantized_x_int = static_cast<int>(xpos * quantize_x);
            const int quantized_y_int = static_cast<int>(ypos * quantize_y);
            const int quantized_x_value = quantized_x_int % quantize_x;
            const int quantized_y_value = quantized_y_int % quantize_y;
            const float quantized_xpos = static_cast<float>(quantized_x_int) / quantize_x;
            const float quantized_ypos = static_cast<float>(quantized_y_int) / quantize_y;

            if (!is_border && out->factor_map[quantized_y_value * quantize_x + quantized_x_value] != 0)
            {
                // Not border pixel and already have coefficient calculated at this quantized position
                meta->coeff_meta = out->factor_map[quantized_y_value * quantize_x + quantized_x_value] - 1;
            }
            else
            {
                // then need computation
                float divider = 0.f;

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
                const size_t new_size = tmp_array_size + coeff_per_pixel;
                if (new_size > tmp_array_capacity)
                {
                    size_t new_capacity = tmp_array_capacity * (1.0 + (initial_growth_factor - 1.0)
                        * (1.0 - static_cast<double>(max(0, static_cast<int>(base_clz - portable_clz(tmp_array_capacity)))) / 32.0));
                    if (new_capacity < new_size)
                        new_capacity = new_size;
                    float* new_tmp = static_cast<float*>(_aligned_malloc(new_capacity * sizeof(float), 64));
                    if (new_tmp==nullptr)
                    {
                        myalignedfree(tmp_array);
						return(false);
                    }
                    memcpy(new_tmp, tmp_array, tmp_array_size * sizeof(float));
                    myalignedfree(tmp_array);
                    tmp_array = new_tmp;
                    tmp_array_capacity = new_capacity;
                }
                memset(tmp_array + tmp_array_size, 0, coeff_per_pixel * sizeof(float));
                int curr_factor_ptr = tmp_array_top;
                tmp_array_size = new_size;

                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; ++lx)
                    {
                        // Euclidean distance to sampling pixel
                        const double dx = (clamp(is_border ? xpos : quantized_xpos,0.0f,static_cast<float>(src_width-1))-window_x)*filter_step_x;
                        const double dy = (clamp(is_border ? ypos : quantized_ypos,0.0f,static_cast<float>(src_height-1))-window_y)*filter_step_y;

						float factor;
						
						switch(kernel_type)
						{
							case SP_JINCSINGLE :
								if (params.bUseLUTkernel)
								{
									//int index = static_cast<int>(llround((samples-1)*(dx*dx+dy*dy)/radius2 + DOUBLE_ROUND_MAGIC_NUMBER));
									int index = static_cast<int>(llround((samples-1)*(dx*dx+dy*dy)/radius2));
									factor = func->GetFactor(index);
								}
								else
									factor = (float)GetFactor2D(dx,dy,radius,params.blur,params.weighting_type);
								break;
							case SP_JINCSUM :
								factor = (float)GetFactor2D_JINCSUM_21(dx,dy,k10,k20,k11,k21,radius2);
								break;
							default : factor = 0.0; break;
						}

                        tmp_array[curr_factor_ptr + static_cast<int64_t>(lx)] = factor;
                        divider += factor;

                        ++window_x;
                    }

                    curr_factor_ptr += out->coeff_stride;

                    window_x = window_begin_x;
                    ++window_y;
                }

                // Second loop to divide the coeff
                curr_factor_ptr = tmp_array_top;
                for (int ly = 0; ly < filter_size; ++ly)
                {
                    for (int lx = 0; lx < filter_size; ++lx)
                    {
                        tmp_array[curr_factor_ptr + static_cast<int64_t>(lx)] /= divider;
                    }

                    curr_factor_ptr += out->coeff_stride;
                }

                // Save factor to table
                if (!is_border)
                    out->factor_map[quantized_y_value * quantize_x + quantized_x_value] = tmp_array_top + 1;

                meta->coeff_meta = tmp_array_top;
                tmp_array_top += coeff_per_pixel;
            }

            xpos += x_step;
        }

        ypos += y_step;
        xpos = start_x;
    }

    // Copy from tmp_array to real array
    out->factor = tmp_array;
	
	return(true);
}


/* Planar resampling with coeff table */
/* 8-16 bit */
//#pragma intel optimization_parameter target_arch=sse
template<typename T>
static void resize_plane_c_1x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx = 0;

	const T *srcp = reinterpret_cast<const T*>(MT_DataGF->src[idx]);
    T *JincMT_RESTRICT dstp = reinterpret_cast<T*>(MT_DataGF->dst[idx]);

    const ptrdiff_t src_stride = MT_DataGF->src_pitch[idx]/sizeof(T);
    const ptrdiff_t dst_stride = MT_DataGF->dst_pitch[idx]/sizeof(T);

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	EWAPixelCoeffMeta *meta_y = coeff->meta + (Y_Min*dst_width);

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

	const float ValMin = Val_Min[idx];
	const float ValMax = Val_Max[idx];

    for (int y = Y_Min; y < Y_Max; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;

        for (int x = 0; x < dst_width; x++)
        {
			const T *src_ptr = srcp + (meta->start_y * src_stride + meta->start_x);
            const float *coeff_ptr = coeff->factor + meta->coeff_meta;

            float result = 0.0f;

            for (int ly = 0; ly < filter_size; ly++)
            {
                for (int lx = 0; lx < filter_size; lx++)
                {
                    result += src_ptr[lx] * coeff_ptr[lx];
                }
                coeff_ptr += coeff_stride;
                src_ptr += src_stride;
            }

            if JincMT_CONSTEXPR (!(std::is_same<T, float>::value))
                dstp[x] = static_cast<T>(lrintf(clamp(result, ValMin, ValMax)));
            else
                dstp[x] = result;

            meta++;
        }
		
		meta_y += dst_width;
        dstp += dst_stride;
    }
}


/* Planar resampling with coeff table */
/* 8-16 bit */
//#pragma intel optimization_parameter target_arch=sse
template<typename T>
static void resize_plane_c_2x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *coeff,
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

	const int Y_Min = (PlaneYMode) ? MT_DataGF->dst_Y_h_min : MT_DataGF->dst_UV_h_min;
	const int Y_Max = (PlaneYMode) ? MT_DataGF->dst_Y_h_max : MT_DataGF->dst_UV_h_max;
	const int dst_width = (PlaneYMode) ? MT_DataGF->dst_Y_w : MT_DataGF->dst_UV_w;

	const float ValMin1 = Val_Min[idx1];
	const float ValMin2 = Val_Min[idx2];
	const float ValMax1 = Val_Max[idx1];
	const float ValMax2 = Val_Max[idx2];

	EWAPixelCoeffMeta *meta_y = coeff->meta + (Y_Min*dst_width);

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; x++)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const float *coeff_ptr = coeff->factor + meta->coeff_meta;

			float result1 = 0.0f;
			float result2 = 0.0f;

			for (int ly = 0; ly < filter_size; ly++)
			{
				for (int lx = 0; lx < filter_size; lx++)
				{
					const float cf = coeff_ptr[lx];

					result1 += src_ptr1[lx] * cf;
					result2 += src_ptr2[lx] * cf;
				}
				coeff_ptr += coeff_stride;
				src_ptr1 += src_pitch1;
				src_ptr2 += src_pitch2;
			}

			if JincMT_CONSTEXPR(!(std::is_same<T, float>::value))
			{
				dst1[x] = static_cast<T>(lrintf(clamp(result1, ValMin1, ValMax1)));
				dst2[x] = static_cast<T>(lrintf(clamp(result2, ValMin2, ValMax2)));
			}
			else
			{
				dst1[x] = result1;
				dst2[x] = result2;
			}

			meta++;
		}

		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
	}
}


/* Planar resampling with coeff table */
/* 8-16 bit */
//#pragma intel optimization_parameter target_arch=sse
template<typename T>
static void resize_plane_c_3x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *coeff,
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

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const float ValMin1 = Val_Min[idx1];
	const float ValMin2 = Val_Min[idx2];
	const float ValMin3 = Val_Min[idx3];
	const float ValMax1 = Val_Max[idx1];
	const float ValMax2 = Val_Max[idx2];
	const float ValMax3 = Val_Max[idx3];

	EWAPixelCoeffMeta *meta_y = coeff->meta + (Y_Min*dst_width);

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; x++)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const float *coeff_ptr = coeff->factor + meta->coeff_meta;

			float result1 = 0.0f;
			float result2 = 0.0f;
			float result3 = 0.0f;

			for (int ly = 0; ly < filter_size; ly++)
			{
				for (int lx = 0; lx < filter_size; lx++)
				{
					const float cf = coeff_ptr[lx];

					result1 += src_ptr1[lx] * cf;
					result2 += src_ptr2[lx] * cf;
					result3 += src_ptr3[lx] * cf;
				}
				coeff_ptr += coeff_stride;
				src_ptr1 += src_pitch1;
				src_ptr2 += src_pitch2;
				src_ptr3 += src_pitch3;
			}

			if JincMT_CONSTEXPR(!(std::is_same<T, float>::value))
			{
				dst1[x] = static_cast<T>(lrintf(clamp(result1, ValMin1, ValMax1)));
				dst2[x] = static_cast<T>(lrintf(clamp(result2, ValMin2, ValMax2)));
				dst3[x] = static_cast<T>(lrintf(clamp(result3, ValMin3, ValMax3)));
			}
			else
			{
				dst1[x] = result1;
				dst2[x] = result2;
				dst3[x] = result3;
			}

			meta++;
		}

		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
	}
}


/* Planar resampling with coeff table */
/* 8-16 bit */
//#pragma intel optimization_parameter target_arch=sse
template<typename T>
static void resize_plane_c_4x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *coeff,
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

	const float ValMin1 = Val_Min[idx1];
	const float ValMin2 = Val_Min[idx2];
	const float ValMin3 = Val_Min[idx3];
	const float ValMin4 = Val_Min[idx4];
	const float ValMax1 = Val_Max[idx1];
	const float ValMax2 = Val_Max[idx2];
	const float ValMax3 = Val_Max[idx3];
	const float ValMax4 = Val_Max[idx4];

	EWAPixelCoeffMeta *meta_y = coeff->meta + (Y_Min*dst_width);

	const int filter_size = coeff->filter_size, coeff_stride = coeff->coeff_stride;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; x++)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const T *src_ptr4 = src4 + (meta->start_y * src_pitch4 + meta->start_x);
			const float *coeff_ptr = coeff->factor + meta->coeff_meta;

			float result1 = 0.0f;
			float result2 = 0.0f;
			float result3 = 0.0f;
			float result4 = 0.0f;

			for (int ly = 0; ly < filter_size; ly++)
			{
				for (int lx = 0; lx < filter_size; lx++)
				{
					const float cf = coeff_ptr[lx];

					result1 += src_ptr1[lx] * cf;
					result2 += src_ptr2[lx] * cf;
					result3 += src_ptr3[lx] * cf;
					result4 += src_ptr4[lx] * cf;
				}
				coeff_ptr += coeff_stride;
				src_ptr1 += src_pitch1;
				src_ptr2 += src_pitch2;
				src_ptr3 += src_pitch3;
				src_ptr4 += src_pitch4;
			}

			if JincMT_CONSTEXPR(!(std::is_same<T, float>::value))
			{
				dst1[x] = static_cast<T>(lrintf(clamp(result1, ValMin1, ValMax1)));
				dst2[x] = static_cast<T>(lrintf(clamp(result2, ValMin2, ValMax2)));
				dst3[x] = static_cast<T>(lrintf(clamp(result3, ValMin3, ValMax3)));
				dst4[x] = static_cast<T>(lrintf(clamp(result4, ValMin4, ValMax4)));
			}
			else
			{
				dst1[x] = result1;
				dst2[x] = result2;
				dst3[x] = result3;
				dst4[x] = result4;
			}

			meta++;
		}

		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
		dst4 += dst_pitch4;
	}
}


uint8_t CreateMTData(MT_Data_Info_JincResizeMT MT_Data[], uint8_t max_threads, uint8_t threads_number, int32_t src_size_x, int32_t src_size_y,
	int32_t dst_size_x, int32_t dst_size_y, uint8_t shift_src_w_UV, uint8_t shift_src_h_UV, uint8_t shift_dst_w_UV, uint8_t shift_dst_h_UV)
{

	if ((max_threads <= 1) || (max_threads > threads_number))
	{
		MT_Data[0].top = true;
		MT_Data[0].bottom = true;
		MT_Data[0].src_Y_h_min = 0;
		MT_Data[0].dst_Y_h_min = 0;
		MT_Data[0].src_Y_h_max = src_size_y;
		MT_Data[0].dst_Y_h_max = dst_size_y;
		MT_Data[0].src_UV_h_min = 0;
		MT_Data[0].dst_UV_h_min = 0;
		MT_Data[0].src_UV_h_max = src_size_y >> shift_src_h_UV;
		MT_Data[0].dst_UV_h_max = dst_size_y >> shift_dst_h_UV;
		MT_Data[0].src_Y_w = src_size_x;
		MT_Data[0].dst_Y_w = dst_size_x;
		MT_Data[0].src_UV_w = src_size_x >> shift_src_w_UV;
		MT_Data[0].dst_UV_w = dst_size_x >> shift_dst_w_UV;
		return(1);
	}

	int32_t _y_min, _dh;
	int32_t src_dh_Y, src_dh_UV, dst_dh_Y, dst_dh_UV;
	int32_t h_y;
	uint8_t i, max_src = 1, max_dst = 1, max;

	dst_dh_Y = (dst_size_y + (uint32_t)max_threads - 1) / (uint32_t)max_threads;
	if (dst_dh_Y < 16) dst_dh_Y = 16;
	if ((dst_dh_Y & 3) != 0) dst_dh_Y = ((dst_dh_Y + 3) >> 2) << 2;

	if (src_size_y == dst_size_y) src_dh_Y = dst_dh_Y;
	else
	{
		src_dh_Y = (src_size_y + (uint32_t)max_threads - 1) / (uint32_t)max_threads;
		if (src_dh_Y < 16) src_dh_Y = 16;
		if ((src_dh_Y & 3) != 0) src_dh_Y = ((src_dh_Y + 3) >> 2) << 2;
	}

	_y_min = src_size_y;
	_dh = src_dh_Y;
	h_y = _dh;
	while (h_y < (_y_min - 16))
	{
		max_src++;
		h_y += _dh;
	}

	_y_min = dst_size_y;
	_dh = dst_dh_Y;
	h_y = _dh;
	while (h_y < (_y_min - 16))
	{
		max_dst++;
		h_y += _dh;
	}

	//max = (max_src < max_dst) ? max_src : max_dst;
	// Split is made on dst size
	max = max_dst;

	if (max == 1)
	{
		MT_Data[0].top = true;
		MT_Data[0].bottom = true;
		MT_Data[0].src_Y_h_min = 0;
		MT_Data[0].dst_Y_h_min = 0;
		MT_Data[0].src_Y_h_max = src_size_y;
		MT_Data[0].dst_Y_h_max = dst_size_y;
		MT_Data[0].src_UV_h_min = 0;
		MT_Data[0].dst_UV_h_min = 0;
		MT_Data[0].src_UV_h_max = src_size_y >> shift_src_h_UV;
		MT_Data[0].dst_UV_h_max = dst_size_y >> shift_dst_h_UV;
		MT_Data[0].src_Y_w = src_size_x;
		MT_Data[0].dst_Y_w = dst_size_x;
		MT_Data[0].src_UV_w = src_size_x >> shift_src_w_UV;
		MT_Data[0].dst_UV_w = dst_size_x >> shift_dst_w_UV;
		return(1);
	}

	src_dh_UV = src_dh_Y >> shift_src_h_UV;
	dst_dh_UV = dst_dh_Y >> shift_dst_h_UV;

	MT_Data[0].top = true;
	MT_Data[0].bottom = false;
	MT_Data[0].src_Y_h_min = 0;
	MT_Data[0].src_Y_h_max = src_dh_Y;
	MT_Data[0].dst_Y_h_min = 0;
	MT_Data[0].dst_Y_h_max = dst_dh_Y;
	MT_Data[0].src_UV_h_min = 0;
	MT_Data[0].src_UV_h_max = src_dh_UV;
	MT_Data[0].dst_UV_h_min = 0;
	MT_Data[0].dst_UV_h_max = dst_dh_UV;

	i = 1;
	while (i < max)
	{
		MT_Data[i].top = false;
		MT_Data[i].bottom = false;
		MT_Data[i].src_Y_h_min = MT_Data[i - 1].src_Y_h_max;
		MT_Data[i].src_Y_h_max = MT_Data[i].src_Y_h_min + src_dh_Y;
		MT_Data[i].dst_Y_h_min = MT_Data[i - 1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max = MT_Data[i].dst_Y_h_min + dst_dh_Y;
		MT_Data[i].src_UV_h_min = MT_Data[i - 1].src_UV_h_max;
		MT_Data[i].src_UV_h_max = MT_Data[i].src_UV_h_min + src_dh_UV;
		MT_Data[i].dst_UV_h_min = MT_Data[i - 1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max = MT_Data[i].dst_UV_h_min + dst_dh_UV;
		i++;
	}

	MT_Data[max - 1].bottom = true;
	MT_Data[max - 1].src_Y_h_max = src_size_y;
	MT_Data[max - 1].dst_Y_h_max = dst_size_y;
	MT_Data[max - 1].src_UV_h_max = src_size_y >> shift_src_h_UV;
	MT_Data[max - 1].dst_UV_h_max = dst_size_y >> shift_dst_h_UV;

	for (i = 0; i < max; i++)
	{
		MT_Data[i].src_Y_w = src_size_x;
		MT_Data[i].dst_Y_w = dst_size_x;
		MT_Data[i].src_UV_w = src_size_x >> shift_src_w_UV;
		MT_Data[i].dst_UV_w = dst_size_x >> shift_dst_w_UV;
	}

	return(max);
}


void JincResizeMT::FreeData(void)
{
	for (int i=0; i<static_cast<int>(out.size()); ++i)
	{
		if (out[i] != nullptr)
		{
			delete_coeff_table(out[i]);
			mydelete(out[i]);
		}
	}
	
	if (init_lut!=nullptr)
	{
		mydeleteT(init_lut->lut);
		mydelete(init_lut);
	}
}

JincResizeMT::JincResizeMT(PClip _child, int target_width, int target_height, double crop_left, double crop_top, double crop_width, double crop_height,
	int quant_x, int quant_y, int tap, double blur, const char *_cplace, uint8_t _threads, int opt, int initial_capacity, bool initial_capacity_def,
	double initial_factor, int _weighting_type, bool _bUseLUTkernel, SP_KERNEL_TYPE _sp_kernel_type,
	float _k10, float _k20, float _k11, float _k21, float _support,
	int range, bool _sleep, bool negativePrefetch, IScriptEnvironment* env)
    : GenericVideoFilter(_child), init_lut(nullptr),has_at_least_v8(false), has_at_least_v11(false),
	avx512(false), avx2(false), sse41(false), subsampled(false), threads (_threads), sleep(_sleep),
	bUseLUTkernel(_bUseLUTkernel),kernel_type(_sp_kernel_type), k10(_k10), k20(_k20), k11(_k11), k21(_k21),
	support(_support)
{
	UserId = 0;

	Jinc_MT = StaticThreadpool;

	for (int16_t i = 0; i < MAX_MT_THREADS; i++)
	{
		MT_Thread[i].pClass = this;
		MT_Thread[i].f_process = 0;
		MT_Thread[i].thread_Id = (uint8_t)i;
		MT_Thread[i].pFunc = Jinc_MT;
	}

    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; };

    has_at_least_v11 = true;
    try { env->CheckVersion(11); }
    catch (const AvisynthError&) { has_at_least_v11 = false; };

	grey = vi.IsY();
	isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
	isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
	bits_per_pixel = (uint8_t)vi.BitsPerComponent();

    if (!vi.IsPlanar())
        env->ThrowError("JincResizeMT: clip must be in planar format.");

    if (kernel_type == SP_JINCSINGLE)
	{
		if ((tap < 1) || (tap > 16))
			env->ThrowError("JincResizeMT: tap must be between 1..16.");
	}

    if ((quant_x < 1) || (quant_x > 256))
        env->ThrowError("JincResizeMT: quant_x must be between 1..256.");

    if ((quant_y < 1) || (quant_y > 256))
        env->ThrowError("JincResizeMT: quant_y must be between 1..256.");

    if ((opt > 3) || (opt < -1))
        env->ThrowError("JincResizeMT: opt must be between -1..3.");

    if ((blur < 0.0) || (blur > 10.0))
        env->ThrowError("JincResizeMT: blur must be between 0.0..10.0.");
	
	if ((_weighting_type < 0) || (_weighting_type > 2))
		env->ThrowError("JincResizeMT: weighting type must be between 0 and 2");

	// Detection of AVX512 & AVX2 doesn't exist on AVS 2.6, so can't be checked.
	// Just can check at least AVX, if there is not even AVX, there is not AVX512 or AVX2.
	if (has_at_least_v8)
	{
		if ((!(env->GetCPUFlags() & CPUF_AVX512F)) && (opt == 3))
			env->ThrowError("JincResizeMT: opt=3 requires AVX-512F.");
		if ((!(env->GetCPUFlags() & CPUF_AVX2)) && (opt == 2))
			env->ThrowError("JincResizeMT: opt=2 requires AVX2.");
	}
	else
	{
		if ((!(env->GetCPUFlags() & CPUF_AVX)) && (opt == 3))
			env->ThrowError("JincResizeMT: opt=3 requires AVX-512F, there is not even AVX.");
		if ((!(env->GetCPUFlags() & CPUF_AVX)) && (opt == 2))
			env->ThrowError("JincResizeMT: opt=2 requires AVX2, there is not even AVX.");
	}
    if ((!(env->GetCPUFlags() & CPUF_SSE4_1)) && (opt == 1))
        env->ThrowError("JincResizeMT: opt=1 requires SSE4.1.");

	if ((range < 0) || (range > 4))
		env->ThrowError("JincResizeMT: range allowed is [0..4].");

	if (initial_factor < 1.0)
		env->ThrowError("JincResizeMT: initial_factor must be >= 1.0.");

	int src_width = vi.width;
	int src_height = vi.height;

	if (initial_capacity_def)
	{
		if (initial_capacity<=0)
			env->ThrowError("JincResizeMT: initial_capacity must be > 0.");
	}
	else
		initial_capacity = max(target_width * target_height, src_width * src_height);

	if ( vi.Is420() && ( ((target_width%2)!=0) || ((target_height%2)!=0) ) )
		env->ThrowError("JincResizeMT: width and height must be multiple of 2 for 4:2:0 chroma subsampling.");
	if ( vi.Is422() && ((target_width%2)!=0) )
		env->ThrowError("JincResizeMT: width must be multiple of 2 for 4:2:2 chroma subsampling.");
	if (vi.IsYV411() && ((target_width%4)!=0) )
		env->ThrowError("JincResizeMT: width must be multiple of 4 for 4:1:1 chroma subsampling.");

    if (crop_width <= 0.0)
        crop_width = src_width - crop_left + crop_width;

    if (crop_height <= 0.0)
        crop_height = src_height - crop_top + crop_height;
	
	planecount = (uint8_t)vi.NumComponents();
	if ((planecount > 1) && !(vi.Is444() || isRGBPfamily)) subsampled = true;
	
	chroma_placement = AVS_CHROMA_UNUSED;

	if (subsampled)
	{
		// placement explicite parameter like in ConvertToXXX or Text
		// input frame properties, if "auto"
		// When called from ConvertToXXX, chroma is not involved.
		auto frame0 = _child->GetFrame(0, env);
		const AVSMap* props = has_at_least_v11 ? env->getFramePropsRO(frame0) : nullptr;
		chromaloc_parse_merge_with_props(vi, _cplace, props, /* ref*/chroma_placement, AVS_CHROMA_LEFT /*default*/, env);

		if ((chroma_placement != AVS_CHROMA_LEFT) && (chroma_placement != AVS_CHROMA_CENTER)
			&& (chroma_placement != AVS_CHROMA_TOP_LEFT))
			env->ThrowError("JincResizeMT: cplace must be MPEG2, MPEG1, topleft/top_left, auto or empty.");
	}

	double radius;
	
	switch(kernel_type)
	{
		case SP_JINCSINGLE :
			radius = jinc_zeros[tap-1];
			break;
		case SP_JINCSUM :
			radius = support;
			break;
		default : radius = 1.0; break; // some non-zero value
	}
		
	
	const int samples = LUT_SIZE_VALUE;  // should be a multiple of 4

	switch((uint8_t)_weighting_type)
	{
		case 0 : weighting_type = SP_WT_NONE; break;
		case 1 : weighting_type = SP_WT_JINC; break;
		case 2 : weighting_type = SP_WT_TRD2; break;
		default : weighting_type = SP_WT_NONE; break;
	}

	// AVX512 benchmarks are worse than AVX2, so for now, only allow AVX512 on request.
#ifdef AVX512_BUILD_POSSIBLE
	//avx512 = ((!!(env->GetCPUFlags() & CPUF_AVX512F)) && (opt < 0)) || (opt == 3);
	avx512 = (opt == 3);
#endif
#ifdef AVX2_BUILD_POSSIBLE
	avx2 = ((!!(env->GetCPUFlags() & CPUF_AVX2)) && (opt < 0)) || (opt == 2) || avx512;
#endif
	sse41 = ((!!(env->GetCPUFlags() & CPUF_SSE4_1)) && (opt < 0)) || (opt == 1) || avx2 || avx512;

	const int mod_align = (avx512) ? 16 : (avx2) ? 8 : (sse41) ? 4 : 0;

	if (bUseLUTkernel)
	{
		init_lut = new Lut();

		if (init_lut == nullptr)
		{
			FreeData();
			env->ThrowError("JincResizeMT: Error creating lut.");
		}

		if (!init_lut->InitLut(samples, radius, blur, weighting_type))
		{
			FreeData();
			env->ThrowError("JincResizeMT: Error allocating lut.");
		}
	}

    out.emplace_back(new EWAPixelCoeff());
    generate_coeff_params params =
    {
        init_lut,
        out[0],
        quant_x,
        quant_y,
        samples,
        src_width,
        src_height,
        target_width,
        target_height,
        radius,
        crop_left,
        crop_top,
        crop_width,
        crop_height,
        initial_capacity,
        initial_factor,
		mod_align,
		bUseLUTkernel,
		blur,
		weighting_type,
		kernel_type,
		k10,
		k20,
		k11,
		k21
    };

	if (!generate_coeff_table_c(params))
	{
		FreeData();
		env->ThrowError("JincResizeMT: Error generating coeff table [0].");
	}

	const int shift_w = (!grey && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
	const int shift_h = (!grey && !isRGBPfamily) ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0;

	if (subsampled)
	{

        out.emplace_back(new EWAPixelCoeff());
        const double div_w = static_cast<double>(1 << shift_w);
        const double div_h = static_cast<double>(1 << shift_h);

        const double crop_left_uv = ((chroma_placement == AVS_CHROMA_LEFT) || (chroma_placement == AVS_CHROMA_TOP_LEFT)) ?
            (0.5 * (1.0 - static_cast<double>(src_width) / static_cast<double>(target_width)) + crop_left) / div_w : crop_left / div_w;
        const double crop_top_uv = (chroma_placement == AVS_CHROMA_TOP_LEFT) ?
            (0.5 * (1.0 - static_cast<double>(src_height) / static_cast<double>(target_height)) + crop_top) / div_h : crop_top / div_h;

        generate_coeff_params params1 = {
            init_lut,
            out[1],
            quant_x,
            quant_y,
            samples,
            src_width >> shift_w,
            src_height >> shift_h,
            target_width >> shift_w,
            target_height >> shift_h,
            radius,
            crop_left_uv,
            crop_top_uv,
            crop_width / div_w,
            crop_height / div_h,
            initial_capacity / (static_cast<int>(div_w) * static_cast<int>(div_h)),
            initial_factor,
			mod_align,
			bUseLUTkernel,
			blur,
			weighting_type,
			kernel_type,
			k10,
			k20,
			k11,
			k21
        };
		if (!generate_coeff_table_c(params1))
		{
			FreeData();
			env->ThrowError("JincResizeMT: Error generating coeff table [1].");
		}
	}

	uint8_t plane_range[4];

	if ((range != 1) && (range != 4))
	{
		if ((!grey) && !isRGBPfamily)
		{
			plane_range[0] = 2;
			plane_range[1] = 3;
			plane_range[2] = 3;
		}
		else
		{
			if (grey)
			{
				for (unsigned char i = 0; i < 3; i++)
					plane_range[i] = (range == 0) ? 2 : range;
			}
			else
			{
				for (unsigned char i = 0; i < 3; i++)
					plane_range[i] = 1;
			}
		}
	}
	else
	{
		if (isRGBPfamily)
			range = 1;

		for (unsigned char i = 0; i < 3; i++)
			plane_range[i] = range;
	}
	plane_range[3] = 1;

	for (unsigned char i = 0; i < 4; i++)
	{
		if (bits_per_pixel <= 16)
		{
			switch (plane_range[i])
			{
				case 2 :
					ValMin[i] = static_cast<float>(16 << (bits_per_pixel - 8));
					ValMax[i] = static_cast<float>(235 << (bits_per_pixel - 8));
					break;
				case 3 :
					ValMin[i] = static_cast<float>(16 << (bits_per_pixel - 8));
					ValMax[i] = static_cast<float>(240 << (bits_per_pixel - 8));
					break;
				case 4 :
					ValMin[i] = static_cast<float>(16 << (bits_per_pixel - 8));
					ValMax[i] = static_cast<float>((1 << bits_per_pixel) - 1);
					break;
				default:
					ValMin[i] = 0.0f;
					ValMax[i] = static_cast<float>((1 << bits_per_pixel) - 1);
					break;
			}
		}
		else
		{
			if ((!grey) && !isRGBPfamily)
			{
				switch (i)
				{
					case 0 :
					case 3 :
						ValMin[i] = 0.0f;
						ValMax[i] = 1.0f;
						break;
					case 1 :
					case 2 :
						ValMin[i] = -0.5f;
						ValMax[i] = 0.5f;
						break;
					default :
						ValMin[i] = 0.0f;
						ValMax[i] = 1.0f;
						break;
				}
			}
			else
			{
				ValMin[i] = 0.0f;
				ValMax[i] = 1.0f;
			}
		}
	}

	if (target_height < 32) threads_number = 1;
	else threads_number = threads;

    if (vi.ComponentSize() == 1)
    {
#ifdef AVX512_BUILD_POSSIBLE
		if (avx512)
		{
			process_frame_1x = resize_plane_avx512_1x<uint8_t>;
			process_frame_2x = resize_plane_avx512_2x<uint8_t>;
			process_frame_3x = resize_plane_avx512_3x<uint8_t>;
			process_frame_4x = resize_plane_avx512_4x<uint8_t>;
		}
		else
#endif
		{
#ifdef AVX2_BUILD_POSSIBLE
			if (avx2)
			{
				process_frame_1x = resize_plane_avx2_1x<uint8_t>;
				process_frame_2x = resize_plane_avx2_2x<uint8_t>;
				process_frame_3x = resize_plane_avx2_3x<uint8_t>;
				process_frame_4x = resize_plane_avx2_4x<uint8_t>;
			}
			else
#endif
			{
				if (sse41)
				{
					process_frame_1x = resize_plane_sse41_1x<uint8_t>;
					process_frame_2x = resize_plane_sse41_2x<uint8_t>;
					process_frame_3x = resize_plane_sse41_3x<uint8_t>;
					process_frame_4x = resize_plane_sse41_4x<uint8_t>;
				}
				else
				{
					process_frame_1x = resize_plane_c_1x<uint8_t>;
					process_frame_2x = resize_plane_c_2x<uint8_t>;
					process_frame_3x = resize_plane_c_3x<uint8_t>;
					process_frame_4x = resize_plane_c_4x<uint8_t>;
				}
			}
		}
    }
    else if (vi.ComponentSize() == 2)
    {
#ifdef AVX512_BUILD_POSSIBLE
		if (avx512)
		{
			process_frame_1x = resize_plane_avx512_1x<uint16_t>;
			process_frame_2x = resize_plane_avx512_2x<uint16_t>;
			process_frame_3x = resize_plane_avx512_3x<uint16_t>;
			process_frame_4x = resize_plane_avx512_4x<uint16_t>;
		}
		else
#endif
		{
#ifdef AVX2_BUILD_POSSIBLE
			if (avx2)
			{
				process_frame_1x = resize_plane_avx2_1x<uint16_t>;
				process_frame_2x = resize_plane_avx2_2x<uint16_t>;
				process_frame_3x = resize_plane_avx2_3x<uint16_t>;
				process_frame_4x = resize_plane_avx2_4x<uint16_t>;
			}
			else
#endif
			{
				if (sse41)
				{
					process_frame_1x = resize_plane_sse41_1x<uint16_t>;
					process_frame_2x = resize_plane_sse41_2x<uint16_t>;
					process_frame_3x = resize_plane_sse41_3x<uint16_t>;
					process_frame_4x = resize_plane_sse41_4x<uint16_t>;
				}
				else
				{
					process_frame_1x = resize_plane_c_1x<uint16_t>;
					process_frame_2x = resize_plane_c_2x<uint16_t>;
					process_frame_3x = resize_plane_c_3x<uint16_t>;
					process_frame_4x = resize_plane_c_4x<uint16_t>;
				}
			}
		}
    }
    else
    {
#ifdef AVX512_BUILD_POSSIBLE
		if (avx512)
		{
			process_frame_1x = resize_plane_avx512_1x<float>;
			process_frame_2x = resize_plane_avx512_2x<float>;
			process_frame_3x = resize_plane_avx512_3x<float>;
			process_frame_4x = resize_plane_avx512_4x<float>;
		}
		else
#endif
		{
#ifdef AVX2_BUILD_POSSIBLE
			if (avx2)
			{
				process_frame_1x = resize_plane_avx2_1x<float>;
				process_frame_2x = resize_plane_avx2_2x<float>;
				process_frame_3x = resize_plane_avx2_3x<float>;
				process_frame_4x = resize_plane_avx2_4x<float>;
			}
			else
#endif
			{
				if (sse41)
				{
					process_frame_1x = resize_plane_sse41_1x<float>;
					process_frame_2x = resize_plane_sse41_2x<float>;
					process_frame_3x = resize_plane_sse41_3x<float>;
					process_frame_4x = resize_plane_sse41_4x<float>;
				}
				else
				{
					process_frame_1x = resize_plane_c_1x<float>;
					process_frame_2x = resize_plane_c_2x<float>;
					process_frame_3x = resize_plane_c_3x<float>;
					process_frame_4x = resize_plane_c_4x<float>;
				}
			}
		}
    }

	threads_number = CreateMTData(MT_Data,threads_number, threads_number, src_width, src_height, target_width, target_height,
		shift_w, shift_h, shift_w, shift_h);

	if (threads_number > 1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("JincResizeMT: Error with the TheadPool while getting UserId!");
		}
		if (!poolInterface->EnableAllowSeveral(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("JincResizeMT: Error with the TheadPool while allowing multiple request on UserId!");
		}
		if (negativePrefetch)
		{
			if (!poolInterface->DisableWaitonRequest(UserId))
			{
				FreeData();
				poolInterface->DeAllocateAllThreads(true);
				env->ThrowError("JincResizeMT: Error with the TheadPool while disabling wait on request on UserId!");
			}
		}
	}

	vi.width = target_width;
	vi.height = target_height;
}

JincResizeMT::~JincResizeMT()
{
	if (threads_number>1) poolInterface->RemoveUserId(UserId);
	
	FreeData();

	if (threads>1) poolInterface->DeAllocateAllThreads(true);
}

int __stdcall JincResizeMT::SetCacheHints(int cachehints, int frame_range)
{
	switch (cachehints)
	{
	case CACHE_GET_MTMODE:
		return MT_NICE_FILTER;
	default:
		return 0;
	}
}

void JincResizeMT::ProcessFrameMT(MT_Data_Info_JincResizeMT *MT_DataGF, uint8_t IdxFn)
{
	switch (IdxFn)
	{
		case 1 : // Y (Grey)
			process_frame_1x(MT_DataGF, true, out[0], ValMin, ValMax);
			break;
		case 2 : // RGB or YUV not subsampled
			process_frame_3x(MT_DataGF, true, out[0], ValMin, ValMax);
			break;
		case 3 : // RGBA or YUVA not subsampled
			process_frame_4x(MT_DataGF, true, out[0], ValMin, ValMax);
			break;
		case 4 : // YUV subsampled
			process_frame_1x(MT_DataGF, true, out[0], ValMin, ValMax);
			process_frame_2x(MT_DataGF, false, out[1], ValMin, ValMax);
			break;
		case 5 : // YUVA subsampled
			process_frame_2x(MT_DataGF, true, out[0], ValMin, ValMax);
			process_frame_2x(MT_DataGF, false, out[1], ValMin, ValMax);
			break;
		default : break;
	}
}

void JincResizeMT::StaticThreadpool(void *ptr)
{
	Public_MT_Data_Thread *data = (Public_MT_Data_Thread *)ptr;
	JincResizeMT *ptrClass = (JincResizeMT *)data->pClass;
	MT_Data_Info_JincResizeMT *MT_DataGF = ((MT_Data_Info_JincResizeMT *)data->pData) + data->thread_Id;

	ptrClass->ProcessFrameMT(MT_DataGF, data->f_process);
}


PVideoFrame __stdcall JincResizeMT::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = (has_at_least_v8) ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi, 64);

	Public_MT_Data_Thread MT_ThreadGF[MAX_MT_THREADS];
	MT_Data_Info_JincResizeMT MT_DataGF[MAX_MT_THREADS];
	int8_t idxPool = -1;

	const ptrdiff_t src_pitch_1 = (!isRGBPfamily) ? (ptrdiff_t)src->GetPitch(PLANAR_Y) : (ptrdiff_t)src->GetPitch(PLANAR_G);
	const ptrdiff_t dst_pitch_1 = (!isRGBPfamily) ? (ptrdiff_t)dst->GetPitch(PLANAR_Y) : (ptrdiff_t)dst->GetPitch(PLANAR_G);
	const BYTE *srcp_1 = (!isRGBPfamily) ? src->GetReadPtr(PLANAR_Y) : src->GetReadPtr(PLANAR_G);
	BYTE *JincMT_RESTRICT dstp_1 = (!isRGBPfamily) ? dst->GetWritePtr(PLANAR_Y) : dst->GetWritePtr(PLANAR_G);

	const ptrdiff_t src_pitch_2 = (!grey && !isRGBPfamily) ? (ptrdiff_t)src->GetPitch(PLANAR_U) : (isRGBPfamily) ? (ptrdiff_t)src->GetPitch(PLANAR_B) : 0;
	const ptrdiff_t dst_pitch_2 = (!grey && !isRGBPfamily) ? (ptrdiff_t)dst->GetPitch(PLANAR_U) : (isRGBPfamily) ? (ptrdiff_t)dst->GetPitch(PLANAR_B) : 0;
	const BYTE *srcp_2 = (!grey && !isRGBPfamily) ? src->GetReadPtr(PLANAR_U) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_B) : nullptr;
	BYTE *JincMT_RESTRICT dstp_2 = (!grey && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_U) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_B) : nullptr;

	const ptrdiff_t src_pitch_3 = (!grey && !isRGBPfamily) ? (ptrdiff_t)src->GetPitch(PLANAR_V) : (isRGBPfamily) ? (ptrdiff_t)src->GetPitch(PLANAR_R) : 0;
	const ptrdiff_t dst_pitch_3 = (!grey && !isRGBPfamily) ? (ptrdiff_t)dst->GetPitch(PLANAR_V) : (isRGBPfamily) ? (ptrdiff_t)dst->GetPitch(PLANAR_R) : 0;
	const BYTE *srcp_3 = (!grey && !isRGBPfamily) ? src->GetReadPtr(PLANAR_V) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_R) : nullptr;
	BYTE *JincMT_RESTRICT dstp_3 = (!grey && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_V) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_R) : nullptr;

	const ptrdiff_t src_pitch_4 = (isAlphaChannel) ? (ptrdiff_t)src->GetPitch(PLANAR_A) : 0;
	const ptrdiff_t dst_pitch_4 = (isAlphaChannel) ? (ptrdiff_t)dst->GetPitch(PLANAR_A) : 0;
	const BYTE *srcp_4 = (isAlphaChannel) ? src->GetReadPtr(PLANAR_A) : nullptr;
	BYTE *JincMT_RESTRICT dstp_4 = (isAlphaChannel) ? dst->GetWritePtr(PLANAR_A) : nullptr;

	memcpy(MT_ThreadGF, MT_Thread, sizeof(MT_ThreadGF));
	memcpy(MT_DataGF, MT_Data, sizeof(MT_Data));

	for (uint8_t i = 0; i < threads_number; i++)
		MT_ThreadGF[i].pData = (void *)MT_DataGF;

	if (threads_number > 1)
	{
		if ((!poolInterface->RequestThreadPool(UserId, idxPool, threads_number, MT_ThreadGF)) || (idxPool == -1))
			env->ThrowError("JincResizeMT: Error with the TheadPool while requesting threadpool!");
	}

	for (uint8_t i = 0; i < threads_number; i++)
	{
		MT_DataGF[i].src[0] = srcp_1;
		MT_DataGF[i].src[1] = srcp_2;
		MT_DataGF[i].src[2] = srcp_3;
		MT_DataGF[i].src[3] = srcp_4;
		MT_DataGF[i].src_pitch[0] = src_pitch_1;
		MT_DataGF[i].src_pitch[1] = src_pitch_2;
		MT_DataGF[i].src_pitch[2] = src_pitch_3;
		MT_DataGF[i].src_pitch[3] = src_pitch_4;
		MT_DataGF[i].dst[0] = dstp_1 + (MT_Data[i].dst_Y_h_min*dst_pitch_1);
		MT_DataGF[i].dst[1] = dstp_2 + (MT_Data[i].dst_UV_h_min*dst_pitch_2);
		MT_DataGF[i].dst[2] = dstp_3 + (MT_Data[i].dst_UV_h_min*dst_pitch_3);
		MT_DataGF[i].dst[3] = dstp_4 + (MT_Data[i].dst_Y_h_min*dst_pitch_4);
		MT_DataGF[i].dst_pitch[0] = dst_pitch_1;
		MT_DataGF[i].dst_pitch[1] = dst_pitch_2;
		MT_DataGF[i].dst_pitch[2] = dst_pitch_3;
		MT_DataGF[i].dst_pitch[3] = dst_pitch_4;
	}

	uint8_t f_proc = 0;

	if (subsampled)
	{
		if (isAlphaChannel) f_proc = 5;
		else f_proc = 4;
	}
	else
	{
		if (grey) f_proc = 1;
		else
		{
			if (isAlphaChannel) f_proc = 3;
			else f_proc = 2;
		}
	}

	if (threads_number > 1)
	{
		for (uint8_t i = 0; i < threads_number; i++)
			MT_ThreadGF[i].f_process = f_proc;
		if (poolInterface->StartThreads(UserId, idxPool)) poolInterface->WaitThreadsEnd(UserId, idxPool);

		poolInterface->ReleaseThreadPool(UserId, sleep, idxPool);
	}
	else
	{
		switch (f_proc)
		{
			case 1: // Y (Grey)
				process_frame_1x(MT_DataGF, true, out[0], ValMin, ValMax);
				break;
			case 2: // RGB or YUV not subsampled
				process_frame_3x(MT_DataGF, true, out[0], ValMin, ValMax);
				break;
			case 3: // RGBA or YUVA not subsampled
				process_frame_4x(MT_DataGF, true, out[0], ValMin, ValMax);
				break;
			case 4: // YUV subsampled
				process_frame_1x(MT_DataGF, true, out[0], ValMin, ValMax);
				process_frame_2x(MT_DataGF, false, out[1], ValMin, ValMax);
				break;
			case 5: // YUVA subsampled
				process_frame_2x(MT_DataGF, true, out[0], ValMin, ValMax);
				process_frame_2x(MT_DataGF, false, out[1], ValMin, ValMax);
				break;
			default: break;
		}
	}

	if (subsampled && has_at_least_v11)
		env->propSetInt(env->getFramePropsRW(dst), "_ChromaLocation", chroma_placement, 0);

    return dst;
}
