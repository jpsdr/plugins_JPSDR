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

#ifndef __Resample_Functions_H__
#define __Resample_Functions_H__

#include <malloc.h>
#include <math.h>
#include <vector>
#include <string.h>

#include "./avisynth.h"
#include "./MatrixClass.h"
#include "./avs/alignment.h"

#define myalignedfree(ptr) if (ptr!=nullptr) { _aligned_free(ptr); ptr=nullptr;}

// Original value: 65536
// 2 bits sacrificed because of 16 bit signed MMX multiplication
// NOTE: Don't change this value. It's hard-coded in SIMD code.
const int FPScale8bits = 14; // fixed point scaler 14 bit
const int FPScale = 1 << FPScale8bits; // fixed point scaler (1<<14)
// for 16 bits: one bit less
const int FPScale16bits = 13;
const int FPScale16 = 1 << FPScale16bits; // fixed point scaler for 10-16 bit SIMD signed operation
const int ALIGN_FLOAT_RESIZER_COEFF_SIZE = 8; // simd friendly
// 16: fits even avx512 32-bit float Horizontal
// but not AVX512 uint8 case, where we process 64 pixels at a time
const int ALIGN_RESIZER_TARGET_SIZE = 64;
typedef struct _SafeLimit
{
  bool overread_possible;
  int source_overread_offset;
  int source_overread_beyond_targetx;
} SafeLimit;


struct ResamplingProgram
{
  IScriptEnvironment *Env;
  int source_size, target_size;
  double crop_start, crop_size;
  int filter_size;
  int filter_size_real; // maybe less than filter_size if dimensions are small
  int filter_size_alignment; // for info, 1 (C), 8 (sse or avx2) or 16 (avx2)
  int target_size_alignment; // coeff table exists (and containt zero coeffs) even beyond target_size. Helps alternative H resizers.

  // Array of Integer indicate starting point of sampling
  std::vector<int> pixel_offset;

  int bits_per_pixel;

  // Array of array of coefficient for each pixel
  // {{pixel[0]_coeff}, {pixel[1]_coeff}, ...}
  short *pixel_coefficient;
  float *pixel_coefficient_float;
  // Array of real kernel size, handles edge cases! <= filter_size
  // for SIMD, coefficients are copied over a padded aligned storage
  std::vector<short> kernel_sizes; 
  // 3.7.4- can be different for each line but then they get equalized and aligned.
  
  bool StatusOk;

  size_t cache_size_L2; // in bytes, for possible use in resizers
  int max_scanlines; // recommended vertical stripe for h resamplers

  // In H resizers, when using SIMD loads for speed, these is a "danger zone".
  // If SIMD load from source pixels (src+offset) over reads beyond the allocated
  // source buffer (image scanline width), it can cause access violation.
  // The following SafeLimit structures indicate such cases for various SIMD sizes.
  // and how much to adjust the source pointer to avoid access violation
  // or stop using large SIMD loads after this limit.
  // Depending on the pixel actual size (1, 2 or 4 bytes) and filter size,
  // the overread size can be calculated and the appropriate SafeLimit entry is used
  // by the actual implementation of the resizer.
  SafeLimit safelimit_filter_size_aligned;
  SafeLimit safelimit_4_pixels;
  SafeLimit safelimit_8_pixels;
  SafeLimit safelimit_16_pixels;
  SafeLimit safelimit_32_pixels;
  SafeLimit safelimit_8_pixels_each8th_target;
  SafeLimit safelimit_16_pixels_each16th_target;
  SafeLimit safelimit_64_pixels_each32th_target;
  SafeLimit safelimit_128_pixels_each64th_target;

  int resampler_h_detect_optimal_scanline(int src_width, int tgt_width, size_t l2_cache_size_bytes, size_t pixel_size);
  bool resize_h_planar_gather_permutex_vstripe_check(int iSamplesInTheGroup, int permutex_index_diff_limit, int kernel_size);

  ResamplingProgram(int filter_size, int source_size, int target_size, double crop_start, double crop_size, int bits_per_pixel, IScriptEnvironment* env)
    : Env(env), source_size(source_size), target_size(target_size), crop_start(crop_start), crop_size(crop_size), filter_size(filter_size), filter_size_real(filter_size),
    bits_per_pixel(bits_per_pixel), pixel_coefficient(nullptr), pixel_coefficient_float(nullptr), cache_size_L2(0), max_scanlines(0)
  {
  // In H resizers, when using SIMD loads for speed, these is a "danger zone".
  // If SIMD load from source pixels (src+offset) over reads beyond the allocated
  // source buffer (image scanline width), it can cause access violation.
  // The following SafeLimit structures indicate such cases for various SIMD sizes.
  // and how much to adjust the source pointer to avoid access violation
  // or stop using large SIMD loads after this limit.
  // Depending on the pixel actual size (1, 2 or 4 bytes) and filter size,
  // the overread size can be calculated and the appropriate SafeLimit entry is used
  // by the actual implementation of the resizer.

	safelimit_filter_size_aligned.overread_possible = false;
	safelimit_filter_size_aligned.source_overread_beyond_targetx = -1;
	safelimit_filter_size_aligned.source_overread_offset = -1;

	safelimit_4_pixels = safelimit_filter_size_aligned;
	safelimit_8_pixels = safelimit_filter_size_aligned;
	safelimit_16_pixels = safelimit_filter_size_aligned;
	safelimit_32_pixels = safelimit_filter_size_aligned;
	safelimit_8_pixels_each8th_target = safelimit_filter_size_aligned;
	safelimit_16_pixels_each16th_target = safelimit_filter_size_aligned;
	safelimit_64_pixels_each32th_target = safelimit_filter_size_aligned;
	safelimit_128_pixels_each64th_target = safelimit_filter_size_aligned;

	StatusOk = true;
	//safelimit_filter_size_aligned = { false, -1, -1 };

    filter_size_alignment = 1;
    // align target_size to 8 units to allow safe up to 8 pixels/cycle in H resizers. modded later.
    target_size_alignment = 1;
    // resize_prepare_coeff can override and realign the size of coefficient table

	if (bits_per_pixel<32)
		pixel_coefficient = (short*) _aligned_malloc(sizeof(short)*target_size*filter_size, 64);
	else
		pixel_coefficient_float = (float*) _aligned_malloc(sizeof(float)*target_size*filter_size, 64);

    pixel_offset.resize(target_size);
    kernel_sizes.resize(target_size);

    if (((bits_per_pixel<32) && (pixel_coefficient==nullptr)) || 
		((bits_per_pixel==32) && (pixel_coefficient_float==nullptr)))
	{
	  myalignedfree(pixel_coefficient_float);
	  myalignedfree(pixel_coefficient);
	  StatusOk = false;
      //env->ThrowError("ResamplingProgram: Could not reserve memory.");
    }

	// If AVS+, try to get L2 cache, otherwise set to 0.
	if (env->FunctionExists("ConvertBits"))
	{
		try { cache_size_L2 = env->GetEnvProperty(AEP_CACHESIZE_L2); } catch (const AvisynthError&) { cache_size_L2=0; }
	}
    const int pixel_size_bytes = (bits_per_pixel + 7) / 8;
    max_scanlines = resampler_h_detect_optimal_scanline(source_size, target_size, cache_size_L2, pixel_size_bytes);

	// Set all values to 0
	if (bits_per_pixel<32) memset(pixel_coefficient,0,sizeof(short)*target_size*filter_size);
	else std::fill_n(pixel_coefficient_float, target_size*filter_size, 0.0f);
  };

  ~ResamplingProgram()
  {
	myalignedfree(pixel_coefficient_float);
	myalignedfree(pixel_coefficient);
  };
};

typedef struct ResamplingProgram ResamplingProgram;

void resize_prepare_coeffs(ResamplingProgram* p, IScriptEnvironment* env, int filter_size_alignment);


/*******************************************
   ***************************************
   **  Helper classes for resample.cpp  **
   ***************************************
 *******************************************/


class ResamplingFunction 
/**
  * Pure virtual base class for resampling functions
  */
{
public:
  virtual double f(double x) = 0;
  virtual double support() = 0;

  virtual ResamplingProgram* GetResamplingProgram(int source_size, double crop_start, double crop_size, int target_size, int bits_per_pixel,
	double center_pos_src, double center_pos_dst, IScriptEnvironment* env);
  virtual ResamplingProgram* GetDesamplingProgram(int source_size, double crop_start, double crop_size, int target_size, int bits_per_pixel,
	double center_pos_src, double center_pos_dst, uint8_t accuracy, int SizeY, uint8_t ShiftC, int &SizeOut,IScriptEnvironment* env);
  virtual int GetDesamplingData(int source_size, double crop_start, double crop_size, int target_size, int bits_per_pixel,
  double center_pos_src, double center_pos_dst, uint8_t ShiftC, IScriptEnvironment* env);
};

class PointFilter : public ResamplingFunction 
/**
  * Nearest neighbour (point sampler), used in PointResize
 **/
{
public:
  double f(double x);  
  double support() { return 0.0; }  // 0.0 crashes it.
  // Pre 3.7.4 : 0.0001. Comment: 0.0 crashes it. 
  // 3.7.4- this 0 is specially handled in GetResamplingProgram
};


class TriangleFilter : public ResamplingFunction 
/**
  * Simple triangle filter, used in BilinearResize
 **/
{
public:
  double f(double x);  
  double support() { return 1.0; }
};


class MitchellNetravaliFilter : public ResamplingFunction 
/**
  * Mitchell-Netraveli filter, used in BicubicResize
 **/
{
public:
  MitchellNetravaliFilter(double b=1.0/3.0, double c=1.0/3.0);
  double f(double x);
  double support() { return 2.0; }

private:
  double p0,p2,p3,q0,q1,q2,q3;
};

class LanczosFilter : public ResamplingFunction
/**
  * Lanczos filter, used in LanczosResize
 **/
{
public:
  LanczosFilter(int _taps=3);
	double f(double x);
	double support() { return taps; };

private:
	double sinc(double value);
  double taps;
};

class BlackmanFilter : public ResamplingFunction
/**
  * Blackman filter, used in BlackmanResize
 **/
{
public:
  BlackmanFilter(int _taps=4);
	double f(double x);
	double support() { return taps; };

private:
  double taps, rtaps;
};

// Spline16
class Spline16Filter : public ResamplingFunction
/**
  * Spline16 of Panorama Tools is a cubic-spline, with derivative set to 0 at the edges (4x4 pixels).
 **/
{
public:
	double f(double x);
	double support() { return 2.0; };

private:
};

// Spline36
class Spline36Filter : public ResamplingFunction
/**
  * Spline36 is like Spline16,  except that it uses 6x6=36 pixels.
 **/
{
public:
	double f(double x);
	double support() { return 3.0; };

private:
};

// Spline64
class Spline64Filter : public ResamplingFunction
/**
  * Spline64 is like Spline36,  except that it uses 8x8=64 pixels.
 **/
{
public:
	double f(double x);
	double support() { return 4.0; };

private:
};


class GaussianFilter : public ResamplingFunction
/**
  * GaussianFilter, from swscale.
 **/
{
public:
  GaussianFilter(double p=30.0, double _b=2.0, double _s=4.0);
	double f(double x);
	double support() { return s; }; // <3.7.4 was fixed at 4.0

private:
 double param;
 double b; // base value since 3.7.4
 double s; // variable support since 3.7.4
};

class SincFilter : public ResamplingFunction
/**
  * Sinc filter, used in SincResize
 **/
{
public:
  SincFilter(int _taps=4);
	double f(double x);
	double support() { return taps; };

private:
  double taps;
};

class SinPowerFilter : public ResamplingFunction
// SinPow kernel, used in SinPowResize
{
public:
SinPowerFilter(double p = 2.5);
double f(double x);
double support() { return 2.0; }; // 2 very important, 4 cause bugs

private:
double param;
};

class SincLin2Filter : public ResamplingFunction
	/**
	* SincLin2 filter, used in SincLin2Resize
	**/
{
public:
	SincLin2Filter(int _taps = 15);
	double f(double x);
	double support() { return taps; };

private:
	double sinc(double value);
	double taps;
};

class UserDefined2Filter : public ResamplingFunction
	/**
	  * User-defined by 2 samples filter, used in UDef2Resize
	 **/
{
public:
	UserDefined2Filter(double _b = 121.0, double _c = 19.0, double _s = 2.3);
	double f(double x);
	double support() { return s; }

private:
	double sinc(double value);
	double a, b, c;
	double s; // variable support
};

#endif  // __Reample_Functions_H__
