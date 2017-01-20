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
#include "avisynth.h"

#define myalignedfree(ptr) if (ptr!=NULL) { _aligned_free(ptr); ptr=NULL;}

// Original value: 65536
// 2 bits sacrificed because of 16 bit signed MMX multiplication
// NOTE: Don't change this value. It's hard-coded in SIMD code.
const int FPScale = 16384; // fixed point scaler

// 09-14-2002 - Vlad59 - Lanczos3Resize - Constant added
#define M_PI 3.14159265358979323846

struct ResamplingProgram {
  int source_size, target_size;
  double crop_start, crop_size;
  int filter_size;

  // Array of Integer indicate starting point of sampling
  int* pixel_offset;

  // Array of array of coefficient for each pixel
  // {{pixel[0]_coeff}, {pixel[1]_coeff}, ...}
  short* pixel_coefficient;
  float* pixel_coefficient_float;

  ResamplingProgram(int filter_size, int source_size, int target_size, double crop_start, double crop_size, IScriptEnvironment* env)
    : filter_size(filter_size), source_size(source_size), target_size(target_size), crop_start(crop_start), crop_size(crop_size),
      pixel_offset(NULL), pixel_coefficient(NULL), pixel_coefficient_float(NULL)
  {
    pixel_offset = (int*) _aligned_malloc(sizeof(int) * target_size, 64); // 64-byte alignment
    pixel_coefficient = (short*) _aligned_malloc(sizeof(short) * target_size * filter_size, 64);
	pixel_coefficient_float = (float*) _aligned_malloc(sizeof(float) * target_size * filter_size, 64);
    if ((pixel_offset==NULL) || (pixel_coefficient==NULL) || (pixel_coefficient_float==NULL)) {
	  myalignedfree(pixel_coefficient_float);
	  myalignedfree(pixel_coefficient);
      myalignedfree(pixel_offset);
      env->ThrowError("ResamplingProgram: Could not reserve memory.");
    }
	
  };

  ~ResamplingProgram() {
	myalignedfree(pixel_coefficient_float);
	myalignedfree(pixel_coefficient);
    myalignedfree(pixel_offset);
  };
};

typedef struct ResamplingProgram ResamplingProgram;


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

  virtual ResamplingProgram* GetResamplingProgram(int source_size, double crop_start, double crop_size, int target_size, IScriptEnvironment* env);
};

class PointFilter : public ResamplingFunction 
/**
  * Nearest neighbour (point sampler), used in PointResize
 **/
{
public:
  double f(double x);  
  double support() { return 0.0001; }  // 0.0 crashes it.
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
  MitchellNetravaliFilter(double b, double c);
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
  LanczosFilter(int taps);
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
  BlackmanFilter(int taps);
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
  GaussianFilter(double p);
	double f(double x);
	double support() { return 4.0; };

private:
 double param;
};

class SincFilter : public ResamplingFunction
/**
  * Sinc filter, used in SincResize
 **/
{
public:
  SincFilter(int taps);
	double f(double x);
	double support() { return taps; };

private:
  double taps;
};


#endif  // __Reample_Functions_H__
