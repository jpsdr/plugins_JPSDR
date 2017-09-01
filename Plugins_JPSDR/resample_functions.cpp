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

#include "resample_functions.h"
#include "./avs/minmax.h"
#include <stdio.h>


/*******************************************
   ***************************************
   **  Helper classes for resample.cpp  **
   ***************************************
 *******************************************/

/***************************
 ***** Point filter *****
 **************************/

double PointFilter::f(double x) {
  return 1.0;
}


/***************************
 ***** Triangle filter *****
 **************************/

double TriangleFilter::f(double x) {
  x = fabs(x);
  return (x<1.0) ? 1.0-x : 0.0;
}





/*********************************
 *** Mitchell-Netravali filter ***
 *********************************/

MitchellNetravaliFilter::MitchellNetravaliFilter (double b=1./3., double c=1./3.) {
  p0 = (   6. -  2.*b            ) / 6.;
  p2 = ( -18. + 12.*b +  6.*c    ) / 6.;
  p3 = (  12. -  9.*b -  6.*c    ) / 6.;
  q0 = (            8.*b + 24.*c ) / 6.;
  q1 = (         - 12.*b - 48.*c ) / 6.;
  q2 = (            6.*b + 30.*c ) / 6.;
  q3 = (      -     b -  6.*c    ) / 6.;
}

double MitchellNetravaliFilter::f (double x) {
  x = fabs(x);
  return (x<1) ? (p0+x*x*(p2+x*p3)) : (x<2) ? (q0+x*(q1+x*(q2+x*q3))) : 0.0;
}


/***********************
 *** Lanczos3 filter ***
 ***********************/
LanczosFilter::LanczosFilter(int t = 3) {
   taps = (double)clamp(t, 1, 100);
}

double LanczosFilter::sinc(double value) {
  if (value > 0.000001) {
    value *= M_PI;
    return sin(value) / value;
  } else {
    return 1.0;
  }
}

double LanczosFilter::f(double value) {
   value = fabs(value);

  if (value < taps) {
    return (sinc(value) * sinc(value / taps));
  } else {
    return 0.0;
  }
}


/***********************
 *** Blackman filter ***
 ***********************/
BlackmanFilter::BlackmanFilter(int t = 4) {
   taps = (double)clamp(t, 1, 100);
   rtaps = 1.0/taps;
}

double BlackmanFilter::f(double value) {
   value = fabs(value);

  if (value < taps) {
    if (value > 0.000001) {
      value *= M_PI;
      return (sin(value) / value) * (0.42 + 0.5*cos(value*rtaps) + 0.08*cos(2*value*rtaps));
    } else {
      return 1.0;
    }
  } else {
    return 0.0;
  }
}


/***********************
 *** Spline16 filter ***
 ***********************/

double Spline16Filter::f(double value) {
  value = fabs(value);

  if (value < 1.0) {
    return ( ( value - 9.0/5.0 ) * value - 1.0/5.0 ) * value + 1.0;
  } else if (value < 2.0) {
    return ( ( -1.0/3.0 * (value-1.0) + 4.0/5.0 ) * (value-1.0) - 7.0/15.0 ) * (value-1.0);
  }
  return 0.0;
}

/***********************
 *** Spline36 filter ***
 ***********************/

double Spline36Filter::f(double value) {
  value = fabs(value);

  if        (value < 1.0) {
    return ( ( 13.0/11.0  * (value    ) - 453.0/ 209.0 ) * (value    ) -   3.0/ 209.0 ) *(value    ) + 1.0;
  } else if (value < 2.0) {
    return ( ( -6.0/11.0  * (value-1.0) + 270.0/ 209.0 ) * (value-1.0) - 156.0/ 209.0 ) *(value-1.0);
  } else if (value < 3.0) {
    return  ( ( 1.0/11.0  * (value-2.0) -  45.0/ 209.0 ) * (value-2.0) +  26.0/ 209.0 ) *(value-2.0);
  }
  return 0.0;
}

/***********************
 *** Spline64 filter ***
 ***********************/

double Spline64Filter::f(double value) {
  value = fabs(value);

  if        (value < 1.0) {
    return (( 49.0/41.0 * (value    ) - 6387.0/2911.0) * (value    ) -    3.0/2911.0) * (value    ) + 1.0;
  } else if (value < 2.0) {
    return ((-24.0/41.0 * (value-1.0) + 4032.0/2911.0) * (value-1.0) - 2328.0/2911.0) * (value-1.0);
  } else if (value < 3.0) {
    return ((  6.0/41.0 * (value-2.0) - 1008.0/2911.0) * (value-2.0) +  582.0/2911.0) * (value-2.0);
  } else if (value < 4.0) {
    return ((- 1.0/41.0 * (value-3.0) +  168.0/2911.0) * (value-3.0) -   97.0/2911.0) * (value-3.0);
  }
  return 0.0;
}

/***********************
 *** Gaussian filter ***
 ***********************/

/* Solve taps from p*value*value < 9 as pow(2.0, -9.0) == 1.0/512.0 i.e 0.5 bit
                     value*value < 9/p       p = param*0.1;
                     value*value < 90/param
                     value*value < 90/{0.1, 22.5, 30.0, 100.0}
                     value*value < {900, 4.0, 3.0, 0.9}
                     value       < {30, 2.0, 1.73, 0.949}         */

GaussianFilter::GaussianFilter(double p = 30.0) {
  param = clamp(p, 0.1, 100.0);
}

double GaussianFilter::f(double value) {
	double p = param*0.1;
	return pow(2.0, - p*value*value);
}

/***********************
 *** Sinc filter ***
 ***********************/
SincFilter::SincFilter(int t = 4) {
   taps = (double)clamp(t, 1, 20);
}

double SincFilter::f(double value) {
   value = fabs(value);

  if (value > 0.000001) {
    value *= M_PI;
    return sin(value)/value;
  } else {
    return 1.0;
  }
}


/******************************
 **** Resampling Patterns  ****
 *****************************/

ResamplingProgram* ResamplingFunction::GetResamplingProgram(int source_size, double crop_start, double crop_size, int target_size, IScriptEnvironment* env)
{
  double filter_scale = double(target_size) / crop_size;
  double filter_step = min(filter_scale, 1.0);
  double filter_support = support() / filter_step;
  int fir_filter_size = int(ceil(filter_support*2));

  ResamplingProgram *program = new ResamplingProgram(fir_filter_size, source_size, target_size, crop_start, crop_size, env);

  // this variable translates such that the image center remains fixed
  double pos;
  double pos_step = crop_size / target_size;

  if (source_size <= filter_support) {
    env->ThrowError("Resize: Source image too small for this resize method. Width=%d, Support=%d",source_size,int(ceil(filter_support)));
  }

  if (fir_filter_size == 1) // PointResize
    pos = crop_start;
  else
    pos = crop_start + ((crop_size - target_size) / (target_size*2)); // TODO this look wrong, gotta check

  for (int i = 0; i < target_size; ++i) {
    // Clamp start and end position such that it does not exceed frame size
    int end_pos = int(pos + filter_support);

    if (end_pos > source_size-1)
      end_pos = source_size-1;

    int start_pos = end_pos - fir_filter_size + 1;

    if (start_pos < 0)
      start_pos = 0;

    program->pixel_offset[i] = start_pos;

    // the following code ensures that the coefficients add to exactly FPScale
    double total = 0.0;

    // Ensure that we have a valid position
    double ok_pos = clamp(pos, 0.0, (double)(source_size-1));

    // Accumulate all coefficients for weighting
    for (int j = 0; j < fir_filter_size; ++j) {
      total += f((start_pos+j - ok_pos) * filter_step);
    }

    if (total == 0.0) {
      // Shouldn't happened for valid positions.
      total = 1.0;
    }

    double value = 0.0;

    // Now we generate real coefficient
    for (int k = 0; k < fir_filter_size; ++k) {
      double new_value = value + f((start_pos+k - ok_pos) * filter_step) / total;
      program->pixel_coefficient[i*fir_filter_size+k] = short(int(new_value*FPScale+0.5) - int(value*FPScale+0.5)); // to make it round across pixels
	  program->pixel_coefficient_float[i*fir_filter_size + k] = float(new_value - value); // no scaling for float
      value = new_value;
    }
    pos += pos_step;
  }
  return program;
}


/******************************
 **** Desampling Patterns  ****
 *****************************/

ResamplingProgram* ResamplingFunction::GetDesamplingProgram(int source_size, double crop_start, double crop_size, int target_size, uint8_t accuracy, IScriptEnvironment* env)
{
  double filter_scale = double(target_size) / crop_size;
  double filter_step = min(filter_scale, 1.0);
  double filter_support = support() / filter_step;
  int fir_filter_size = int(ceil(filter_support*2));

  ResamplingProgram *program = new ResamplingProgram(fir_filter_size, source_size, target_size, crop_start, crop_size, env);

  // this variable translates such that the image center remains fixed
  double pos;
  double pos_step = crop_size / target_size;

  if ((source_size<=filter_support) || (target_size<=filter_support)) {
    env->ThrowError("Resize: Source or target image too small for this resize method. Width=%d,%d, Support=%d",source_size,target_size,int(ceil(filter_support)));
  }

  if (fir_filter_size == 1) // PointResize
    pos = crop_start;
  else
    pos = crop_start + ((crop_size - target_size) / (target_size*2)); // TODO this look wrong, gotta check

  for (int i = 0; i < target_size; ++i) {
    // Clamp start and end position such that it does not exceed frame size
    int end_pos = int(pos + filter_support);

    if (end_pos > source_size-1)
      end_pos = source_size-1;

    int start_pos = end_pos - fir_filter_size + 1;

    if (start_pos < 0)
      start_pos = 0;

    program->pixel_offset[i] = start_pos;

    // the following code ensures that the coefficients add to exactly FPScale
    double total = 0.0;

    // Ensure that we have a valid position
    double ok_pos = clamp(pos, 0.0, (double)(source_size-1));

    // Accumulate all coefficients for weighting
    for (int j = 0; j < fir_filter_size; ++j) {
      total += f((start_pos+j - ok_pos) * filter_step);
    }

    if (total == 0.0) {
      // Shouldn't happened for valid positions.
      total = 1.0;
    }

    double value = 0.0;

    // Now we generate real coefficient
    for (int k = 0; k < fir_filter_size; ++k) {
      double new_value = value + f((start_pos+k - ok_pos) * filter_step) / total;
	  program->pixel_coefficient_float[i*fir_filter_size + k] = float(new_value - value); // no scaling for float
      value = new_value;
    }
    pos += pos_step;
  }

  Matrix_Compute A(target_size,source_size,DATA_FLOAT),B(source_size,source_size,DATA_FLOAT),C(source_size,target_size,DATA_FLOAT);

  if ((env->GetCPUFlags()&CPUF_SSE2)!=0) A.EnableSSE2();
  if ((env->GetCPUFlags()&CPUF_SSE2)!=0) B.EnableSSE2();
  if ((env->GetCPUFlags()&CPUF_SSE2)!=0) C.EnableSSE2();

  if ((env->GetCPUFlags()&CPUF_AVX)!=0) A.EnableAVX();
  if ((env->GetCPUFlags()&CPUF_AVX)!=0) B.EnableAVX();
  if ((env->GetCPUFlags()&CPUF_AVX)!=0) C.EnableAVX();

  A.FillZero();

  for (int i=0; i<target_size; i++)
  {
	  for (int j=0; j<fir_filter_size; j++)
		  A.SetF(i,program->pixel_offset[i]+j,program->pixel_coefficient_float[i*fir_filter_size+j]);
  }

  delete program;

  B.Product_tAA(A);
  B.Inverse();
  C.Product_AtB(B,A);

  switch(accuracy)
  {
	case 0 :
		for (int i=0; i<source_size; i++)
		{
			for (int j=0; j<target_size; j++)
				if (((int16_t)floor(0.5+C.GetF(i,j)*16384))==0) C.SetF(i,j,0.0);
		}
		break;
	case 1 :
		for (int i=0; i<source_size; i++)
		{
			float maxf=0.0;

			for (int j=0; j<target_size; j++)
				if (fabs(C.GetF(i,j))>maxf) maxf=(float)fabs(C.GetF(i,j));

			for (int j=0; j<target_size; j++)
				if (fabs(C.GetF(i,j)/maxf)<0.001) C.SetF(i,j,0.0);
		}
		break;
	case 2 :
		for (int i=0; i<source_size; i++)
		{
			for (int j=0; j<target_size; j++)
				if (fabs(C.GetF(i,j))<1e-6) C.SetF(i,j,0.0);
		}
		break;
	default :
		for (int i=0; i<source_size; i++)
		{
			for (int j=0; j<target_size; j++)
				if (((int16_t)floor(0.5+C.GetF(i,j)*16384))==0) C.SetF(i,j,0.0);
		}
		break;
  }

  fir_filter_size=0;
  for (int i=0; i<source_size; i++)
  {
	  int j1,j2,j=0;

	  while ((j<target_size) && (C.GetF(i,j)==0.0)) j++;
	  j1=j;
	  if (j1<target_size)
	  {
		  j=target_size-1;
		  while ((j>0) && (C.GetF(i,j)==0.0)) j--;
		  j2=j;
	  }
	  if ((j1<target_size) && ((j2-j1+1)>fir_filter_size)) fir_filter_size=j2-j1+1;
  }

  program = new ResamplingProgram(fir_filter_size, target_size, source_size, 0, source_size, env);

  for (int i=0; i<source_size; i++)
  {
	  int start_pos=0;

	   while ((start_pos<target_size) && (C.GetF(i,start_pos)==0.0)) start_pos++;

	  int end_pos = start_pos+(fir_filter_size-1);

	  if (end_pos>=target_size) start_pos=target_size-fir_filter_size;

	  program->pixel_offset[i] = start_pos;

	  for (int j=0; j<fir_filter_size; j++)
	  {
		  program->pixel_coefficient_float[i*fir_filter_size+j]=C.GetF(i,start_pos+j);
		  program->pixel_coefficient[i*fir_filter_size+j] = short(int(C.GetF(i,start_pos+j)*FPScale+0.5));
	  }
  }

  return program;
}
