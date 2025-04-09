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

#include "./resample.h"

#define myfree(ptr) if (ptr!=NULL) { free(ptr); ptr=NULL;}
#define mydelete(ptr) if (ptr!=NULL) { delete ptr; ptr=NULL;}
#define mydelete2(ptr) if (ptr!=NULL) { delete[] ptr; ptr=NULL;}

#include <type_traits>

#include "./resample_sse.h"

#if _MSC_VER >= 1900
  #define AVX2_BUILD_POSSIBLE
#endif

#ifdef AVX2_BUILD_POSSIBLE
#include "./resample_avx2.h"
#endif

extern ThreadPoolInterface *poolInterface;


static bool is_paramstring_empty_or_auto(const char* param)
{
	if (param==nullptr) return true;
	return (strcoll(param,"auto")==0); // true is match
}

static bool getChromaLocation(const char* chromaloc_name, IScriptEnvironment* env, ChromaLocation_e& _ChromaLocation)
{
	ChromaLocation_e index=AVS_CHROMA_UNUSED;

	if (strcoll(chromaloc_name,"left")==0) index=AVS_CHROMA_LEFT;
	if (strcoll(chromaloc_name,"center")==0) index=AVS_CHROMA_CENTER;
	if (strcoll(chromaloc_name,"top_left")==0) index=AVS_CHROMA_TOP_LEFT;
	if (strcoll(chromaloc_name,"top")==0) index=AVS_CHROMA_TOP; // not used in Avisynth
	if (strcoll(chromaloc_name,"bottom_left")==0) index=AVS_CHROMA_BOTTOM_LEFT; // not used in Avisynth
	if (strcoll(chromaloc_name,"bottom")==0) index=AVS_CHROMA_BOTTOM; // not used in Avisynth
	if (strcoll(chromaloc_name,"dv")==0) index=AVS_CHROMA_DV; // Special to Avisynth
	// compatibility
	if (strcoll(chromaloc_name,"mpeg1")==0) index=AVS_CHROMA_CENTER;
	if (strcoll(chromaloc_name,"mpeg2")==0) index=AVS_CHROMA_LEFT;
	if (strcoll(chromaloc_name,"jpeg")==0) index=AVS_CHROMA_CENTER;

	if (index!=AVS_CHROMA_UNUSED)
	{
		_ChromaLocation = index;
		return true;
	}

	env->ThrowError("Unknown chroma placement");
	// empty
	return false;
}

static void chromaloc_parse_merge_with_props(const VideoInfo& vi, const char* chromaloc_name, const AVSMap* props, ChromaLocation_e& _ChromaLocation, ChromaLocation_e _ChromaLocation_Default, IScriptEnvironment* env)
{
	if (props!=nullptr)
	{
		if (vi.Is420() || vi.Is422() || vi.IsYV411())
		{ // yes, YV411 can also have valid _ChromaLocation, if 'left'-ish one is given
			if (env->propNumElements(props,"_ChromaLocation")>0)
				_ChromaLocation_Default = (ChromaLocation_e)env->propGetIntSaturated(props,"_ChromaLocation",0,nullptr);
		}
		else
		{
			// Theoretically RGB and not subsampled formats must not have chroma location
			if (env->propNumElements(props,"_ChromaLocation")>0)
			{
				// Uncommented for a while, just ignore when there is any
				// env->ThrowError("Error: _ChromaLocation property found at a non-subsampled source.");
			}
		}
	}

	if (is_paramstring_empty_or_auto(chromaloc_name) || !getChromaLocation(chromaloc_name, env, _ChromaLocation))
		_ChromaLocation = _ChromaLocation_Default;
}

// Borrowed from fmtconv
// ChromaPlacement.cpp
// Author : Laurent de Soras, 2015

// Fixes the vertical chroma placement when the picture is interlaced.
// ofs = ordinate to skip between TFF and BFF, relative to the chroma grid. A
// single line of full-res picture is 0.25.
static inline void ChromaPlacement_fix_itl(double& cp_v, bool interlaced_flag, bool top_flag, double ofs = 0.5)
{
  assert(cp_v >= 0);

  if (interlaced_flag)
  {
    cp_v *= 0.5;
    if (!top_flag) cp_v += ofs;
  }
}
/*
ss_h and ss_v are log2(subsampling)
rgb_flag actually means that chroma subsampling doesn't apply.

http://www.mir.com/DMG/chroma.html

cp_* is the position of the sampling point relative to the frame
top/left border, in the plane coordinates. For reference, the border
of the frame is at 0.5 units of luma from the first luma sampling point.
I. e., the luma sampling point is at the pixel's center.
*/

// PF added BOTTOM, BOTTOM_LEFT, TOP
// Pass ChromaLocation_e::AVS_CHROMA_UNUSED for defaults
// plane index 0:Y, 1:U, 2:V
// cplace is a ChromaLocation_e constant
static void ChromaPlacement_compute_cplace(double& cp_h, double& cp_v, ChromaLocation_e cplace, int plane_index, int ss_h, int ss_v, bool rgb_flag, bool interlaced_flag, bool top_flag)
{
  assert(cplace >= 0 || cplace == AVS_CHROMA_UNUSED);
  assert(cplace < AVS_CHROMA_DV);
  assert(ss_h >= 0);
  assert(ss_v >= 0);
  assert(plane_index >= 0);

  // Generic case for luma, non-subsampled chroma and center (MPEG-1) chroma.
  cp_h = 0.5;
  cp_v = 0.5;
  ChromaPlacement_fix_itl(cp_v, interlaced_flag, top_flag);

  // Subsampled chroma
  if (!rgb_flag && plane_index > 0)
  {
    if (ss_h > 0) // horizontal subsampling 420 411
    {
      if (cplace == AVS_CHROMA_LEFT // mpeg2
        || cplace == AVS_CHROMA_DV
        || cplace == AVS_CHROMA_TOP_LEFT
        || cplace == AVS_CHROMA_BOTTOM_LEFT
        )
      {
        cp_h = 0.5 / (1 << ss_h);
      }
    }

    if (ss_v == 1) // vertical subsampling 420, 422
    {
      if (cplace == AVS_CHROMA_LEFT)
      {
        cp_v = 0.5;
        ChromaPlacement_fix_itl(cp_v, interlaced_flag, top_flag);
      }
      else if (cplace == AVS_CHROMA_DV
        || cplace == AVS_CHROMA_TOP_LEFT
        || cplace == AVS_CHROMA_TOP
        )
      {
        cp_v = 0.25;
        ChromaPlacement_fix_itl(cp_v, interlaced_flag, top_flag, 0.25);

        if (cplace == AVS_CHROMA_DV && plane_index == 2) // V
        {
          cp_v += 0.5;
        }
      }
      else if (cplace == AVS_CHROMA_BOTTOM_LEFT
        || cplace == AVS_CHROMA_BOTTOM
        )
      {
        cp_v = 0.75;
        ChromaPlacement_fix_itl(cp_v, interlaced_flag, top_flag, 0.25);
      }
    }  // ss_v == 1
  }
}

// returns the requested horizontal or vertical pixel center position
static void GetCenterShiftForResizers(double& center_pos_luma, double& center_pos_chroma, bool preserve_center, ChromaLocation_e chroma_placement, VideoInfo &vi, bool for_horizontal)
{
  double center_pos_h_luma = 0;
  double center_pos_v_luma = 0;
  // if not needed, these won't be used
  double center_pos_h_chroma = 0;
  double center_pos_v_chroma = 0;

  // chroma, only if applicable
  if (vi.IsPlanar() && vi.NumComponents() > 1 && !vi.IsRGB())
  {
    double cp_s_h = 0;
    double cp_s_v = 0;

    if (preserve_center)
	{
      // same for source and destination
      int plane_index = 1; // U
      int src_ss_h = vi.GetPlaneWidthSubsampling(PLANAR_U);
      int src_ss_v = vi.GetPlaneHeightSubsampling(PLANAR_U);

      ChromaLocation_e chromaplace = AVS_CHROMA_CENTER; // MPEG1

      ChromaPlacement_compute_cplace(
        cp_s_h, cp_s_v, chroma_placement, plane_index, src_ss_h, src_ss_v,
        vi.IsRGB(),
        false, // interlacing flag, we don't handle it here
        false  // top_flag, we don't handle it here
      );
    }

    center_pos_h_chroma = cp_s_h;
    center_pos_v_chroma = cp_s_v;
  }

  // luma/rgb planes
  if (preserve_center)
  {
    center_pos_h_luma = 0.5;
    center_pos_v_luma = 0.5;
  }
  else
  {
    center_pos_h_luma = 0.0;
    center_pos_v_luma = 0.0;
  }

  // fill return ref values
  if (for_horizontal)
  {
    center_pos_luma = center_pos_h_luma;
    center_pos_chroma = center_pos_h_chroma;
  }
  else
  {
    // vertical
    center_pos_luma = center_pos_v_luma;
    center_pos_chroma = center_pos_v_chroma;
  }
}

/***************************************
 ***** Vertical Resizer Assembly *******
 ***************************************/

template<typename pixel_t>
static void resize_v_planar_pointresize(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  pixel_t *src0 = (pixel_t *)src;
  pixel_t *dst0 = (pixel_t *)dst;
  dst_pitch/=sizeof(pixel_t);

  for (int y = MinY; y < MaxY; y++)
  {
	const pixel_t *src_ptr = src0 + pitch_table[program->pixel_offset[y]];
    
	memcpy(dst0,src_ptr,width*sizeof(pixel_t));

    dst0+=dst_pitch;
  }
}


static void resize_v_c_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
	const int filter_size = program->filter_size;
	const int kernel_size = program->filter_size_real;
	const short *current_coeff = program->pixel_coefficient+filter_size*MinY;

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int Offset = 1 << (FPScale8bits-1);

	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		const int TabMax[4] = {235,240,235,240};

		for (int y = MinY; y < MaxY; y++)
		{
			const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

			for (int x = 0; x < width; x++)
			{
				int result = 0;

				for (int i = 0; i < kernel_size; i++)
					result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

				result = (result+Offset) >> FPScale8bits;
				result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x3] : (result<16) ? 16 : result;
				dst[x] = (BYTE) result;
			}

			dst += dst_pitch;
			current_coeff += filter_size;
		}
	}
	else
	{
		for (int y = MinY; y < MaxY; y++)
		{
			const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

			for (int x = 0; x < width; x++)
			{
				int result = 0;

				for (int i = 0; i < kernel_size; i++)
					result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

				result = (result+Offset) >> FPScale8bits;
				result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
				dst[x] = (BYTE) result;
			}

			dst += dst_pitch;
			current_coeff += filter_size;
		}
	}
}


static void resize_v_c_planar_f(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const int kernel_size = program->filter_size_real;
  const float *current_coeff = program->pixel_coefficient_float+filter_size*MinY;

  const float *src0 = (float *)src;
  float *dst0 = (float *)dst;

  dst_pitch>>=2;

  for (int y = MinY; y < MaxY; y++)
  {
	const float *src_ptr = src0 + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < width; x++)
	{
      float result = 0;

      for (int i = 0; i < kernel_size; i++)
		result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

      dst0[x] = result;
    }

    dst0 += dst_pitch;
    current_coeff += filter_size;
  }
}


static void resize_v_c_planar_s(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
	const int filter_size = program->filter_size;
	const int kernel_size = program->filter_size_real;
	const short *current_coeff = program->pixel_coefficient+filter_size*MinY;

	const uint16_t *src0 = (uint16_t *)src;
	uint16_t *dst0 = (uint16_t *)dst;
	const __int64 val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const __int64 val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
		((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));

	const __int64 Offset = 1 << (FPScale16bits-1);

	dst_pitch>>=1;

	for (int y = MinY; y < MaxY; y++)
	{
		const uint16_t *src_ptr = src0 + pitch_table[program->pixel_offset[y]];

		for (int x = 0; x < width; x++)
		{
			__int64 result = 0;

			for (int i = 0; i < kernel_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

			result = (result+Offset) >> FPScale16bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst0[x] = (uint16_t) result;
		}
		dst0 += dst_pitch;
		current_coeff += filter_size;
	}
}

__forceinline static void resize_v_create_pitch_table(int* table, int pitch, int height, uint8_t pixel_size)
{
  switch(pixel_size)
  {
	case 2 : pitch>>=1; break;
	case 4 : pitch>>=2; break;
	default : ;
  }
  for (int i=0; i<height; i++)
    table[i]=i*pitch;
}


/***************************************
 ********* Horizontal Resizer** ********
 ***************************************/

static void resize_h_c_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const int kernel_size = program->filter_size_real;
  int y_src_pitch=0,y_dst_pitch=0;
  
	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int Offset = 1 << (FPScale8bits-1);

  // external loop y is much faster

	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		const int TabMax[4] = {235,240,235,240};

		for (int y = 0; y < height; y++)
		{
			const short *current_coeff=program->pixel_coefficient;
	  
			for (int x = 0; x < width; x++)
			{
				const int begin = program->pixel_offset[x];
				int result = 0;
		  
				for (int i = 0; i < kernel_size; i++)
	    			result+=(src+y_src_pitch)[(begin+i)]*current_coeff[i];
		
				result = (result + Offset) >> FPScale8bits;
				result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<16) ? 16 : result;
				(dst + y_dst_pitch)[x] = (BYTE)result;		  		  
				current_coeff+=filter_size;
			}
			y_dst_pitch+=dst_pitch;
			y_src_pitch+=src_pitch;	  
		}
	}
	else
	{
		for (int y = 0; y < height; y++)
		{
			const short *current_coeff=program->pixel_coefficient;
	  
			for (int x = 0; x < width; x++)
			{
				const int begin = program->pixel_offset[x];
				int result = 0;
		  
				for (int i = 0; i < kernel_size; i++)
	    			result+=(src+y_src_pitch)[(begin+i)]*current_coeff[i];
		
				result = (result + Offset) >> FPScale8bits;
				result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
				(dst + y_dst_pitch)[x] = (BYTE)result;		  		  
				current_coeff+=filter_size;
			}
			y_dst_pitch+=dst_pitch;
			y_src_pitch+=src_pitch;	  
		}
	}
 
}


static void resize_h_c_planar_s(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const int kernel_size = program->filter_size_real;
  int y_src_pitch=0,y_dst_pitch=0;
  const uint16_t *src0 = (uint16_t *)src;
  uint16_t *dst0 = (uint16_t *)dst;
	const __int64 val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const __int64 val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
		((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));

	const __int64 Offset = 1 << (FPScale16bits-1);


  src_pitch>>=1;
  dst_pitch>>=1;
  
	for (int y = 0; y < height; y++)
	{
		const short *current_coeff=program->pixel_coefficient;
  
		for (int x = 0; x < width; x++)
		{
			const int begin = program->pixel_offset[x];
			__int64 result = 0;
		  
			for (int i = 0; i < kernel_size; i++)
				result+=(src0+y_src_pitch)[(begin+i)]*current_coeff[i];
		  
			result = (result + Offset) >> FPScale16bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			(dst0 + y_dst_pitch)[x] = (uint16_t)result;
			current_coeff+=filter_size;
		}
		y_dst_pitch+=dst_pitch;
		y_src_pitch+=src_pitch;
	}
}


static void resize_h_c_planar_f(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const int kernel_size = program->filter_size_real;
  int y_src_pitch=0,y_dst_pitch=0;
  const float *src0=(float *)src;
  float *dst0=(float *)dst;

  src_pitch>>=2;
  dst_pitch>>=2;

  for (int y = 0; y < height; y++)
  {
	  const float *current_coeff=program->pixel_coefficient_float;
	  
	  for (int x = 0; x < width; x++)
	  {
		  const int begin = program->pixel_offset[x];
		  float result = 0;
		  
		  for (int i = 0; i < kernel_size; i++)
			  result+=(src0+y_src_pitch)[(begin+i)]*current_coeff[i];
		  
		  (dst0 + y_dst_pitch)[x] = result;
		  current_coeff+=filter_size;
	  }
	  y_dst_pitch+=dst_pitch;
	  y_src_pitch+=src_pitch;
  }
}

/********************************************************************
***** Declare index of new filters for Avisynth's filter engine *****
********************************************************************/


FilteredResizeH::FilteredResizeH( PClip _child, double subrange_left, double subrange_width,int target_width,
	uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,bool negativePrefetch,
	bool _avsp,bool preserve_center,ChromaLocation_e chroma_placement,ResamplingFunction* func,IScriptEnvironment* env )
  : GenericVideoFilter(_child),
  resampling_program_luma(NULL), resampling_program_chroma(NULL),
  filter_storage_luma(NULL), filter_storage_chroma(NULL),threads(_threads),sleep(_sleep),
  avsp(_avsp)
{
  src_height = vi.height;
  dst_height = vi.height;
  
  pixelsize = (uint8_t)vi.ComponentSize(); // AVS16
  grey = vi.IsY();
  bits_per_pixel = (uint8_t)vi.BitsPerComponent();
  isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
  isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
  mode_YUY2 = vi.IsYUY2();

  Enable_MMX = (env->GetCPUFlags() & CPUF_MMX)!=0;
  Enable_SSE2 = (env->GetCPUFlags() & CPUF_SSE2)!=0;
  Enable_SSE3 = (env->GetCPUFlags() & CPUF_SSE3)!=0;
  Enable_SSSE3 = (env->GetCPUFlags() & CPUF_SSSE3)!=0;
  Enable_SSE4_1 = (env->GetCPUFlags() & CPUF_SSE4_1)!=0;
  Enable_AVX2 = false;
#ifdef AVX2_BUILD_POSSIBLE
  Enable_AVX2 = avsp && ((env->GetCPUFlags() & CPUF_AVX2)!=0);
#endif

  double center_pos_h_luma;
  double center_pos_h_chroma;
  GetCenterShiftForResizers(center_pos_h_luma, center_pos_h_chroma, preserve_center, chroma_placement, vi, true); // True for horizontal

  src_width = vi.IsPlanar() ? vi.width : vi.BytesFromPixels(vi.width)/pixelsize;
  dst_width = vi.IsPlanar() ? target_width : vi.BytesFromPixels(target_width)/pixelsize;

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;

	int16_t i;

	ResampleH_MT=StaticThreadpoolH;

	for (i=0; i<MAX_MT_THREADS; i++)
	{
		MT_Thread[i].pClass=this;
		MT_Thread[i].f_process=0;
		MT_Thread[i].thread_Id=(uint8_t)i;
		MT_Thread[i].pFunc=ResampleH_MT;
	}

	UserId=0;
	
	const int shift_w = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
	const int shift_h = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0;

	if (vi.height<32) threads_number=1;
	else threads_number=threads;

  // Main resampling program
  int SizeH;

  if (desample) resampling_program_luma = func->GetDesamplingProgram(target_width, subrange_left, subrange_width, vi.width, bits_per_pixel,
	  center_pos_h_luma, center_pos_h_luma, // for resizing it's the same for source and dest
	  accuracy, 0, shift_w, SizeH, env);
  else
  {
	  resampling_program_luma = func->GetResamplingProgram(vi.width, subrange_left, subrange_width, target_width, bits_per_pixel,
		  center_pos_h_luma, center_pos_h_luma, // for resizing it's the same for source and dest
		  env);
	  SizeH=dst_width;
  }
  
  if (resampling_program_luma==NULL)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  if (desample) env->ThrowError("ResizeHMT: Error while GetDesamplingProgram for luma!");
	  else env->ThrowError("ResizeHMT: Error while GetResamplingProgram for luma!");
  }

  if (desample && ((SizeH>vi.width) || (SizeH==-1)))
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  if (SizeH>vi.width) env->ThrowError("ResizeHMT: Desampling can only downscale!");
	  else env->ThrowError("ResizeHMT: Matrix can't be reversed!");
  }

  if (vi.IsPlanar() && !grey && !isRGBPfamily)
  {
    const int div   = 1 << shift_w;

	if (desample)
	{
		int SizeOut;

	    resampling_program_chroma = func->GetDesamplingProgram(
		  target_width   >> shift_w,
	      subrange_left   / div,
		  subrange_width  / div,
	      vi.width   >> shift_w,
		  bits_per_pixel,
		  center_pos_h_chroma, center_pos_h_chroma, // horizontal
		  accuracy,SizeH,shift_w,SizeOut,
		  env);
		if (SizeOut==-1)
		{
			FreeData();
			if (threads>1) poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeHMT: Matrix can't be reversed!");
		}
	}
	else
	{
	    resampling_program_chroma = func->GetResamplingProgram(
		  vi.width       >> shift_w,
	      subrange_left   / div,
		  subrange_width  / div,
	      target_width   >> shift_w,
		  bits_per_pixel,
		  center_pos_h_chroma, center_pos_h_chroma, // horizontal
		  env);
	}

	if (resampling_program_chroma==NULL)
	{
		FreeData();
		if (threads>1) poolInterface->DeAllocateAllThreads(true);
		if (desample) env->ThrowError("ResizeHMT: Error while GetDesamplingProgram for chroma!");
		else env->ThrowError("ResizeHMT: Error while GetResamplingProgram for chroma!");
	}

	const int w_UV=vi.width >> shift_w;
	
	resampler_h_chroma = GetResampler(true,resampling_program_chroma,env);
  }
  
  resampler_h_luma = GetResampler(true,resampling_program_luma,env);
  
  threads_number=CreateMTData(threads_number,src_width,vi.height,SizeH,vi.height,shift_w,shift_h);

	if (threads_number>1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeHMT: Error with the TheadPool while getting UserId!");
		}
		if (!poolInterface->EnableAllowSeveral(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeHMT: Error with the TheadPool while allowing multiple request on UserId!");
		}
		if (negativePrefetch)
		{
			if (!poolInterface->DisableWaitonRequest(UserId))
			{
				FreeData();
				poolInterface->DeAllocateAllThreads(true);
				env->ThrowError("ResizeHMT: Error with the TheadPool while disabling wait on request on UserId!");
			}
		}
	}

	has_at_least_v8=true;
	try { env->CheckVersion(8); } catch (const AvisynthError&) { has_at_least_v8=false; }

  // Change target video info size
  vi.width =SizeH;
}


int __stdcall FilteredResizeH::SetCacheHints(int cachehints,int frame_range)
{
  switch (cachehints)
  {
  case CACHE_GET_MTMODE :
    return MT_NICE_FILTER;
  default :
    return 0;
  }
}


uint8_t FilteredResizeH::CreateMTData(uint8_t max_threads,int32_t src_size_x,int32_t src_size_y,int32_t dst_size_x,int32_t dst_size_y,int UV_w,int UV_h)
{
	if ((max_threads<=1) || (max_threads>threads_number))
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	int32_t _y_min,_dh;
	int32_t src_dh_Y,src_dh_UV,dst_dh_Y,dst_dh_UV;
	int32_t h_y;
	uint8_t i,max_src=1,max_dst=1,max;

	dst_dh_Y=(dst_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
	if (dst_dh_Y<16) dst_dh_Y=16;
	if ((dst_dh_Y & 3)!=0) dst_dh_Y=((dst_dh_Y+3) >> 2) << 2;

	if (src_size_y==dst_size_y) src_dh_Y=dst_dh_Y;
	else
	{
		src_dh_Y=(src_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
		if (src_dh_Y<16) src_dh_Y=16;
		if ((src_dh_Y & 3)!=0) src_dh_Y=((src_dh_Y+3) >> 2) << 2;
	}

	_y_min=src_size_y;
	_dh=src_dh_Y;
	h_y=_dh;
	while (h_y<(_y_min-16))
	{
		max_src++;
		h_y+=_dh;
	}

	_y_min=dst_size_y;
	_dh=dst_dh_Y;
	h_y=_dh;
	while (h_y<(_y_min-16))
	{
		max_dst++;
		h_y+=_dh;
	}

	max=(max_src<max_dst) ? max_src:max_dst;

	if (max==1)
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	src_dh_UV= (UV_h>0) ? src_dh_Y>>UV_h : src_dh_Y;
	dst_dh_UV= (UV_h>0) ? dst_dh_Y>>UV_h : dst_dh_Y;

	MT_Data[0].top=true;
	MT_Data[0].bottom=false;
	MT_Data[0].src_Y_h_min=0;
	MT_Data[0].src_Y_h_max=src_dh_Y;
	MT_Data[0].dst_Y_h_min=0;
	MT_Data[0].dst_Y_h_max=dst_dh_Y;
	MT_Data[0].src_UV_h_min=0;
	MT_Data[0].src_UV_h_max=src_dh_UV;
	MT_Data[0].dst_UV_h_min=0;
	MT_Data[0].dst_UV_h_max=dst_dh_UV;

	i=1;
	while (i<max)
	{
		MT_Data[i].top=false;
		MT_Data[i].bottom=false;
		MT_Data[i].src_Y_h_min=MT_Data[i-1].src_Y_h_max;
		MT_Data[i].src_Y_h_max=MT_Data[i].src_Y_h_min+src_dh_Y;
		MT_Data[i].dst_Y_h_min=MT_Data[i-1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max=MT_Data[i].dst_Y_h_min+dst_dh_Y;
		MT_Data[i].src_UV_h_min=MT_Data[i-1].src_UV_h_max;
		MT_Data[i].src_UV_h_max=MT_Data[i].src_UV_h_min+src_dh_UV;
		MT_Data[i].dst_UV_h_min=MT_Data[i-1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max=MT_Data[i].dst_UV_h_min+dst_dh_UV;
		i++;
	}

	MT_Data[max-1].bottom=true;
	MT_Data[max-1].src_Y_h_max=src_size_y;
	MT_Data[max-1].dst_Y_h_max=dst_size_y;
	if (UV_h>0)
	{
		MT_Data[max-1].src_UV_h_max=src_size_y >> UV_h;
		MT_Data[max-1].dst_UV_h_max=dst_size_y >> UV_h;
	}
	else
	{
		MT_Data[max-1].src_UV_h_max=src_size_y;
		MT_Data[max-1].dst_UV_h_max=dst_size_y;
	}

	for (i=0; i<max; i++)
	{
		MT_Data[i].src_Y_w=src_size_x;
		MT_Data[i].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[i].src_UV_w=src_size_x >> UV_w;
			MT_Data[i].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[i].src_UV_w=src_size_x;
			MT_Data[i].dst_UV_w=dst_size_x;
		}
	}

	return(max);
}


void FilteredResizeH::ResamplerLumaMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_h_luma(MT_DataGF->dst1,MT_DataGF->src1,MT_DataGF->dst_pitch1,MT_DataGF->src_pitch1,
		MT_DataGF->resampling_program_luma,MT_DataGF->dst_Y_w,MT_DataGF->dst_Y_h_max-MT_DataGF->dst_Y_h_min,
		bits_per_pixel,plane_range[0],mode_YUY2);
}


void FilteredResizeH::ResamplerLumaMT2(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_h_luma(MT_DataGF->dst2,MT_DataGF->src2,MT_DataGF->dst_pitch2,MT_DataGF->src_pitch2,
		MT_DataGF->resampling_program_luma,MT_DataGF->dst_Y_w,MT_DataGF->dst_Y_h_max-MT_DataGF->dst_Y_h_min,
		bits_per_pixel,plane_range[1],mode_YUY2);
}


void FilteredResizeH::ResamplerLumaMT3(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_h_luma(MT_DataGF->dst3,MT_DataGF->src3,MT_DataGF->dst_pitch3,MT_DataGF->src_pitch3,
		MT_DataGF->resampling_program_luma,MT_DataGF->dst_Y_w,MT_DataGF->dst_Y_h_max-MT_DataGF->dst_Y_h_min,
		bits_per_pixel,plane_range[2],mode_YUY2);
}

void FilteredResizeH::ResamplerLumaMT4(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_h_luma(MT_DataGF->dst4,MT_DataGF->src4,MT_DataGF->dst_pitch4,MT_DataGF->src_pitch4,
		MT_DataGF->resampling_program_luma,MT_DataGF->dst_Y_w,MT_DataGF->dst_Y_h_max-MT_DataGF->dst_Y_h_min,
		bits_per_pixel,plane_range[3],mode_YUY2);
}

void FilteredResizeH::ResamplerUChromaMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_h_chroma(MT_DataGF->dst2,MT_DataGF->src2,MT_DataGF->dst_pitch2,MT_DataGF->src_pitch2,
		MT_DataGF->resampling_program_chroma,MT_DataGF->dst_UV_w,MT_DataGF->dst_UV_h_max-MT_DataGF->dst_UV_h_min,
		bits_per_pixel,plane_range[1],mode_YUY2);
}


void FilteredResizeH::ResamplerVChromaMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_h_chroma(MT_DataGF->dst3,MT_DataGF->src3,MT_DataGF->dst_pitch3,MT_DataGF->src_pitch3,
		MT_DataGF->resampling_program_chroma,MT_DataGF->dst_UV_w,MT_DataGF->dst_UV_h_max-MT_DataGF->dst_UV_h_min,
		bits_per_pixel,plane_range[2],mode_YUY2);
}


void FilteredResizeH::StaticThreadpoolH(void *ptr)
{
	Public_MT_Data_Thread *data=(Public_MT_Data_Thread *)ptr;
	FilteredResizeH *ptrClass=(FilteredResizeH *)data->pClass;
	MT_Data_Info_ResampleMT *MT_DataGF=((MT_Data_Info_ResampleMT *)data->pData)+data->thread_Id;

	switch(data->f_process)
	{
		case 1 : ptrClass->ResamplerLumaMT(MT_DataGF);
			break;
		case 2 : ptrClass->ResamplerUChromaMT(MT_DataGF);
			break;
		case 3 : ptrClass->ResamplerVChromaMT(MT_DataGF);
			break;
		case 4 : ptrClass->ResamplerLumaMT2(MT_DataGF);
			break;
		case 5 : ptrClass->ResamplerLumaMT3(MT_DataGF);
			break;
		case 6 : ptrClass->ResamplerLumaMT4(MT_DataGF);
			break;		
		default : ;
	}
}


PVideoFrame __stdcall FilteredResizeH::GetFrame(int n, IScriptEnvironment* env)
{
  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = (has_at_least_v8)?env->NewVideoFrameP(vi,&src):env->NewVideoFrame(vi,64);
  
  const int src_pitch_1 = src->GetPitch();
  const int dst_pitch_1 = dst->GetPitch();
  const BYTE *srcp_1 = src->GetReadPtr();
        BYTE *dstp_1 = dst->GetWritePtr();

	const int src_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_U) : (isRGBPfamily) ? src->GetPitch(PLANAR_B) : 0;
	const int dst_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_U) : (isRGBPfamily) ? dst->GetPitch(PLANAR_B) : 0;
	const BYTE *srcp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_U) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_B) : NULL;
	BYTE *dstp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_U) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_B) : NULL;

	const int src_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_V) : (isRGBPfamily) ? src->GetPitch(PLANAR_R) : 0;
	const int dst_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_V) : (isRGBPfamily) ? dst->GetPitch(PLANAR_R) : 0;
	const BYTE *srcp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_V) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_R) : NULL;
	BYTE *dstp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_V) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_R) : NULL;
	
	const int src_pitch_4 = (isAlphaChannel) ? src->GetPitch(PLANAR_A) : 0;
	const int dst_pitch_4 = (isAlphaChannel) ? dst->GetPitch(PLANAR_A) : 0;
	const BYTE *srcp_4 = (isAlphaChannel) ? src->GetReadPtr(PLANAR_A) : NULL;
	BYTE *dstp_4 = (isAlphaChannel) ? dst->GetWritePtr(PLANAR_A) : NULL;

	Public_MT_Data_Thread MT_ThreadGF[MAX_MT_THREADS];
	MT_Data_Info_ResampleMT MT_DataGF[MAX_MT_THREADS];

  memcpy(MT_ThreadGF,MT_Thread,sizeof(MT_Thread));
  memcpy(MT_DataGF,MT_Data,sizeof(MT_Data));

  for(uint8_t i=0; i<threads_number; i++)
	MT_ThreadGF[i].pData=(void *)MT_DataGF;
	
  int8_t idxPool=-1;

  if (threads_number>1)
  {
	if ((!poolInterface->RequestThreadPool(UserId,idxPool,threads_number,MT_ThreadGF)) || (idxPool==-1))
		env->ThrowError("ResizeHMT: Error with the TheadPool while requesting threadpool!");
  }
  
	for(uint8_t i=0; i<threads_number; i++)
	{
		MT_DataGF[i].src1=srcp_1+(MT_Data[i].src_Y_h_min*src_pitch_1);
		MT_DataGF[i].src2=srcp_2+(MT_Data[i].src_UV_h_min*src_pitch_2);
		MT_DataGF[i].src3=srcp_3+(MT_Data[i].src_UV_h_min*src_pitch_3);
		MT_DataGF[i].src4=srcp_4+(MT_Data[i].src_Y_h_min*src_pitch_4);
		MT_DataGF[i].src_pitch1=src_pitch_1;
		MT_DataGF[i].src_pitch2=src_pitch_2;
		MT_DataGF[i].src_pitch3=src_pitch_3;
		MT_DataGF[i].src_pitch4=src_pitch_4;
		MT_DataGF[i].dst1=dstp_1+(MT_Data[i].dst_Y_h_min*dst_pitch_1);
		MT_DataGF[i].dst2=dstp_2+(MT_Data[i].dst_UV_h_min*dst_pitch_2);
		MT_DataGF[i].dst3=dstp_3+(MT_Data[i].dst_UV_h_min*dst_pitch_3);
		MT_DataGF[i].dst4=dstp_4+(MT_Data[i].dst_Y_h_min*dst_pitch_4);
		MT_DataGF[i].dst_pitch1=dst_pitch_1;
		MT_DataGF[i].dst_pitch2=dst_pitch_2;
		MT_DataGF[i].dst_pitch3=dst_pitch_3;
		MT_DataGF[i].dst_pitch4=dst_pitch_4;
		MT_DataGF[i].filter_storage_luma=filter_storage_luma;
		MT_DataGF[i].resampling_program_luma=resampling_program_luma;
		MT_DataGF[i].resampling_program_chroma=resampling_program_chroma;
		MT_DataGF[i].filter_storage_chromaU=filter_storage_chroma;
		MT_DataGF[i].filter_storage_chromaV=filter_storage_chroma;
	}

	if (threads_number>1)
	{
		for(uint8_t i=0; i<threads_number; i++)
			MT_ThreadGF[i].f_process=1;
		if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);

		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=2;
			if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=3;
			if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);
		}
		else
		{
			if (isRGBPfamily)
			{
				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=4;
				if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);

				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=5;
				if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);								
			}
		}

		if (isAlphaChannel)
		{
			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=6;
			if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);												
		}

		for(uint8_t i=0; i<threads_number; i++)
			MT_ThreadGF[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep,idxPool);
	}
	else
	{
		// Do resizing
		ResamplerLumaMT(MT_DataGF);
    
		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			// Plane U resizing   
			ResamplerUChromaMT(MT_DataGF);
			// Plane V resizing
			ResamplerVChromaMT(MT_DataGF);
		}
		else
		{
			if (isRGBPfamily)
			{
				// Plane B resizing
				ResamplerLumaMT2(MT_DataGF);
				// Plane R resizing
				ResamplerLumaMT3(MT_DataGF);
			}
		}
		// Plane A resizing
		if (isAlphaChannel) ResamplerLumaMT4(MT_DataGF);
	}

  return dst;
}


ResamplerH FilteredResizeH::GetResampler(bool aligned, ResamplingProgram* program, IScriptEnvironment* env)
{
	int simd_coeff_count_padding = 8;

	if (Enable_SSSE3)
	{
		// both 8 and 16 bit SSSE3 and AVX2 horizontal resizer benefits from 16 pixels/cycle
		// float is also using 32 bytes, but as 32/sizeof(float) = 8, then don't need 16
		if (pixelsize == 1 || pixelsize == 2) simd_coeff_count_padding = 16;
   }

	resize_prepare_coeffs(program, env, simd_coeff_count_padding);

	if (pixelsize==1)
	{
		if (Enable_SSSE3)
		{
#ifdef AVX2_BUILD_POSSIBLE				
			if (Enable_AVX2) return resizer_h_avx2_generic_uint8_t;
			else
#endif			
			{
				if (program->filter_size_real>8) return resizer_h_ssse3_generic;
				else return resizer_h_ssse3_8;
			}
		}
		else return resize_h_c_planar;
	}
	else if (pixelsize==2)
	{ 
		if (Enable_SSSE3)
		{
#ifdef AVX2_BUILD_POSSIBLE				
			if (Enable_AVX2)
			{
				if(bits_per_pixel<16) return resizer_h_avx2_generic_uint16_t<true>;
				else return resizer_h_avx2_generic_uint16_t<false>;
			}
			else
#endif
			{
				if (Enable_SSE4_1) 
				{
					if (bits_per_pixel<16) return resizer_h_sse41_generic_uint16_t<true>;
					else return resizer_h_sse41_generic_uint16_t<false>;
				}
				else // SSSE3 needed
				{
					if (bits_per_pixel<16) return resizer_h_ssse3_generic_uint16_t<true>;
					else return resizer_h_ssse3_generic_uint16_t<false>;
				}
			}
		}
		else return resize_h_c_planar_s;
	}
	else
	{ //if (pixelsize == 4)
		if (Enable_SSSE3)
		{
			const int filtersizealign8 = AlignNumber(program->filter_size_real,8);
			const int filtersizemod8 = program->filter_size_real & 7;

#ifdef AVX2_BUILD_POSSIBLE
			if (Enable_AVX2)
			{
				if (filtersizealign8==8)
				{
					switch (filtersizemod8)
					{
						case 0 : return resizer_h_avx2_generic_float<1,0>; break;
						case 1 : return resizer_h_avx2_generic_float<1,1>; break;
						case 2 : return resizer_h_avx2_generic_float<1,2>; break;
						case 3 : return resizer_h_avx2_generic_float<1,3>; break;
						case 4 : return resizer_h_avx2_generic_float<1,4>; break;
						case 5 : return resizer_h_avx2_generic_float<1,5>; break;
						case 6 : return resizer_h_avx2_generic_float<1,6>; break;
						case 7 : return resizer_h_avx2_generic_float<1,7>; break;
						default : return NULL; break;
					}
				}
				else
				{
					if (filtersizealign8==16)
					{
						switch (filtersizemod8)
						{
							case 0 : return resizer_h_avx2_generic_float<2,0>; break;
							case 1 : return resizer_h_avx2_generic_float<2,1>; break;
							case 2 : return resizer_h_avx2_generic_float<2,2>; break;
							case 3 : return resizer_h_avx2_generic_float<2,3>; break;
							case 4 : return resizer_h_avx2_generic_float<2,4>; break;
							case 5 : return resizer_h_avx2_generic_float<2,5>; break;
							case 6 : return resizer_h_avx2_generic_float<2,6>; break;
							case 7 : return resizer_h_avx2_generic_float<2,7>; break;
							default : return NULL; break;
						}
					}
					else
					{
						switch (filtersizemod8)
						{
							case 0 : return resizer_h_avx2_generic_float<-1,0>; break;
							case 1 : return resizer_h_avx2_generic_float<-1,1>; break;
							case 2 : return resizer_h_avx2_generic_float<-1,2>; break;
							case 3 : return resizer_h_avx2_generic_float<-1,3>; break;
							case 4 : return resizer_h_avx2_generic_float<-1,4>; break;
							case 5 : return resizer_h_avx2_generic_float<-1,5>; break;
							case 6 : return resizer_h_avx2_generic_float<-1,6>; break;
							case 7 : return resizer_h_avx2_generic_float<-1,7>; break;
							default : return NULL; break;
						}
					}
				}
			}
			else
#endif		
			// SSSE3
			{
				if (filtersizealign8==8) 
				{
					switch (filtersizemod8)
					{
						case 0 : return resizer_h_ssse3_generic_float<1,0>; break;
						case 1 : return resizer_h_ssse3_generic_float<1,1>; break;
						case 2 : return resizer_h_ssse3_generic_float<1,2>; break;
						case 3 : return resizer_h_ssse3_generic_float<1,3>; break;
						case 4 : return resizer_h_ssse3_generic_float<1,4>; break;
						case 5 : return resizer_h_ssse3_generic_float<1,5>; break;
						case 6 : return resizer_h_ssse3_generic_float<1,6>; break;
						case 7 : return resizer_h_ssse3_generic_float<1,7>; break;
						default : return NULL; break;
					}
				}
				else
				{
					if (filtersizealign8==16)
					{
						switch (filtersizemod8)
						{
							case 0 : return resizer_h_ssse3_generic_float<2,0>; break;
							case 1 : return resizer_h_ssse3_generic_float<2,1>; break;
							case 2 : return resizer_h_ssse3_generic_float<2,2>; break;
							case 3 : return resizer_h_ssse3_generic_float<2,3>; break;
							case 4 : return resizer_h_ssse3_generic_float<2,4>; break;
							case 5 : return resizer_h_ssse3_generic_float<2,5>; break;
							case 6 : return resizer_h_ssse3_generic_float<2,6>; break;
							case 7 : return resizer_h_ssse3_generic_float<2,7>; break;
							default : return NULL; break;
						}
					}
					else
					{
						switch (filtersizemod8)
						{
							case 0 : return resizer_h_ssse3_generic_float<-1,0>; break;
							case 1 : return resizer_h_ssse3_generic_float<-1,1>; break;
							case 2 : return resizer_h_ssse3_generic_float<-1,2>; break;
							case 3 : return resizer_h_ssse3_generic_float<-1,3>; break;
							case 4 : return resizer_h_ssse3_generic_float<-1,4>; break;
							case 5 : return resizer_h_ssse3_generic_float<-1,5>; break;
							case 6 : return resizer_h_ssse3_generic_float<-1,6>; break;
							case 7 : return resizer_h_ssse3_generic_float<-1,7>; break;
							default : return NULL; break;
						}
					}
				}
			}
		}
		else return resize_h_c_planar_f;
	}
}


void FilteredResizeH::FreeData(void) 
{
	mydelete(resampling_program_luma);
	mydelete(resampling_program_chroma);
	
  myalignedfree(filter_storage_luma);
  myalignedfree(filter_storage_chroma);
}

FilteredResizeH::~FilteredResizeH(void)
{
	if (threads_number>1) poolInterface->RemoveUserId(UserId);
	FreeData();
	if (threads>1) poolInterface->DeAllocateAllThreads(true);
}


/***************************************
 ***** Filtered Resize - Vertical ******
 ***************************************/


FilteredResizeV::FilteredResizeV( PClip _child, double subrange_top, double subrange_height, int target_height,
	uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC, bool negativePrefetch,
	bool _avsp, bool preserve_center, ChromaLocation_e chroma_placement, bool ResizeH, ResamplingFunction* func, IScriptEnvironment* env )
  : GenericVideoFilter(_child),
    resampling_program_luma(NULL), resampling_program_chroma(NULL),
    src_pitch_table_luma(NULL), src_pitch_table_chromaU(NULL), src_pitch_table_chromaV(NULL),
    src_pitch_luma(-1), src_pitch_chromaU(-1), src_pitch_chromaV(-1),
    filter_storage_luma_aligned(NULL), filter_storage_luma_unaligned(NULL),
    filter_storage_chroma_aligned(NULL), filter_storage_chroma_unaligned(NULL),
	sleep(_sleep),threads(_threads),avsp(_avsp)
{
	int16_t i;

    pixelsize = (uint8_t)vi.ComponentSize(); // AVS16
	grey = vi.IsY();
	bits_per_pixel = (uint8_t)vi.BitsPerComponent();
	isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
	isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
	mode_YUY2 = vi.IsYUY2();

  Enable_MMX = (env->GetCPUFlags() & CPUF_MMX)!=0;
  Enable_SSE2 = (env->GetCPUFlags() & CPUF_SSE2)!=0;
  Enable_SSE3 = (env->GetCPUFlags() & CPUF_SSE3)!=0;
  Enable_SSSE3 = (env->GetCPUFlags() & CPUF_SSSE3)!=0;
  Enable_SSE4_1 = (env->GetCPUFlags() & CPUF_SSE4_1)!=0;
  Enable_AVX2 = false;
#ifdef AVX2_BUILD_POSSIBLE
  Enable_AVX2 = avsp && ((env->GetCPUFlags() & CPUF_AVX2)!=0);
#endif

  double center_pos_v_luma;
  double center_pos_v_chroma;
  GetCenterShiftForResizers(center_pos_v_luma, center_pos_v_chroma, preserve_center, chroma_placement, vi, ResizeH); // False for vertical

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;
	
    ResampleV_MT=StaticThreadpoolV;
	
	for (i=0; i<MAX_MT_THREADS; i++)
	{
		MT_Thread[i].pClass=this;
		MT_Thread[i].f_process=0;
		MT_Thread[i].thread_Id=(uint8_t)i;
		MT_Thread[i].pFunc=ResampleV_MT;
	}
	UserId=0;
	
	const int shift_w = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
	const int shift_h = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0;

	const int work_width = vi.IsPlanar() ? vi.width : vi.BytesFromPixels(vi.width)/pixelsize;
	
	if (vi.height<32) threads_number=1;
	else threads_number=threads;

  // Create resampling program and pitch table
  int SizeV;

  if (ShiftC==0) ShiftC=shift_h;
																																						
  if (desample) resampling_program_luma  = func->GetDesamplingProgram(target_height, subrange_top, subrange_height, vi.height, bits_per_pixel,
	  center_pos_v_luma, center_pos_v_luma, // for resizing it's the same for source and dest
	  accuracy, ChromaS, ShiftC, SizeV, env);
  else
  {
	  resampling_program_luma  = func->GetResamplingProgram(vi.height, subrange_top, subrange_height, target_height, bits_per_pixel,
		  center_pos_v_luma, center_pos_v_luma, // for resizing it's the same for source and dest
		  env);
	  SizeV=target_height;
  }

  if (resampling_program_luma==NULL)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  if (desample) env->ThrowError("ResizeVMT: Error while GetDesamplingProgram for luma!");
	  else env->ThrowError("ResizeVMT: Error while GetResamplingProgram for luma!");
  }

  if (desample && ((SizeV>vi.height) || (SizeV==-1)))
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  if (SizeV>vi.height) env->ThrowError("ResizeVMT: Desampling can only downscale!");
	  else env->ThrowError("ResizeVMT: Matrix can't be reversed!");
  }

  src_pitch_table_luma = (int *)_aligned_malloc(sizeof(int) * vi.height, 64);
  if (src_pitch_table_luma==NULL)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  env->ThrowError("ResizeVMT: Could not reserve memory in a resampler.");
  }
  
  if (vi.IsPlanar() && !grey && !isRGBPfamily)
  {
    const int div   = 1 << shift_h;

	if (desample)
	{
		int SizeOut;

	    resampling_program_chroma = func->GetDesamplingProgram(
		                              target_height  >> shift_h,
			                          subrange_top    / div,
				                      subrange_height / div,
					                  vi.height  >> shift_h,
									  bits_per_pixel,
									  center_pos_v_chroma, center_pos_v_chroma, // for resizing it's the same for source and dest
									  accuracy,SizeV,shift_h,SizeOut,
						              env);
		if (SizeOut==-1)
		{
			FreeData();
			if (threads>1) poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeVMT: Matrix can't be reversed!");
		}
	}
	else
	{
	    resampling_program_chroma = func->GetResamplingProgram(
		                              vi.height      >> shift_h,
			                          subrange_top    / div,
				                      subrange_height / div,
					                  target_height  >> shift_h,
									  bits_per_pixel,
									  center_pos_v_chroma, center_pos_v_chroma, // for resizing it's the same for source and dest
						              env);
	}
	if (resampling_program_chroma==NULL)
	{
		FreeData();
		if (threads>1) poolInterface->DeAllocateAllThreads(true);
		if (desample) env->ThrowError("ResizeVMT: Error while GetDesamplingProgram for chroma!");
		else env->ThrowError("ResizeVMT: Error while GetResamplingProgram for chroma!");
	}

	src_pitch_table_chromaU = (int *)_aligned_malloc(sizeof(int) * (vi.height >> shift_h), 64);
	src_pitch_table_chromaV = (int *)_aligned_malloc(sizeof(int) * (vi.height >> shift_h), 64);
	if ((src_pitch_table_chromaU==NULL) || (src_pitch_table_chromaV==NULL))
	{
		FreeData();
		if (threads>1) poolInterface->DeAllocateAllThreads(true);
		env->ThrowError("ResizeVMT: Could not reserve memory in a resampler.");
	}	

	const int h_UV=vi.height >> shift_h;

    resampler_chroma_aligned = GetResampler(true,resampling_program_chroma,env);
    resampler_chroma_unaligned = GetResampler(false,resampling_program_chroma,env);
  }

  resampler_luma_aligned   = GetResampler(true,resampling_program_luma,env);
  resampler_luma_unaligned = GetResampler(false,resampling_program_luma,env);

  threads_number=CreateMTData(threads_number,work_width,vi.height,work_width,SizeV,shift_w,shift_h);

	if (threads_number>1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeVMT: Error with the TheadPool while getting UserId!");
		}
		if (!poolInterface->EnableAllowSeveral(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeVMT: Error with the TheadPool while allowing multiple request on UserId!");
		}
		if (negativePrefetch)
		{
			if (!poolInterface->DisableWaitonRequest(UserId))
			{
				FreeData();
				poolInterface->DeAllocateAllThreads(true);
				env->ThrowError("ResizeVMT: Error with the TheadPool while disabling wait on request on UserId!");
			}
		}
	}

	has_at_least_v8=true;
	try { env->CheckVersion(8); } catch (const AvisynthError&) { has_at_least_v8=false; }

  // Change target video info size
  vi.height = SizeV;
}


int __stdcall FilteredResizeV::SetCacheHints(int cachehints,int frame_range)
{
  switch (cachehints)
  {
  case CACHE_GET_MTMODE :
    return MT_NICE_FILTER;
  default :
    return 0;
  }
}


uint8_t FilteredResizeV::CreateMTData(uint8_t max_threads,int32_t src_size_x,int32_t src_size_y,int32_t dst_size_x,int32_t dst_size_y,int UV_w,int UV_h)
{
	if ((max_threads<=1) || (max_threads>threads_number))
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	int32_t _y_min,_dh;
	int32_t src_dh_Y,src_dh_UV,dst_dh_Y,dst_dh_UV;
	int32_t h_y;
	uint8_t i,max_src=1,max_dst=1,max;

	dst_dh_Y=(dst_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
	if (dst_dh_Y<16) dst_dh_Y=16;
	if ((dst_dh_Y & 3)!=0) dst_dh_Y=((dst_dh_Y+3) >> 2) << 2;

	if (src_size_y==dst_size_y) src_dh_Y=dst_dh_Y;
	else
	{
		src_dh_Y=(src_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
		if (src_dh_Y<16) src_dh_Y=16;
		if ((src_dh_Y & 3)!=0) src_dh_Y=((src_dh_Y+3) >> 2) << 2;
	}

	_y_min=src_size_y;
	_dh=src_dh_Y;
	h_y=_dh;
	while (h_y<(_y_min-16))
	{
		max_src++;
		h_y+=_dh;
	}

	_y_min=dst_size_y;
	_dh=dst_dh_Y;
	h_y=_dh;
	while (h_y<(_y_min-16))
	{
		max_dst++;
		h_y+=_dh;
	}

	max=(max_src<max_dst) ? max_src:max_dst;

	if (max==1)
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	src_dh_UV= (UV_h>0) ? src_dh_Y>>UV_h : src_dh_Y;
	dst_dh_UV= (UV_h>0) ? dst_dh_Y>>UV_h : dst_dh_Y;

	MT_Data[0].top=true;
	MT_Data[0].bottom=false;
	MT_Data[0].src_Y_h_min=0;
	MT_Data[0].src_Y_h_max=src_dh_Y;
	MT_Data[0].dst_Y_h_min=0;
	MT_Data[0].dst_Y_h_max=dst_dh_Y;
	MT_Data[0].src_UV_h_min=0;
	MT_Data[0].src_UV_h_max=src_dh_UV;
	MT_Data[0].dst_UV_h_min=0;
	MT_Data[0].dst_UV_h_max=dst_dh_UV;

	i=1;
	while (i<max)
	{
		MT_Data[i].top=false;
		MT_Data[i].bottom=false;
		MT_Data[i].src_Y_h_min=MT_Data[i-1].src_Y_h_max;
		MT_Data[i].src_Y_h_max=MT_Data[i].src_Y_h_min+src_dh_Y;
		MT_Data[i].dst_Y_h_min=MT_Data[i-1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max=MT_Data[i].dst_Y_h_min+dst_dh_Y;
		MT_Data[i].src_UV_h_min=MT_Data[i-1].src_UV_h_max;
		MT_Data[i].src_UV_h_max=MT_Data[i].src_UV_h_min+src_dh_UV;
		MT_Data[i].dst_UV_h_min=MT_Data[i-1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max=MT_Data[i].dst_UV_h_min+dst_dh_UV;
		i++;
	}

	MT_Data[max-1].bottom=true;
	MT_Data[max-1].src_Y_h_max=src_size_y;
	MT_Data[max-1].dst_Y_h_max=dst_size_y;
	if (UV_h>0)
	{
		MT_Data[max-1].src_UV_h_max=src_size_y >> UV_h;
		MT_Data[max-1].dst_UV_h_max=dst_size_y >> UV_h;
	}
	else
	{
		MT_Data[max-1].src_UV_h_max=src_size_y;
		MT_Data[max-1].dst_UV_h_max=dst_size_y;
	}

	for (i=0; i<max; i++)
	{
		MT_Data[i].src_Y_w=src_size_x;
		MT_Data[i].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[i].src_UV_w=src_size_x >> UV_w;
			MT_Data[i].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[i].src_UV_w=src_size_x;
			MT_Data[i].dst_UV_w=dst_size_x;
		}
	}

	return(max);
}


void FilteredResizeV::ResamplerLumaAlignedMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_aligned(MT_DataGF->dst1,MT_DataGF->src1,MT_DataGF->dst_pitch1,MT_DataGF->src_pitch1,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma,plane_range[0],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_unaligned(MT_DataGF->dst1,MT_DataGF->src1,MT_DataGF->dst_pitch1,MT_DataGF->src_pitch1,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma,plane_range[0],mode_YUY2);
}

void FilteredResizeV::ResamplerLumaAlignedMT2(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_aligned(MT_DataGF->dst2,MT_DataGF->src2,MT_DataGF->dst_pitch2,MT_DataGF->src_pitch2,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma2,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT2(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_unaligned(MT_DataGF->dst2,MT_DataGF->src2,MT_DataGF->dst_pitch2,MT_DataGF->src_pitch2,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma2,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaAlignedMT3(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_aligned(MT_DataGF->dst3,MT_DataGF->src3,MT_DataGF->dst_pitch3,MT_DataGF->src_pitch3,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma3,plane_range[2],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT3(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_unaligned(MT_DataGF->dst3,MT_DataGF->src3,MT_DataGF->dst_pitch3,MT_DataGF->src_pitch3,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma3,plane_range[2],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaAlignedMT4(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_aligned(MT_DataGF->dst4,MT_DataGF->src4,MT_DataGF->dst_pitch4,MT_DataGF->src_pitch4,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma4,plane_range[3],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT4(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_luma_unaligned(MT_DataGF->dst4,MT_DataGF->src4,MT_DataGF->dst_pitch4,MT_DataGF->src_pitch4,
		MT_DataGF->resampling_program_luma,MT_DataGF->src_Y_w,bits_per_pixel,MT_DataGF->dst_Y_h_min,MT_DataGF->dst_Y_h_max,
		MT_DataGF->src_pitch_table_luma,MT_DataGF->filter_storage_luma4,plane_range[3],mode_YUY2);
}


void FilteredResizeV::ResamplerUChromaAlignedMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_chroma_aligned(MT_DataGF->dst2,MT_DataGF->src2,MT_DataGF->dst_pitch2,MT_DataGF->src_pitch2,
		MT_DataGF->resampling_program_chroma,MT_DataGF->src_UV_w,bits_per_pixel,MT_DataGF->dst_UV_h_min,MT_DataGF->dst_UV_h_max,
		MT_DataGF->src_pitch_table_chromaU,MT_DataGF->filter_storage_chromaU,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerUChromaUnalignedMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_chroma_unaligned(MT_DataGF->dst2,MT_DataGF->src2,MT_DataGF->dst_pitch2,MT_DataGF->src_pitch2,
		MT_DataGF->resampling_program_chroma,MT_DataGF->src_UV_w,bits_per_pixel,MT_DataGF->dst_UV_h_min,MT_DataGF->dst_UV_h_max,
		MT_DataGF->src_pitch_table_chromaU,MT_DataGF->filter_storage_chromaU,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerVChromaAlignedMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_chroma_aligned(MT_DataGF->dst3,MT_DataGF->src3,MT_DataGF->dst_pitch3,MT_DataGF->src_pitch3,
		MT_DataGF->resampling_program_chroma,MT_DataGF->src_UV_w,bits_per_pixel,MT_DataGF->dst_UV_h_min,MT_DataGF->dst_UV_h_max,
		MT_DataGF->src_pitch_table_chromaV,MT_DataGF->filter_storage_chromaV,plane_range[2],mode_YUY2);
}


void FilteredResizeV::ResamplerVChromaUnalignedMT(MT_Data_Info_ResampleMT *MT_DataGF)
{
	resampler_chroma_unaligned(MT_DataGF->dst3,MT_DataGF->src3,MT_DataGF->dst_pitch3,MT_DataGF->src_pitch3,
		MT_DataGF->resampling_program_chroma,MT_DataGF->src_UV_w,bits_per_pixel,MT_DataGF->dst_UV_h_min,MT_DataGF->dst_UV_h_max,
		MT_DataGF->src_pitch_table_chromaV,MT_DataGF->filter_storage_chromaV,plane_range[2],mode_YUY2);
}


void FilteredResizeV::StaticThreadpoolV(void *ptr)
{
	Public_MT_Data_Thread *data=(Public_MT_Data_Thread *)ptr;
	FilteredResizeV *ptrClass=(FilteredResizeV *)data->pClass;
	MT_Data_Info_ResampleMT *MT_DataGF=((MT_Data_Info_ResampleMT *)data->pData)+data->thread_Id;
	
	switch(data->f_process)
	{
		case 1 : ptrClass->ResamplerLumaAlignedMT(MT_DataGF);
			break;
		case 2 : ptrClass->ResamplerLumaUnalignedMT(MT_DataGF);
			break;
		case 3 : ptrClass->ResamplerUChromaAlignedMT(MT_DataGF);
			break;
		case 4 : ptrClass->ResamplerUChromaUnalignedMT(MT_DataGF);
			break;
		case 5 : ptrClass->ResamplerVChromaAlignedMT(MT_DataGF);
			break;
		case 6 : ptrClass->ResamplerVChromaUnalignedMT(MT_DataGF);
			break;
		case 7 : ptrClass->ResamplerLumaAlignedMT2(MT_DataGF);
			break;
		case 8 : ptrClass->ResamplerLumaUnalignedMT2(MT_DataGF);
			break;			
		case 9 : ptrClass->ResamplerLumaAlignedMT3(MT_DataGF);
			break;
		case 10 : ptrClass->ResamplerLumaUnalignedMT3(MT_DataGF);
			break;			
		case 11 : ptrClass->ResamplerLumaAlignedMT4(MT_DataGF);
			break;
		case 12 : ptrClass->ResamplerLumaUnalignedMT4(MT_DataGF);
			break;			
		default : ;
	}
}


PVideoFrame __stdcall FilteredResizeV::GetFrame(int n, IScriptEnvironment* env)
{
  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = (has_at_least_v8)?env->NewVideoFrameP(vi,&src):env->NewVideoFrame(vi,64);
  
  const int src_pitch_1 = src->GetPitch();
  const int dst_pitch_1 = dst->GetPitch();
  const BYTE *srcp_1 = src->GetReadPtr();
        BYTE *dstp_1 = dst->GetWritePtr();

	const int src_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_U) : (isRGBPfamily) ? src->GetPitch(PLANAR_B) : 0;
	const int dst_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_U) : (isRGBPfamily) ? dst->GetPitch(PLANAR_B) : 0;
	const BYTE *srcp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_U) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_B) : NULL;
	BYTE *dstp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_U) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_B) : NULL;

	const int src_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_V) : (isRGBPfamily) ? src->GetPitch(PLANAR_R) : 0;
	const int dst_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_V) : (isRGBPfamily) ? dst->GetPitch(PLANAR_R) : 0;
	const BYTE *srcp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_V) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_R) : NULL;
	BYTE *dstp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_V) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_R) : NULL;
	
	const int src_pitch_4 = (isAlphaChannel) ? src->GetPitch(PLANAR_A) : 0;
	const int dst_pitch_4 = (isAlphaChannel) ? dst->GetPitch(PLANAR_A) : 0;
	const BYTE *srcp_4 = (isAlphaChannel) ? src->GetReadPtr(PLANAR_A) : NULL;
	BYTE *dstp_4 = (isAlphaChannel) ? dst->GetWritePtr(PLANAR_A) : NULL;  

	Public_MT_Data_Thread MT_ThreadGF[MAX_MT_THREADS];
	MT_Data_Info_ResampleMT MT_DataGF[MAX_MT_THREADS];

  memcpy(MT_ThreadGF,MT_Thread,sizeof(MT_Thread));
  memcpy(MT_DataGF,MT_Data,sizeof(MT_Data));

  for(uint8_t i=0; i<threads_number; i++)
	MT_ThreadGF[i].pData=(void *)MT_DataGF;

  // Create pitch table
  if (src_pitch_luma != src->GetPitch())
  {
    resize_v_create_pitch_table(src_pitch_table_luma, src->GetPitch(), src->GetHeight(),pixelsize);
	src_pitch_luma = src->GetPitch();
  }

  if (!grey && vi.IsPlanar() && !isRGBPfamily)
  {
	if (src_pitch_chromaU != src->GetPitch(PLANAR_U))
	{
		resize_v_create_pitch_table(src_pitch_table_chromaU, src->GetPitch(PLANAR_U), src->GetHeight(PLANAR_U),pixelsize);
		src_pitch_chromaU = src->GetPitch(PLANAR_U);
	}	  
	if (src_pitch_chromaV != src->GetPitch(PLANAR_V))
	{
		resize_v_create_pitch_table(src_pitch_table_chromaV, src->GetPitch(PLANAR_V), src->GetHeight(PLANAR_V),pixelsize);
		src_pitch_chromaV = src->GetPitch(PLANAR_V);
	}	
  }

  int8_t idxPool=-1;

  if (threads_number>1)
  {
	if ((!poolInterface->RequestThreadPool(UserId,idxPool,threads_number,MT_ThreadGF)) || (idxPool==-1))
		env->ThrowError("ResizeVMT: Error with the TheadPool while requesting threadpool!");
  }

	for(uint8_t i=0; i<threads_number; i++)
	{		
		MT_DataGF[i].src1=srcp_1;
		MT_DataGF[i].src2=srcp_2;
		MT_DataGF[i].src3=srcp_3;
		MT_DataGF[i].src4=srcp_4;
		MT_DataGF[i].src_pitch1=src_pitch_1;
		MT_DataGF[i].src_pitch2=src_pitch_2;
		MT_DataGF[i].src_pitch3=src_pitch_3;
		MT_DataGF[i].src_pitch4=src_pitch_4;
		MT_DataGF[i].dst1=dstp_1+(MT_Data[i].dst_Y_h_min*dst_pitch_1);
		MT_DataGF[i].dst2=dstp_2+(MT_Data[i].dst_UV_h_min*dst_pitch_2);
		MT_DataGF[i].dst3=dstp_3+(MT_Data[i].dst_UV_h_min*dst_pitch_3);
		MT_DataGF[i].dst4=dstp_4+(MT_Data[i].dst_Y_h_min*dst_pitch_4);
		MT_DataGF[i].dst_pitch1=dst_pitch_1;
		MT_DataGF[i].dst_pitch2=dst_pitch_2;
		MT_DataGF[i].dst_pitch3=dst_pitch_3;
		MT_DataGF[i].dst_pitch4=dst_pitch_4;
		if (IsPtrAligned(srcp_1, 16) && ((src_pitch_1 & 15) == 0))
			MT_DataGF[i].filter_storage_luma=filter_storage_luma_aligned;
		else
			MT_DataGF[i].filter_storage_luma=filter_storage_luma_unaligned;
		if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0))
			MT_DataGF[i].filter_storage_luma4=filter_storage_luma_aligned;
		else
			MT_DataGF[i].filter_storage_luma4=filter_storage_luma_unaligned;
		MT_DataGF[i].src_pitch_table_luma=src_pitch_table_luma;
		MT_DataGF[i].src_pitch_table_chromaU=src_pitch_table_chromaU;
		MT_DataGF[i].src_pitch_table_chromaV=src_pitch_table_chromaV;
		MT_DataGF[i].resampling_program_luma=resampling_program_luma;
		MT_DataGF[i].resampling_program_chroma=resampling_program_chroma;
		if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0))
		{
			MT_DataGF[i].filter_storage_chromaU=filter_storage_chroma_aligned;
			MT_DataGF[i].filter_storage_luma2=filter_storage_luma_aligned;
		}
		else
		{
			MT_DataGF[i].filter_storage_chromaU=filter_storage_chroma_unaligned;
			MT_DataGF[i].filter_storage_luma2=filter_storage_luma_unaligned;
		}
		if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0))
		{
			MT_DataGF[i].filter_storage_chromaV=filter_storage_chroma_aligned;
			MT_DataGF[i].filter_storage_luma3=filter_storage_luma_aligned;
		}
		else
		{
			MT_DataGF[i].filter_storage_chromaV=filter_storage_chroma_unaligned;
			MT_DataGF[i].filter_storage_luma3=filter_storage_luma_unaligned;
		}
	}

	if (threads_number>1)
	{
		uint8_t f_proc;

		if (IsPtrAligned(srcp_1, 16) && ((src_pitch_1 & 15) == 0)) f_proc=1;
		else f_proc=2;

		for(uint8_t i=0; i<threads_number; i++)
			MT_ThreadGF[i].f_process=f_proc;
		if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);

		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0)) f_proc=3;
			else f_proc=4;

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);

			if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0)) f_proc=5;
			else f_proc=6;

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);
		}
		else
		{
			if (isRGBPfamily)
			{
				if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0)) f_proc=7;
				else f_proc=8;

				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=f_proc;
				if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);

				if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0)) f_proc=9;
				else f_proc=10;

				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=f_proc;
				if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);							
			}
		}
		
		if (isAlphaChannel)
		{
			if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0)) f_proc=11;
			else f_proc=12;

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId,idxPool)) poolInterface->WaitThreadsEnd(UserId,idxPool);			
		}

		for(uint8_t i=0; i<threads_number; i++)
			MT_ThreadGF[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep,idxPool);
	}
	else
	{
		// Do resizing
		if (IsPtrAligned(srcp_1, 16) && ((src_pitch_1 & 15) == 0))
			ResamplerLumaAlignedMT(MT_DataGF);
		else
			ResamplerLumaUnalignedMT(MT_DataGF);
    
		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			// Plane U resizing   
			if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0))
				ResamplerUChromaAlignedMT(MT_DataGF);
			else
				ResamplerUChromaUnalignedMT(MT_DataGF);

			// Plane V resizing
			if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0))
				ResamplerVChromaAlignedMT(MT_DataGF);
			else
				ResamplerVChromaUnalignedMT(MT_DataGF);
		}
		else
		{
			if (isRGBPfamily)
			{
				if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0))
					ResamplerLumaAlignedMT2(MT_DataGF);
				else
					ResamplerLumaUnalignedMT2(MT_DataGF);		
				
				if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0))
					ResamplerLumaAlignedMT3(MT_DataGF);
				else
					ResamplerLumaUnalignedMT3(MT_DataGF);								
			}			
		}
		
		if (isAlphaChannel)
		{
			if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0))
				ResamplerLumaAlignedMT4(MT_DataGF);
			else
				ResamplerLumaUnalignedMT4(MT_DataGF);	
		}
	}

  return dst;
}


ResamplerV FilteredResizeV::GetResampler(bool aligned, ResamplingProgram* program, IScriptEnvironment* env)
{
  resize_prepare_coeffs(program, env, 8); 

  if (program->filter_size_real==1)
  {
    // Fast pointresize
    switch (pixelsize) // AVS16
    {
    case 1: return resize_v_planar_pointresize<uint8_t>;
    case 2: return resize_v_planar_pointresize<uint16_t>;
    default: // case 4:
      return resize_v_planar_pointresize<float>;
    }
  }
  else
  {
    // Other resizers
    if (pixelsize==1)
    {
      if (Enable_SSSE3)
	  {
#ifdef AVX2_BUILD_POSSIBLE
		  if (aligned && Enable_AVX2) return resize_v_avx2_planar_uint8_t;
		  else
#endif
		  {
			if (aligned && Enable_SSE4_1)
			{
				return resize_v_sse41_planar;
			}
			else if (aligned)
			{ // SSSE3 aligned
				return resize_v_ssse3_planarT<simd_load_aligned>;
			}
			else if (Enable_SSE3)
			{ // SSE3 lddqu
				return resize_v_ssse3_planarT<simd_load_unaligned_sse3>;
			}
			else
			{ // SSSE3 unaligned
				return resize_v_ssse3_planarT<simd_load_unaligned>;
			}
		  }
      }
      else if (Enable_SSE2)
	  {
        if (aligned && Enable_SSE4_1)
		{ // SSE4.1 movntdqa constantly provide ~2% performance increase in my testing
          return resize_v_sse2_planar;
        }
        else if (aligned)
		{ // SSE2 aligned
          return resize_v_sse2_planarT<simd_load_aligned>;
        }
        else if (Enable_SSE3)
		{ // SSE2 lddqu
          return resize_v_ssse3_planar;
        }
        else
		{ // SSE2 unaligned
          return resize_v_sse2_planarT<simd_load_unaligned>;
        }
#ifdef X86_32
      }
      else if (Enable_MMX)
	  {
        return resize_v_mmx_planar;
#endif
      }
      else { // C version
        return resize_v_c_planar;
      }
    } 
    else if (pixelsize==2)
	{
#ifdef AVX2_BUILD_POSSIBLE		
		if (aligned && Enable_AVX2)
		{
			if(bits_per_pixel<16) return resize_v_avx2_planar_uint16_t<true>;
			else return resize_v_avx2_planar_uint16_t<false>;
		}
		else
#endif			
		if (aligned && Enable_SSE4_1)
		{
			if (bits_per_pixel<16) return resize_v_sse41_planar_uint16_t<true>;
			else return resize_v_sse41_planar_uint16_t<false>;
		}
		else if (aligned && Enable_SSE2)
		{
			if (bits_per_pixel<16) return resize_v_sse2_planar_uint16_t<true>;
			else return resize_v_sse2_planar_uint16_t<false>;
		}
		else
		{ // C version
			return resize_v_c_planar_s;
		}
    }
    else
	{ // if (pixelsize== 4) 
#ifdef AVX2_BUILD_POSSIBLE			
		if (aligned && Enable_AVX2) return resize_v_avx2_planar_float;
		else
#endif			
		if (aligned && Enable_SSE2) return resize_v_sse2_planar_float;
		else return resize_v_c_planar_f;
    }
  }
}


void FilteredResizeV::FreeData(void) 
{
	mydelete(resampling_program_luma);
	mydelete(resampling_program_chroma);
	myalignedfree(src_pitch_table_luma);
	myalignedfree(src_pitch_table_chromaU);
	myalignedfree(src_pitch_table_chromaV);

  myalignedfree(filter_storage_luma_aligned);
  myalignedfree(filter_storage_luma_unaligned);
  myalignedfree(filter_storage_chroma_aligned);
  myalignedfree(filter_storage_chroma_unaligned);
}


FilteredResizeV::~FilteredResizeV(void)
{
	if (threads_number>1) poolInterface->RemoveUserId(UserId);
	FreeData();
	if (threads>1) poolInterface->DeAllocateAllThreads(true);
}


/**********************************************
 *******   Resampling Factory Methods   *******
 **********************************************/

PClip FilteredResizeMT::CreateResizeH(PClip clip, double subrange_left, double subrange_width, int target_width,
                    uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,bool negativePrefetch,
					bool _avsp,bool preserve_center,ChromaLocation_e chroma_placement,
					ResamplingFunction* func,IScriptEnvironment* env)
{
	return new FilteredResizeH(clip, subrange_left, subrange_width, target_width,_threads,_sleep,range_mode,desample,
		accuracy,negativePrefetch,_avsp,preserve_center,chroma_placement,func, env);
}


PClip FilteredResizeMT::CreateResizeV(PClip clip, double subrange_top, double subrange_height, int target_height,
                    uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC,
					bool negativePrefetch,bool _avsp,bool preserve_center,ChromaLocation_e chroma_placement,
					bool ResizeH,ResamplingFunction* func,IScriptEnvironment* env)
{
	return new FilteredResizeV(clip, subrange_top, subrange_height, target_height,_threads,_sleep,range_mode,desample,
		accuracy,ChromaS,ShiftC,negativePrefetch,_avsp,preserve_center,chroma_placement,ResizeH,func,env);
}


PClip FilteredResizeMT::CreateResize(PClip clip, int target_width, int target_height, int force,int _threads,
	bool _LogicalCores,bool _MaxPhysCores, bool _SetAffinity,bool _sleep,int prefetch,
	int range_mode,bool desample,int accuracy,int order,int thread_level,
	const AVSValue* args,ResamplingFunction* f,
	bool preserve_center,const char *placement_name,ChromaLocation_e forced_chroma_placement,
	IScriptEnvironment* env)
{
  // 0 - return unchanged if no resize needed
  // 1 - force H
  // 2 - force V
  // 3 - force H and V
  const bool force_H = (force==1) || (force==3);
  const bool force_V = (force==2) || (force==3);

  const VideoInfo& vi = clip->GetVideoInfo();
  double subrange_left = args[0].AsFloat(0), subrange_top = args[1].AsFloat(0);
  const bool negativePrefetch=(prefetch<0)?true:false;
  
  if (target_height <= 0)
    env->ThrowError("ResizeMT: Height must be greater than 0.");

  if (target_width <= 0) {
    env->ThrowError("ResizeMT: Width must be greater than 0.");
  }

  if ((range_mode<0) || (range_mode>4)) env->ThrowError("ResizeMT: [range] must be between 0 and 4.");

  prefetch=abs(prefetch);
  if (prefetch>MAX_THREAD_POOL) env->ThrowError("ResizeMT: [prefetch] can't be higher than %d.",MAX_THREAD_POOL);
  if (prefetch==0) prefetch=1;
  if ((_threads<0) || (_threads>MAX_MT_THREADS)) env->ThrowError("ResizeMT: [threads] must be between 0 and %d.",MAX_MT_THREADS);
  if ((accuracy<0) || (accuracy>2)) env->ThrowError("ResizeMT: [accuracy] must be between 0 and 2.");
  if ((order<0) || (order>2)) env->ThrowError("ResizeMT: [order] must be between 0 and 2.");
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ResizeMT: [ThreadLevel] must be between 1 and 7.");

  const bool avsp=env->FunctionExists("ConvertBits");  
  const bool isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
  const bool grey = vi.IsY();  
  const bool isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
  const uint8_t bits_per_pixel = (uint8_t)vi.BitsPerComponent();


	bool has_at_least_v11=true;
	try { env->CheckVersion(11); } catch (const AvisynthError&) { has_at_least_v11=false; }

  // use forced_chroma_placement >= 0 and placement_name == nullptr together
  ChromaLocation_e chroma_placement = forced_chroma_placement >= 0 ? forced_chroma_placement : AVS_CHROMA_UNUSED;
  if (placement_name)
  {
    // no format-oriented defaults
    if (vi.IsYV411() || vi.Is420() || vi.Is422())
	{
      // placement explicite parameter like in ConvertToXXX or Text
      // input frame properties, if "auto"
      // When called from ConvertToXXX, chroma is not involved.
      auto frame0 = clip->GetFrame(0, env);
      const AVSMap* props = has_at_least_v11 ? env->getFramePropsRO(frame0):nullptr;
      chromaloc_parse_merge_with_props(vi, placement_name, props, /* ref*/chroma_placement, AVS_CHROMA_UNUSED /*default*/, env);
    }
  }


  if (vi.IsPlanar() && !grey && !isRGBPfamily)
  {
    int  mask;
	
	mask = (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1;
    if ((target_height & mask)!=0)
      env->ThrowError("ResizeMT: Planar destination height must be a multiple of %d.", mask+1);
  
    mask = (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1;
    if ((target_width & mask)!=0)
      env->ThrowError("ResizeMT: Planar destination width must be a multiple of %d.", mask+1);
  
  }

  double subrange_width=(desample)?args[2].AsDblDef(target_width):args[2].AsDblDef(vi.width);
  double subrange_height=(desample)?args[3].AsDblDef(target_height):args[3].AsDblDef(vi.height);

  // Crop style syntax
  if (desample)
  {
	if (subrange_width  <= 0.0) subrange_width  = target_width  - subrange_left + subrange_width;
	if (subrange_height <= 0.0) subrange_height = target_height - subrange_top  + subrange_height;
  }
  else
  {
	if (subrange_width  <= 0.0) subrange_width  = vi.width  - subrange_left + subrange_width;
	if (subrange_height <= 0.0) subrange_height = vi.height - subrange_top  + subrange_height;
  }

  // Packed RGB is bottom-top reversed
  if (vi.IsRGB() && !isRGBPfamily)
  {
	  if (desample)
		subrange_top = target_height - subrange_top - subrange_height;
	  else
		subrange_top = vi.height - subrange_top - subrange_height;
  }

  const int shift = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
  int SizeH;

  if (desample) SizeH=f->GetDesamplingData(target_width, subrange_left, subrange_width, vi.width,bits_per_pixel, 0.5, 0.5,shift,env);
  else SizeH=target_width;

  if (SizeH==-1) env->ThrowError("ResizeMT: Error while GetDesamplingData");

  bool fast_resize=((env->GetCPUFlags() & CPUF_SSSE3) == CPUF_SSSE3 ) && vi.IsPlanar();

  PClip result;
  // ensure that the intermediate area is maximal
  const double area_FirstH = (desample)?subrange_height*vi.width:subrange_height*target_width;
  const double area_FirstV = (desample)?subrange_width*vi.height:subrange_width*target_height;

  bool VFirst;

  if (desample)
  {
	  switch(order)
	  {
		case 0 : VFirst=(bits_per_pixel==32)?(area_FirstH<area_FirstV):(area_FirstH>=area_FirstV); break;
		case 1 : VFirst=true; break;
		case 2 : VFirst=false; break;
		default : VFirst=(bits_per_pixel==32)?(area_FirstH<area_FirstV):(area_FirstH>=area_FirstV); break;
	  }
  }
  else VFirst=(bits_per_pixel==32)?(area_FirstH>=area_FirstV):(area_FirstH<area_FirstV);

  const bool FTurnL=(!avsp) && (env->FunctionExists("FTurnLeft") && ((env->GetCPUFlags() & CPUF_SSE2)!=0)) && (!vi.IsRGB());
  const bool FTurnR=(!avsp) && (env->FunctionExists("FTurnRight") && ((env->GetCPUFlags() & CPUF_SSE2)!=0)) && (!vi.IsRGB());

  auto turnRightFunction = (FTurnR) ? "FTurnRight" : "TurnRight";
  auto turnLeftFunction =  (FTurnL) ? "FTurnLeft" : "TurnLeft";

  uint8_t plane_range[4];

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;

	bool step1,step2;

  bool CropV,CropH;

  if (desample)
  {
	  CropV=((subrange_top==int(subrange_top)) && (subrange_height==vi.height)
		&& (subrange_top>=0) && ((subrange_top+subrange_height)<= target_height));
	  CropH=((subrange_left==int(subrange_left)) && (subrange_width==vi.width)
		  && (subrange_left>=0) && ((subrange_left+subrange_width)<=target_width));
  }
  else
  {
	  CropV=((subrange_top==int(subrange_top)) && (subrange_height==target_height)
		&& (subrange_top>=0) && ((subrange_top+subrange_height)<= vi.height));
	  CropH=((subrange_left==int(subrange_left)) && (subrange_width==target_width)
		  && (subrange_left>=0) && ((subrange_left+subrange_width)<=vi.width));
  }


  if (VFirst)
  {
	if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)) step1=false;
	else
	{
		if (CropV)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;

			if (((int(subrange_top) | int(subrange_height)) & mask) == 0) step1=false;
			else step1=true;
		}
		else step1=true;
	}
	if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)))
	{
		if (CropH)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

		    if (((int(subrange_left) | int(subrange_width)) & mask) == 0) step2=false;
			else step2=true;
		}
		else step2=true;
	}
	else step2=false;
  }
  else
  {
	if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)) step1=false;
	else
	{
		if (CropH)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

		    if (((int(subrange_left) | int(subrange_width)) & mask) == 0) step1=false;
			else step1=true;
		}
		else step1=true;
	}
	if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)))
	{
		if (CropV)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
			if (((int(subrange_top) | int(subrange_height)) & mask) == 0) step2=false;
			else step2=true;
		}
		else step2=true;
	}
	else step2=false;
  }


	uint8_t threads_number=1;

	if ((_threads!=1) && (step1 || step2))
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ResizeMT: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(_threads,_LogicalCores);

		if (threads_number==0) env->ThrowError("ResizeMT: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (_SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,_MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ResizeMT: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,_MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ResizeMT: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,_MaxPhysCores,_SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ResizeMT: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

  if (!fast_resize)
  {
	  if (VFirst)
	  {
		if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height) && !force_V) result=clip;
		else
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;

				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,0,int(subrange_top),vi.width,int(subrange_height),true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
					preserve_center,chroma_placement,false,f,env);
			}
			else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
				preserve_center,chroma_placement,false,f,env);
		}
		if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width) && !force_H))
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,int(subrange_left),0,int(subrange_width),target_height,true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else
				{
					if (!vi.IsRGB() || isRGBPfamily)
					{
						if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
						{
							const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
							const int div   = 1 << shift;

							AVSValue v,vv,vu,va;
							
							if (avsp)
							{
								AVSValue sargs[2] = {result,"Y"};
								
								v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="U";
								vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="V";
								vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								if (isAlphaChannel)
								{
									sargs[1]="A";
									va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								}
							}
							else
							{
								vu = env->Invoke("UtoY8",result).AsClip();
								vv = env->Invoke("VtoY8",result).AsClip();
								v = env->Invoke("ConvertToY8",result).AsClip();								
							}
							
							v = env->Invoke(turnRightFunction,v).AsClip();
							vu = env->Invoke(turnRightFunction,vu).AsClip();
							vv = env->Invoke(turnRightFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
							v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[0],desample,accuracy,0,shift,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);

							VideoInfo vR = v.AsClip()->GetVideoInfo();
							int ChromaS=vR.height >> shift;

							vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[1],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[2],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[3],desample,accuracy,0,shift,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							v = env->Invoke(turnLeftFunction,v).AsClip();
							vu = env->Invoke(turnLeftFunction,vu).AsClip();
							vv = env->Invoke(turnLeftFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();

						    AVSValue ytouvargs[4] = {vu,vv,v,va};
						    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
							else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
						    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
						}
						else
						{
							result=env->Invoke(turnRightFunction,result).AsClip();
							result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							result=env->Invoke(turnLeftFunction,result).AsClip();
						}
					}
					else
					{
						result=env->Invoke(turnLeftFunction,result).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						result=env->Invoke(turnRightFunction,result).AsClip();
					}
				}
			}
			else
			{
				if (!vi.IsRGB() || isRGBPfamily)
				{
				    if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
					{
						const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
						const int div   = 1 << shift;

						AVSValue v,vv,vu,va;
						
						if (avsp)
						{
							AVSValue sargs[2] = {result,"Y"};
								
							v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="U";
							vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="V";
							vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							if (isAlphaChannel)
							{
								sargs[1]="A";
								va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							}
						}
						else
						{
							vu = env->Invoke("UtoY8",result).AsClip();
							vv = env->Invoke("VtoY8",result).AsClip();
							v = env->Invoke("ConvertToY8",result).AsClip();								
						}
			
						v = env->Invoke(turnRightFunction,v).AsClip();
						vu = env->Invoke(turnRightFunction,vu).AsClip();
						vv = env->Invoke(turnRightFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
						v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[0],desample,accuracy,0,shift,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);

						VideoInfo vR = v.AsClip()->GetVideoInfo();
						int ChromaS=vR.height >> shift;

						vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[1],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[2],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[3],desample,accuracy,0,shift,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						v = env->Invoke(turnLeftFunction,v).AsClip();
						vu = env->Invoke(turnLeftFunction,vu).AsClip();
						vv = env->Invoke(turnLeftFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();
							
					    AVSValue ytouvargs[4] = {vu,vv,v,va};
					    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
						else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
					    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
					}
					else
					{
						result=env->Invoke(turnRightFunction,result).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						result=env->Invoke(turnLeftFunction,result).AsClip();
					}
				}
				else
				{
					result=env->Invoke(turnLeftFunction,result).AsClip();
					result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
						preserve_center,chroma_placement,true,f,env);
					result=env->Invoke(turnRightFunction,result).AsClip();
				}
			}
		}
	  }
	  else
	  {
		if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width) && !force_H) result=clip;
		else
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,int(subrange_left),0,int(subrange_width),vi.height,true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else
				{
					if (!vi.IsRGB() || isRGBPfamily)
					{
					    if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
						{
							const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
							const int div   = 1 << shift;

							AVSValue v,vv,vu,va;
							
							if (avsp)
							{
								AVSValue sargs[2] = {clip,"Y"};
								
								v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="U";
								vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="V";
								vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								if (isAlphaChannel)
								{
									sargs[1]="A";
									va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								}
							}
							else
							{
								vu = env->Invoke("UtoY8",clip).AsClip();
								vv = env->Invoke("VtoY8",clip).AsClip();
								v = env->Invoke("ConvertToY8",clip).AsClip();								
							}		

							v = env->Invoke(turnRightFunction,v).AsClip();
							vu = env->Invoke(turnRightFunction,vu).AsClip();
							vv = env->Invoke(turnRightFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
							v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[0],desample,accuracy,0,shift,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);

							VideoInfo vR = v.AsClip()->GetVideoInfo();
							int ChromaS=vR.height >> shift;

							vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[1],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[2],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[3],desample,accuracy,0,shift,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							v = env->Invoke(turnLeftFunction,v).AsClip();
							vu = env->Invoke(turnLeftFunction,vu).AsClip();
							vv = env->Invoke(turnLeftFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();

						    AVSValue ytouvargs[4] = {vu,vv,v,va};
						    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
							else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
						    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
						}
						else
						{
							result=env->Invoke(turnRightFunction,clip).AsClip();
							result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
								preserve_center,chroma_placement,true,f,env);
							result=env->Invoke(turnLeftFunction,result).AsClip();
						}
					}
					else
					{
						result=env->Invoke(turnLeftFunction,clip).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						result=env->Invoke(turnRightFunction,result).AsClip();
					}
				}
			}
			else
			{
				if (!vi.IsRGB() || isRGBPfamily)
				{
					if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
					{
						const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
						const int div   = 1 << shift;

						AVSValue v,vv,vu,va;
						
						if (avsp)
						{
							AVSValue sargs[2] = {clip,"Y"};
								
							v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="U";
							vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="V";
							vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							if (isAlphaChannel)
							{
								sargs[1]="A";
								va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							}
						}
						else
						{
							vu = env->Invoke("UtoY8",clip).AsClip();
							vv = env->Invoke("VtoY8",clip).AsClip();
							v = env->Invoke("ConvertToY8",clip).AsClip();								
						}
						
						v = env->Invoke(turnRightFunction,v).AsClip();
						vu = env->Invoke(turnRightFunction,vu).AsClip();
						vv = env->Invoke(turnRightFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
						v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[0],desample,accuracy,0,shift,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);

						VideoInfo vR = v.AsClip()->GetVideoInfo();
						int ChromaS=vR.height >> shift;

						vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[1],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[2],desample,accuracy,ChromaS,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[3],desample,accuracy,0,shift,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						v = env->Invoke(turnLeftFunction,v).AsClip();
						vu = env->Invoke(turnLeftFunction,vu).AsClip();
						vv = env->Invoke(turnLeftFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();

					    AVSValue ytouvargs[4] = {vu,vv,v,va};
					    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
						else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
					    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
					}
					else
					{
						result=env->Invoke(turnRightFunction,clip).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
							preserve_center,chroma_placement,true,f,env);
						result=env->Invoke(turnLeftFunction,result).AsClip();
					}
				}
				else
				{
					result=env->Invoke(turnLeftFunction,clip).AsClip();
					result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
						preserve_center,chroma_placement,true,f,env);
					result=env->Invoke(turnRightFunction,result).AsClip();
				}
			}
		}
		if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height) && !force_V))
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,0,int(subrange_top),target_width,int(subrange_height),true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
					preserve_center,chroma_placement,false,f,env);
			}
			else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
				preserve_center,chroma_placement,false,f,env);
		}
	  }
  }
  else
  {	  
	  if (VFirst)
	  {
		if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height) && !force_V) result=clip;
		else
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;

				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,0,int(subrange_top),vi.width,int(subrange_height),true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
					preserve_center,chroma_placement,false,f,env);
			}
			else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
				preserve_center,chroma_placement,false,f,env);
		}
		if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width) && !force_H))
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,int(subrange_left),0,int(subrange_width),target_height,true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeH(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,negativePrefetch,avsp,
					preserve_center,chroma_placement,f,env);
			}
			else result = CreateResizeH(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,negativePrefetch,avsp,
				preserve_center,chroma_placement,f,env);
		}		
	  }
	  else
	  {
		if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width) && !force_H) result=clip;
		else
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,int(subrange_left),0,int(subrange_width),vi.height,true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeH(clip, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,negativePrefetch,avsp,
					preserve_center,chroma_placement,f,env);
			}
			else result = CreateResizeH(clip, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,negativePrefetch,avsp,
				preserve_center,chroma_placement,f,env);
		}
		if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height) && !force_V))
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,0,int(subrange_top),target_width,int(subrange_height),true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
					preserve_center,chroma_placement,false,f,env);
			}
			else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,negativePrefetch,avsp,
				preserve_center,chroma_placement,false,f,env);
		}
	  }
  }
  
  return result;
}


AVSValue __cdecl FilteredResizeMT::Create_PointResize(AVSValue args, void*, IScriptEnvironment* env)
{
  PointFilter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_BilinearResize(AVSValue args, void*, IScriptEnvironment* env)
{
  TriangleFilter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_BicubicResize(AVSValue args, void*, IScriptEnvironment* env)
{
  MitchellNetravaliFilter f(args[3].AsDblDef(1./3.), args[4].AsDblDef(1./3.));
  const int force = args[9].AsInt(0);

  bool preserve_center = args[10].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[11].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 12;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_LanczosResize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(args[7].AsInt(3));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Lanczos4Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(4);
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_BlackmanResize(AVSValue args, void*, IScriptEnvironment* env)
{
  BlackmanFilter f(args[7].AsInt(4));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Spline16Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline16Filter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Spline36Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline36Filter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Spline64Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline64Filter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_GaussianResize(AVSValue args, void*, IScriptEnvironment* env)
{
  GaussianFilter f(args[7].AsFloat(30.0f),args[8].AsFloat(2.0f),args[9].AsFloat(4.0f));
  const int force = args[10].AsInt(0);

  bool preserve_center = args[11].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[12].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 13;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_SincResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincFilter f(args[7].AsInt(4));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_SinPowerResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SinPowerFilter f(args[7].AsFloat(2.5f));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_SincLin2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincLin2Filter f(args[7].AsInt(15));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_UserDefined2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  UserDefined2Filter f(args[3].AsFloat(121.0f),args[4].AsFloat(19.0f), args[5].AsFloat(2.3f));
  const int force = args[10].AsInt(0);

  bool preserve_center = args[11].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[12].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 13;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),false,0,0,args[Offset_Arg+7].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

// Desample functions

AVSValue __cdecl FilteredResizeMT::Create_DeBilinearResize(AVSValue args, void*, IScriptEnvironment* env)
{
  TriangleFilter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeBicubicResize(AVSValue args, void*, IScriptEnvironment* env)
{
  MitchellNetravaliFilter f(args[3].AsDblDef(1./3.), args[4].AsDblDef(1./3.));
  const int force = args[9].AsInt(0);

  bool preserve_center = args[10].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[11].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 12;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeLanczosResize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(args[7].AsInt(3));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeLanczos4Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(4);
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeBlackmanResize(AVSValue args, void*, IScriptEnvironment* env)
{
  BlackmanFilter f(args[7].AsInt(4));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline16Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline16Filter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline36Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline36Filter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline64Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline64Filter f;
  const int force = args[7].AsInt(0);

  bool preserve_center = args[8].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[9].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 10;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeGaussianResize(AVSValue args, void*, IScriptEnvironment* env)
{
  GaussianFilter f(args[7].AsFloat(30.0f),args[8].AsFloat(2.0f),args[9].AsFloat(4.0f));
  const int force = args[10].AsInt(0);

  bool preserve_center = args[11].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[12].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 13;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSincResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincFilter f(args[7].AsInt(4));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSinPowerResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SinPowerFilter f(args[7].AsFloat(2.5f));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSincLin2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincLin2Filter f(args[7].AsInt(15));
  const int force = args[8].AsInt(0);

  bool preserve_center = args[9].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[10].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 11;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeUserDefined2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  UserDefined2Filter f(args[3].AsFloat(121.0f),args[4].AsFloat(19.0f),args[5].AsFloat(2.3f));
  const int force = args[10].AsInt(0);

  bool preserve_center = args[11].AsBool(true); // [keep_center] default Avisynth
  const char* placement_name = args[12].AsString("auto"); // [placement]s
  const ChromaLocation_e forced_chroma_placement = AVS_CHROMA_UNUSED; // no force, used internally

  const unsigned short Offset_Arg = 13;

  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),force,args[Offset_Arg].AsInt(0),
	  args[Offset_Arg+1].AsBool(true),args[Offset_Arg+2].AsBool(true),args[Offset_Arg+3].AsBool(false),args[Offset_Arg+4].AsBool(false),
	  args[Offset_Arg+5].AsInt(0),args[Offset_Arg+6].AsInt(1),true,args[Offset_Arg+7].AsInt(0),args[Offset_Arg+8].AsInt(0),args[Offset_Arg+9].AsInt(6),&args[3],&f,
	  preserve_center,placement_name,forced_chroma_placement,env);
}
