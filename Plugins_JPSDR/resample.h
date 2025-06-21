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

#ifndef __Resample_H__
#define __Resample_H__

//#include <stdint.h>
#include <windows.h>
#include "./avisynth.h"
#include "./resample_functions.h"
#include "./ThreadPoolInterface.h"

#define RESAMPLE_MT_VERSION "ResampleMT 2.6.0 JPSDR"

typedef enum ChromaLocation_e
{
  AVS_CHROMA_UNUSED = -1,
  AVS_CHROMA_LEFT = 0,
  AVS_CHROMA_CENTER = 1,
  AVS_CHROMA_TOP_LEFT = 2,
  AVS_CHROMA_TOP = 3,
  AVS_CHROMA_BOTTOM_LEFT = 4,
  AVS_CHROMA_BOTTOM = 5,
  AVS_CHROMA_DV = 6 // Special to Avisynth
} ChromaLocation_e;

// Resizer function pointer
typedef void (*ResamplerV)(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2);
typedef void (*ResamplerH)(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int target_height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);


typedef struct _MT_Data_Info_ResampleMT
{
	const BYTE*src1,*src2,*src3,*src4;
	BYTE *dst1,*dst2,*dst3,*dst4;
	int src_pitch1,src_pitch2,src_pitch3,src_pitch4;
	int dst_pitch1,dst_pitch2,dst_pitch3,dst_pitch4;
	int32_t src_Y_h_min,src_Y_h_max,src_Y_w;
	int32_t src_UV_h_min,src_UV_h_max,src_UV_w;
	int32_t dst_Y_h_min,dst_Y_h_max,dst_Y_w;
	int32_t dst_UV_h_min,dst_UV_h_max,dst_UV_w;
	void *filter_storage_luma,*filter_storage_luma2,*filter_storage_luma3,*filter_storage_luma4;
	void *filter_storage_chromaU,*filter_storage_chromaV;
	int *src_pitch_table_luma,*src_pitch_table_chromaU,*src_pitch_table_chromaV;
	ResamplingProgram *resampling_program_luma,*resampling_program_chroma;
	bool top,bottom;
} MT_Data_Info_ResampleMT;



/**
  * Class to resize in the horizontal direction using a specified sampling filter
  * Helper for resample functions
 **/
class FilteredResizeH : public GenericVideoFilter
{
public:
  FilteredResizeH( PClip _child, double subrange_left, double subrange_width, int target_width, uint8_t _threads,
	  bool _sleep,int range_mode,bool desample,int accuracy, bool negativePrefetch,
	  bool _avsp,bool preserve_center,ChromaLocation_e chroma_placement,
	  ResamplingFunction* func,IScriptEnvironment* env );
  virtual ~FilteredResizeH(void);
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);

  int __stdcall SetCacheHints(int cachehints, int frame_range);

  //static ResamplerH GetResampler(int CPU, bool aligned, int pixelsize, int bits_per_pixel, ResamplingProgram* program, IScriptEnvironment* env);
  ResamplerH GetResampler(bool aligned, ResamplingProgram* program, IScriptEnvironment* env);

private:
	Public_MT_Data_Thread MT_Thread[MAX_MT_THREADS];
	MT_Data_Info_ResampleMT MT_Data[MAX_MT_THREADS];
	uint8_t threads,threads_number;
	bool sleep;
	uint32_t UserId;
	
	ThreadPoolFunction ResampleH_MT;

	static void StaticThreadpoolH(void *ptr);

	uint8_t CreateMTData(uint8_t max_threads,int32_t src_size_x,int32_t src_size_y,int32_t dst_size_x,int32_t dst_size_y, int UV_w, int UV_h);

	void FreeData(void);

	void ResamplerLumaMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaMT2(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaMT3(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaMT4(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerUChromaMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerVChromaMT(MT_Data_Info_ResampleMT *MT_DataGF);
	

  // Resampling
  ResamplingProgram *resampling_program_luma;
  ResamplingProgram *resampling_program_chroma;

  // Note: these pointer are currently not used; they are used to pass data into run-time resampler.
  // They are kept because this may be needed later (like when we implemented actual horizontal resizer.)
  void* filter_storage_luma;
  void* filter_storage_chroma;

  int src_width, src_height, dst_width, dst_height;
  bool grey,avsp,isRGBPfamily,isAlphaChannel,has_at_least_v8;
  uint8_t pixelsize; // AVS16
  uint8_t bits_per_pixel;
  uint8_t plane_range[4];
  bool mode_YUY2;
  bool Enable_MMX,Enable_SSE2,Enable_SSE3,Enable_SSSE3,Enable_SSE4_1,Enable_AVX2;

  ResamplerH resampler_h_luma;
  ResamplerH resampler_h_chroma;
};


/**
  * Class to resize in the vertical direction using a specified sampling filter
  * Helper for resample functions
 **/
class FilteredResizeV : public GenericVideoFilter
{
public:
  FilteredResizeV( PClip _child, double subrange_top, double subrange_height, int target_height, uint8_t _threads,
	  bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC,bool negativePrefetch,
	  bool _avsp,bool preserve_center,ChromaLocation_e chroma_placement,
	  bool ResizeH,ResamplingFunction* func,IScriptEnvironment* env);
  virtual ~FilteredResizeV(void);
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);

	int __stdcall SetCacheHints(int cachehints, int frame_range);

  //static ResamplerV GetResampler(int CPU, bool aligned,int pixelsize, int bits_per_pixel, void*& storage, ResamplingProgram* program);
  //ResamplerV GetResampler(bool aligned,void*& storage, ResamplingProgram* program);
  ResamplerV GetResampler(bool aligned, ResamplingProgram* program, IScriptEnvironment* env);

private:
	Public_MT_Data_Thread MT_Thread[MAX_MT_THREADS];
	MT_Data_Info_ResampleMT MT_Data[MAX_MT_THREADS];
	uint8_t threads,threads_number;
	bool sleep;
	uint32_t UserId;

	ThreadPoolFunction ResampleV_MT;

	static void StaticThreadpoolV(void *ptr);

	uint8_t CreateMTData(uint8_t max_threads,int32_t src_size_x,int32_t src_size_y,int32_t dst_size_x,int32_t dst_size_y, int UV_w, int UV_h);

	void FreeData(void);

	void ResamplerLumaAlignedMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaUnalignedMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaAlignedMT2(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaUnalignedMT2(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaAlignedMT3(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaUnalignedMT3(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaAlignedMT4(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerLumaUnalignedMT4(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerUChromaAlignedMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerUChromaUnalignedMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerVChromaAlignedMT(MT_Data_Info_ResampleMT *MT_DataGF);
	void ResamplerVChromaUnalignedMT(MT_Data_Info_ResampleMT *MT_DataGF);

  bool grey,avsp,isRGBPfamily,isAlphaChannel,has_at_least_v8;
  uint8_t pixelsize; // AVS16
  uint8_t bits_per_pixel;
  uint8_t plane_range[4];
  bool mode_YUY2;
  bool Enable_MMX,Enable_SSE2,Enable_SSE3,Enable_SSSE3,Enable_SSE4_1,Enable_AVX2;
	
  ResamplingProgram *resampling_program_luma;
  ResamplingProgram *resampling_program_chroma;
  int *src_pitch_table_luma;
  int *src_pitch_table_chromaU;
  int *src_pitch_table_chromaV;
  int src_pitch_luma;
  int src_pitch_chromaU;
  int src_pitch_chromaV;

  // Note: these pointer are currently not used; they are used to pass data into run-time resampler.
  // They are kept because this may be needed later (like when we implemented actual horizontal resizer.)
  void* filter_storage_luma_aligned;
  void* filter_storage_luma_unaligned;
  void* filter_storage_chroma_aligned;
  void* filter_storage_chroma_unaligned;

  ResamplerV resampler_luma_aligned;
  ResamplerV resampler_luma_unaligned;
  ResamplerV resampler_chroma_aligned;
  ResamplerV resampler_chroma_unaligned;
};



class FilteredResizeMT
{
public:
static PClip CreateResizeH( PClip clip, double subrange_left, double subrange_width, int target_width, uint8_t _threads,
	                         bool _sleep,int range_mode,bool desample,int accuracy, bool negativePrefetch,
							 bool _avsp, bool preserve_center,ChromaLocation_e chroma_placement,
							 ResamplingFunction* func,IScriptEnvironment* env );
static PClip CreateResizeV( PClip clip, double subrange_top, double subrange_height, int target_height, uint8_t _threads,
	                         bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC, bool negativePrefetch,
							 bool _avsp,bool preserve_center,ChromaLocation_e chroma_placement,
							 bool ResizeH,ResamplingFunction* func,IScriptEnvironment* env );

static PClip CreateResize( PClip clip, int target_width, int target_height, int force, int _threads,
	bool _LogicalCores,bool _MaxPhysCores, bool _SetAffinity,bool _sleep,int prefetch,int range_mode,
	bool desample,int accuracy,int order,int thread_level,
	const AVSValue* args,ResamplingFunction* f,
	bool preserve_center,const char *placement_name,ChromaLocation_e forced_chroma_placement,
	IScriptEnvironment* env );

static AVSValue __cdecl Create_PointResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_BilinearResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_BicubicResize(AVSValue args, void*, IScriptEnvironment* env);

// 09-14-2002 - Vlad59 - Lanczos3Resize - 
static AVSValue __cdecl Create_LanczosResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_Lanczos4Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_BlackmanResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_Spline16Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_Spline36Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_Spline64Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_GaussianResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_SincResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_SinPowerResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_SincLin2Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_UserDefined2Resize(AVSValue args, void*, IScriptEnvironment* env);

// Desample functions

static AVSValue __cdecl Create_DeBilinearResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeBicubicResize(AVSValue args, void*, IScriptEnvironment* env);

// 09-14-2002 - Vlad59 - Lanczos3Resize - 
static AVSValue __cdecl Create_DeLanczosResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeLanczos4Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeBlackmanResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeSpline16Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeSpline36Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeSpline64Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeGaussianResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeSincResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeSinPowerResize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeSincLin2Resize(AVSValue args, void*, IScriptEnvironment* env);

static AVSValue __cdecl Create_DeUserDefined2Resize(AVSValue args, void*, IScriptEnvironment* env);
};


#endif // __Resample_H__


