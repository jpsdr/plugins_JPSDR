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

#ifndef __JINCRESIZEMT_H__
#define __JINCRESIZEMT_H__

#include <stdint.h>
#include <vector>
#include <immintrin.h>

#include "avisynth.h"
#include "ThreadPoolInterface.h"

#define JINCRESIZEMT_VERSION "JincResizeMT 1.1.0 JPSDR"

#define JincMT_RESTRICT __restrict

// VS 2017 v15.3
#if _MSC_VER >= 1911
#define JincMT_CONSTEXPR constexpr
#else
#define JincMT_CONSTEXPR
#endif

#ifndef __CHROMALOCATION__
#define __CHROMALOCATION__
typedef enum _ChromaLocation_e
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
#endif

typedef struct _MT_Data_Info_JincResizeMT
{
	const BYTE *src[4];
	BYTE *JincMT_RESTRICT dst[4];
	ptrdiff_t src_pitch[4];
	ptrdiff_t dst_pitch[4];
	int32_t src_Y_h_min, src_Y_h_max, src_Y_w;
	int32_t src_UV_h_min, src_UV_h_max, src_UV_w;
	int32_t dst_Y_h_min, dst_Y_h_max, dst_Y_w;
	int32_t dst_UV_h_min, dst_UV_h_max, dst_UV_w;
	bool top, bottom;
} MT_Data_Info_JincResizeMT;

struct EWAPixelCoeffMeta
{
    int start_x;
    int start_y;
    int coeff_meta;
};

struct EWAPixelCoeff
{
    float *factor;
    EWAPixelCoeffMeta *meta;
    int *factor_map;
    int filter_size;
    int coeff_stride;
	
	EWAPixelCoeff() : factor(nullptr), meta(nullptr), factor_map(nullptr) {}
};

#define LUT_SIZE_VALUE 1024

class Lut
{
    int lut_size;

public:
    Lut();
	virtual ~Lut();
	bool InitLut(int lutsize, double radius, double blur);
    float GetFactor(int index);

    double* lut;
};

typedef void (*JincResizeMT_Process)(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *coeff,
	const float Val_Min[], const float Val_Max[]);

class JincResizeMT : public GenericVideoFilter
{
    Lut *init_lut;
	std::vector<EWAPixelCoeff*> out;
    bool avx512,avx2,sse41;
    uint8_t planecount;
    bool has_at_least_v8,has_at_least_v11;
	bool grey,isRGBPfamily,isAlphaChannel;
	uint8_t bits_per_pixel;
	bool subsampled;
    float ValMin[4],ValMax[4];

	ChromaLocation_e chroma_placement;

	JincResizeMT_Process process_frame_1x, process_frame_2x, process_frame_3x, process_frame_4x;

	Public_MT_Data_Thread MT_Thread[MAX_MT_THREADS];
	MT_Data_Info_JincResizeMT MT_Data[MAX_MT_THREADS];
	uint8_t threads,threads_number;
	bool sleep;
	uint32_t UserId;
	
	ThreadPoolFunction Jinc_MT;

	static void StaticThreadpool(void *ptr);

	void ProcessFrameMT(MT_Data_Info_JincResizeMT *MT_DataGF, uint8_t idxFn);

	void FreeData(void);

public:
	JincResizeMT(PClip _child, int target_width, int target_height, double crop_left, double crop_top, double crop_width, double crop_height,
		int quant_x, int quant_y, int tap, double blur, const char *_cplace, uint8_t _threads, int opt, int initial_capacity, bool initial_capacity_def,
		double initial_factor, int range, bool _sleep, bool negativePrefetch,IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment *env);
    virtual ~JincResizeMT();
	int __stdcall SetCacheHints(int cachehints, int frame_range);
};

#endif
