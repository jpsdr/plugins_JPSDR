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
// Intrinsics for SSE4.1, SSSE3, SSE3, SSE2, ISSE and MMX
#include <emmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <algorithm>

#if _MSC_VER >= 1900
  #define AVX2_BUILD_POSSIBLE
#endif

#ifdef AVX2_BUILD_POSSIBLE
#include "./resample_avx2.h"
#endif

extern ThreadPoolInterface *poolInterface;


/***************************************
 ********* Templated SSE Loader ********
 ***************************************/

typedef __m128i (SSELoader)(const __m128i*);
typedef __m128 (SSELoader_ps)(const float*);

__forceinline __m128i simd_load_aligned(const __m128i* adr)
{
  return _mm_load_si128(adr);
}

__forceinline __m128i simd_load_unaligned(const __m128i* adr)
{
  return _mm_loadu_si128(adr);
}

#if defined(CLANG)
__attribute__((__target__("sse3")))
#endif
__forceinline __m128i simd_load_unaligned_sse3(const __m128i* adr)
{
  return _mm_lddqu_si128(adr);
}

#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
__forceinline __m128i simd_load_streaming(const __m128i* adr)
{
  return _mm_stream_load_si128(const_cast<__m128i*>(adr));
}

// float loaders
__forceinline __m128 simd_loadps_aligned(const float * adr)
{
  return _mm_load_ps(adr);
}

__forceinline __m128 simd_loadps_unaligned(const float* adr)
{
  return _mm_loadu_ps(adr);
}

// useful SIMD helpers

// sse2 replacement of _mm_mullo_epi32 in SSE4.1
// use it after speed test, may have too much overhead and C is faster
__forceinline __m128i _MM_MULLO_EPI32(const __m128i &a, const __m128i &b)
{
  // for SSE 4.1: return _mm_mullo_epi32(a, b);
  __m128i tmp1 = _mm_mul_epu32(a,b); // mul 2,0
  __m128i tmp2 = _mm_mul_epu32( _mm_srli_si128(a,4), _mm_srli_si128(b,4)); // mul 3,1
  // shuffle results to [63..0] and pack. a2->a1, a0->a0
  return _mm_unpacklo_epi32(_mm_shuffle_epi32(tmp1, _MM_SHUFFLE (0,0,2,0)), _mm_shuffle_epi32(tmp2, _MM_SHUFFLE (0,0,2,0)));
}

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4309)
#endif
// fake _mm_packus_epi32 (orig is SSE4.1 only)
__forceinline __m128i _MM_PACKUS_EPI32( __m128i a, __m128i b )
{
  const static __m128i val_32 = _mm_set1_epi32(0x8000);
  const static __m128i val_16 = _mm_set1_epi16(0x8000);

  a = _mm_sub_epi32(a, val_32);
  b = _mm_sub_epi32(b, val_32);
  a = _mm_packs_epi32(a, b);
  a = _mm_add_epi16(a, val_16);
  return a;
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif
// fake _mm_packus_epi32 (orig is SSE4.1 only)
// only for packing 00000000..0000FFFF range integers, does not clamp properly above that, e.g. 00010001
__forceinline __m128i _MM_PACKUS_EPI32_SRC_TRUEWORD(__m128i a, __m128i b)
{
  a = _mm_slli_epi32 (a, 16);
  a = _mm_srai_epi32 (a, 16);
  b = _mm_slli_epi32 (b, 16);
  b = _mm_srai_epi32 (b, 16);
  a = _mm_packs_epi32 (a, b);
  return a;
}

__forceinline __m128i _MM_CMPLE_EPU16(__m128i x, __m128i y)
{
  // Returns 0xFFFF where x <= y:
  return _mm_cmpeq_epi16(_mm_subs_epu16(x, y), _mm_setzero_si128());
}

__forceinline __m128i _MM_BLENDV_SI128(__m128i x, __m128i y, __m128i mask)
{
  // Replace bit in x with bit in y when matching bit in mask is set:
  return _mm_or_si128(_mm_andnot_si128(mask, x), _mm_and_si128(mask, y));
}

// sse2 simulation of SSE4's _mm_min_epu16
__forceinline __m128i _MM_MIN_EPU16(__m128i x, __m128i y)
{
  // Returns x where x <= y, else y:
  return _MM_BLENDV_SI128(y, x, _MM_CMPLE_EPU16(x, y));
}

// sse2 simulation of SSE4's _mm_max_epu16
__forceinline __m128i _MM_MAX_EPU16(__m128i x, __m128i y)
{
  // Returns x where x >= y, else y:
  return _MM_BLENDV_SI128(x, y, _MM_CMPLE_EPU16(x, y));
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

				for (int i = 0; i < filter_size; i++)
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

				for (int i = 0; i < filter_size; i++)
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

      for (int i = 0; i < filter_size; i++)
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

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

			result = (result+Offset) >> FPScale16bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst0[x] = (uint16_t) result;
		}
		dst0 += dst_pitch;
		current_coeff += filter_size;
	}
}


#ifdef X86_32
static void resize_v_mmx_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3;
  const int sizeMod2 = (filter_size >> 1) << 1;
  const bool notMod2 = sizeMod2 < filter_size;

  const int Offset = 1 << (FPScale8bits-1);

  const __m64 zero = _mm_setzero_si64();

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m64 val_min_m64 = _mm_set1_pi16((short)((val_min << 8)|val_min));
	const __m64 val_max_m64 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_pi16((short)(((int)240 << 8)|235)) : _mm_set1_pi16((short)((val_max << 8)|val_max));

  
  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x += 8)
	{
      __m64 result_1 = _mm_set1_pi32(Offset); // Init. with rounder (16384/2 = 8192)
      __m64 result_2 = result_1;
      __m64 result_3 = result_1;
      __m64 result_4 = result_1;

      for (int i = 0; i < sizeMod2; i += 2)
	  {
        __m64 src_p1 = *(reinterpret_cast<const __m64*>(src_ptr+pitch_table[i]+x));   // For detailed explanation please see SSE2 version.
        __m64 src_p2 = *(reinterpret_cast<const __m64*>(src_ptr+pitch_table[i+1]+x));

        __m64 src_l = _mm_unpacklo_pi8(src_p1, src_p2);                                   
        __m64 src_h = _mm_unpackhi_pi8(src_p1, src_p2);                                   

        __m64 src_1 = _mm_unpacklo_pi8(src_l, zero);                                      
        __m64 src_2 = _mm_unpackhi_pi8(src_l, zero);                                      
        __m64 src_3 = _mm_unpacklo_pi8(src_h, zero);                                      
        __m64 src_4 = _mm_unpackhi_pi8(src_h, zero);                                      

        __m64 coeff = _mm_cvtsi32_si64(*reinterpret_cast<const int*>(current_coeff+i));   
        coeff = _mm_unpacklo_pi32(coeff, coeff);                                               

        __m64 dst_1 = _mm_madd_pi16(src_1, coeff);                                        
        __m64 dst_2 = _mm_madd_pi16(src_2, coeff);                                        
        __m64 dst_3 = _mm_madd_pi16(src_3, coeff);
        __m64 dst_4 = _mm_madd_pi16(src_4, coeff);

        result_1 = _mm_add_pi32(result_1, dst_1);
        result_2 = _mm_add_pi32(result_2, dst_2);
        result_3 = _mm_add_pi32(result_3, dst_3);
        result_4 = _mm_add_pi32(result_4, dst_4);
      }

      if (notMod2)
	  { // do last odd row
        __m64 src_p = *(reinterpret_cast<const __m64*>(src_ptr+pitch_table[sizeMod2]+x));

        __m64 src_l = _mm_unpacklo_pi8(src_p, zero);
        __m64 src_h = _mm_unpackhi_pi8(src_p, zero);

        __m64 coeff = _mm_set1_pi16(current_coeff[sizeMod2]);

        __m64 dst_ll = _mm_mullo_pi16(src_l, coeff);   // Multiply by coefficient
        __m64 dst_lh = _mm_mulhi_pi16(src_l, coeff);
        __m64 dst_hl = _mm_mullo_pi16(src_h, coeff);
        __m64 dst_hh = _mm_mulhi_pi16(src_h, coeff);

        __m64 dst_1 = _mm_unpacklo_pi16(dst_ll, dst_lh); // Unpack to 32-bit integer
        __m64 dst_2 = _mm_unpackhi_pi16(dst_ll, dst_lh);
        __m64 dst_3 = _mm_unpacklo_pi16(dst_hl, dst_hh);
        __m64 dst_4 = _mm_unpackhi_pi16(dst_hl, dst_hh);

        result_1 = _mm_add_pi32(result_1, dst_1);
        result_2 = _mm_add_pi32(result_2, dst_2);
        result_3 = _mm_add_pi32(result_3, dst_3);
        result_4 = _mm_add_pi32(result_4, dst_4);
      }

      // Divide by 16348 (FPRound)
      result_1  = _mm_srai_pi32(result_1, FPScale8bits);
      result_2  = _mm_srai_pi32(result_2, FPScale8bits);
      result_3  = _mm_srai_pi32(result_3, FPScale8bits);
      result_4  = _mm_srai_pi32(result_4, FPScale8bits);

      // Pack and store
      __m64 result_l = _mm_packs_pi32(result_1, result_2);
      __m64 result_h = _mm_packs_pi32(result_3, result_4);
      __m64 result   = _mm_packs_pu16(result_l, result_h);

		result = _mm_max_pu8(result,val_min_m64);
		result = _mm_min_pu8(result,val_max_m64);

      *(reinterpret_cast<__m64*>(dst+x)) = result;
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod8; x < width; x++)
		{
			int result = 0;

		    for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod8; x < width; x++)
		{
			int result = 0;

		    for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	dst += dst_pitch;
	current_coeff += filter_size;

  }

  _mm_empty();
}
#endif

template<SSELoader load>
static void resize_v_sse2_planarT(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;
  const int sizeMod2 = (filter_size >> 1) << 1;
  const bool notMod2 = sizeMod2 < filter_size;
  const int Offset = 1 << (FPScale8bits-1);

  const __m128i zero = _mm_setzero_si128();


	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x += 16)
	{
      __m128i result_1 = _mm_set1_epi32(Offset); // Init. with rounder (16384/2 = 8192)
      __m128i result_2 = result_1;
      __m128i result_3 = result_1;
      __m128i result_4 = result_1;
      
      for (int i = 0; i < sizeMod2; i += 2)
	  {
        __m128i src_p1 = load(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[i]+x));   // p|o|n|m|l|k|j|i|h|g|f|e|d|c|b|a
        __m128i src_p2 = load(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[i+1]+x)); // P|O|N|M|L|K|J|I|H|G|F|E|D|C|B|A
         
        __m128i src_l = _mm_unpacklo_epi8(src_p1, src_p2);                                   // Hh|Gg|Ff|Ee|Dd|Cc|Bb|Aa
        __m128i src_h = _mm_unpackhi_epi8(src_p1, src_p2);                                   // Pp|Oo|Nn|Mm|Ll|Kk|Jj|Ii

        __m128i src_1 = _mm_unpacklo_epi8(src_l, zero);                                      // .D|.d|.C|.c|.B|.b|.A|.a
        __m128i src_2 = _mm_unpackhi_epi8(src_l, zero);                                      // .H|.h|.G|.g|.F|.f|.E|.e
        __m128i src_3 = _mm_unpacklo_epi8(src_h, zero);                                      // etc.
        __m128i src_4 = _mm_unpackhi_epi8(src_h, zero);                                      // etc.

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));   // XX|XX|XX|XX|XX|XX|CO|co
        coeff = _mm_shuffle_epi32(coeff, 0);                                                 // CO|co|CO|co|CO|co|CO|co
        
        __m128i dst_1 = _mm_madd_epi16(src_1, coeff);                                         // CO*D+co*d | CO*C+co*c | CO*B+co*b | CO*A+co*a
        __m128i dst_2 = _mm_madd_epi16(src_2, coeff);                                         // etc.
        __m128i dst_3 = _mm_madd_epi16(src_3, coeff);
        __m128i dst_4 = _mm_madd_epi16(src_4, coeff);

        result_1 = _mm_add_epi32(result_1, dst_1);
        result_2 = _mm_add_epi32(result_2, dst_2);
        result_3 = _mm_add_epi32(result_3, dst_3);
        result_4 = _mm_add_epi32(result_4, dst_4);
      }
      
      if (notMod2)
	  { // do last odd row
        __m128i src_p = load(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[sizeMod2]+x));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        __m128i coeff = _mm_set1_epi16(current_coeff[sizeMod2]);

        __m128i dst_ll = _mm_mullo_epi16(src_l, coeff);   // Multiply by coefficient
        __m128i dst_lh = _mm_mulhi_epi16(src_l, coeff);
        __m128i dst_hl = _mm_mullo_epi16(src_h, coeff);
        __m128i dst_hh = _mm_mulhi_epi16(src_h, coeff);

        __m128i dst_1 = _mm_unpacklo_epi16(dst_ll, dst_lh); // Unpack to 32-bit integer
        __m128i dst_2 = _mm_unpackhi_epi16(dst_ll, dst_lh);
        __m128i dst_3 = _mm_unpacklo_epi16(dst_hl, dst_hh);
        __m128i dst_4 = _mm_unpackhi_epi16(dst_hl, dst_hh);

        result_1 = _mm_add_epi32(result_1, dst_1);
        result_2 = _mm_add_epi32(result_2, dst_2);
        result_3 = _mm_add_epi32(result_3, dst_3);
        result_4 = _mm_add_epi32(result_4, dst_4);
      }
      
      // Divide by 16348 (FPRound)
      result_1  = _mm_srai_epi32(result_1, FPScale8bits);
      result_2  = _mm_srai_epi32(result_2, FPScale8bits);
      result_3  = _mm_srai_epi32(result_3, FPScale8bits);
      result_4  = _mm_srai_epi32(result_4, FPScale8bits);

      // Pack and store
      __m128i result_l = _mm_packs_epi32(result_1, result_2);
      __m128i result_h = _mm_packs_epi32(result_3, result_4);
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}


#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
static void resize_v_sse2_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;
  const int sizeMod2 = (filter_size >> 1) << 1;
  const bool notMod2 = sizeMod2 < filter_size;
  const int Offset = 1 << (FPScale8bits-1);

  const __m128i zero = _mm_setzero_si128();


	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x += 16)
	{
      __m128i result_1 = _mm_set1_epi32(Offset); // Init. with rounder (16384/2 = 8192)
      __m128i result_2 = result_1;
      __m128i result_3 = result_1;
      __m128i result_4 = result_1;
      
      for (int i = 0; i < sizeMod2; i += 2)
	  {
        __m128i src_p1 = _mm_stream_load_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[i]+x)));   // p|o|n|m|l|k|j|i|h|g|f|e|d|c|b|a
        __m128i src_p2 = _mm_stream_load_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[i+1]+x))); // P|O|N|M|L|K|J|I|H|G|F|E|D|C|B|A
         
        __m128i src_l = _mm_unpacklo_epi8(src_p1, src_p2);                                   // Hh|Gg|Ff|Ee|Dd|Cc|Bb|Aa
        __m128i src_h = _mm_unpackhi_epi8(src_p1, src_p2);                                   // Pp|Oo|Nn|Mm|Ll|Kk|Jj|Ii

        __m128i src_1 = _mm_unpacklo_epi8(src_l, zero);                                      // .D|.d|.C|.c|.B|.b|.A|.a
        __m128i src_2 = _mm_unpackhi_epi8(src_l, zero);                                      // .H|.h|.G|.g|.F|.f|.E|.e
        __m128i src_3 = _mm_unpacklo_epi8(src_h, zero);                                      // etc.
        __m128i src_4 = _mm_unpackhi_epi8(src_h, zero);                                      // etc.

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));   // XX|XX|XX|XX|XX|XX|CO|co
        coeff = _mm_shuffle_epi32(coeff, 0);                                                 // CO|co|CO|co|CO|co|CO|co
        
        __m128i dst_1 = _mm_madd_epi16(src_1, coeff);                                         // CO*D+co*d | CO*C+co*c | CO*B+co*b | CO*A+co*a
        __m128i dst_2 = _mm_madd_epi16(src_2, coeff);                                         // etc.
        __m128i dst_3 = _mm_madd_epi16(src_3, coeff);
        __m128i dst_4 = _mm_madd_epi16(src_4, coeff);

        result_1 = _mm_add_epi32(result_1, dst_1);
        result_2 = _mm_add_epi32(result_2, dst_2);
        result_3 = _mm_add_epi32(result_3, dst_3);
        result_4 = _mm_add_epi32(result_4, dst_4);
      }
      
      if (notMod2)
	  { // do last odd row
        __m128i src_p = _mm_stream_load_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[sizeMod2]+x)));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        __m128i coeff = _mm_set1_epi16(current_coeff[sizeMod2]);

        __m128i dst_ll = _mm_mullo_epi16(src_l, coeff);   // Multiply by coefficient
        __m128i dst_lh = _mm_mulhi_epi16(src_l, coeff);
        __m128i dst_hl = _mm_mullo_epi16(src_h, coeff);
        __m128i dst_hh = _mm_mulhi_epi16(src_h, coeff);

        __m128i dst_1 = _mm_unpacklo_epi16(dst_ll, dst_lh); // Unpack to 32-bit integer
        __m128i dst_2 = _mm_unpackhi_epi16(dst_ll, dst_lh);
        __m128i dst_3 = _mm_unpacklo_epi16(dst_hl, dst_hh);
        __m128i dst_4 = _mm_unpackhi_epi16(dst_hl, dst_hh);

        result_1 = _mm_add_epi32(result_1, dst_1);
        result_2 = _mm_add_epi32(result_2, dst_2);
        result_3 = _mm_add_epi32(result_3, dst_3);
        result_4 = _mm_add_epi32(result_4, dst_4);
      }
      
      // Divide by 16348 (FPRound)
      result_1  = _mm_srai_epi32(result_1, FPScale8bits);
      result_2  = _mm_srai_epi32(result_2, FPScale8bits);
      result_3  = _mm_srai_epi32(result_3, FPScale8bits);
      result_4  = _mm_srai_epi32(result_4, FPScale8bits);

      // Pack and store
      __m128i result_l = _mm_packs_epi32(result_1, result_2);
      __m128i result_h = _mm_packs_epi32(result_3, result_4);
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}


template<SSELoader load>
#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void resize_v_ssse3_planarT(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;
  const int Offset = 1 << (FPScale8bits-1);

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  const __m128i zero = _mm_setzero_si128();
  const __m128i coeff_unpacker = _mm_set_epi8(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0);

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE* src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x+=16)
	{
      __m128i result_l = _mm_set1_epi16(32); // Init. with rounder ((1 << 6)/2 = 32)
      __m128i result_h = result_l;

      const BYTE* src2_ptr = src_ptr+x;
      
      for (int i = 0; i < filter_size; i++)
	  {
        __m128i src_p = load(reinterpret_cast<const __m128i*>(src2_ptr));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        src_l = _mm_slli_epi16(src_l, 7);
        src_h = _mm_slli_epi16(src_h, 7);

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));
                coeff = _mm_shuffle_epi8(coeff, coeff_unpacker);

        __m128i dst_l = _mm_mulhrs_epi16(src_l, coeff);   // Multiply by coefficient (SSSE3)
        __m128i dst_h = _mm_mulhrs_epi16(src_h, coeff);

        result_l = _mm_add_epi16(result_l, dst_l);
        result_h = _mm_add_epi16(result_h, dst_h);

        src2_ptr += src_pitch;
      }

      // Divide by 64
      result_l  = _mm_srai_epi16(result_l, 6);
      result_h  = _mm_srai_epi16(result_h, 6);

      // Pack and store
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}


#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void resize_v_ssse3_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;
  const int Offset = 1 << (FPScale8bits-1);

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  const __m128i zero = _mm_setzero_si128();
  const __m128i coeff_unpacker = _mm_set_epi8(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0);

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE* src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x+=16)
	{
      __m128i result_l = _mm_set1_epi16(32); // Init. with rounder ((1 << 6)/2 = 32)
      __m128i result_h = result_l;

      const BYTE* src2_ptr = src_ptr+x;
      
      for (int i = 0; i < filter_size; i++)
	  {
        __m128i src_p = _mm_lddqu_si128(reinterpret_cast<const __m128i*>(src2_ptr));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        src_l = _mm_slli_epi16(src_l, 7);
        src_h = _mm_slli_epi16(src_h, 7);

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));
                coeff = _mm_shuffle_epi8(coeff, coeff_unpacker);

        __m128i dst_l = _mm_mulhrs_epi16(src_l, coeff);   // Multiply by coefficient (SSSE3)
        __m128i dst_h = _mm_mulhrs_epi16(src_h, coeff);

        result_l = _mm_add_epi16(result_l, dst_l);
        result_h = _mm_add_epi16(result_h, dst_h);

        src2_ptr += src_pitch;
      }

      // Divide by 64
      result_l  = _mm_srai_epi16(result_l, 6);
      result_h  = _mm_srai_epi16(result_h, 6);

      // Pack and store
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}


// The only difference between resize_v_sse41_planar and resize_v_ssse3_planar is the load operation
#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
static void resize_v_sse41_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;
  const int Offset = 1 << (FPScale8bits-1);

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  const __m128i zero = _mm_setzero_si128();
  const __m128i coeff_unpacker = _mm_set_epi8(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0);

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE* src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x+=16)
	{
      __m128i result_l = _mm_set1_epi16(32); // Init. with rounder ((1 << 6)/2 = 32)
      __m128i result_h = result_l;

      const BYTE* src2_ptr = src_ptr+x;
      
      for (int i = 0; i < filter_size; i++)
	  {
        __m128i src_p = _mm_stream_load_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src2_ptr)));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        src_l = _mm_slli_epi16(src_l, 7);
        src_h = _mm_slli_epi16(src_h, 7);

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));
                coeff = _mm_shuffle_epi8(coeff, coeff_unpacker);

        __m128i dst_l = _mm_mulhrs_epi16(src_l, coeff);   // Multiply by coefficient (SSSE3)
        __m128i dst_h = _mm_mulhrs_epi16(src_h, coeff);

        result_l = _mm_add_epi16(result_l, dst_l);
        result_h = _mm_add_epi16(result_h, dst_h);

        src2_ptr += src_pitch;
      }

      // Divide by 64
      result_l  = _mm_srai_epi16(result_l, 6);
      result_h  = _mm_srai_epi16(result_h, 6);

      // Pack and store
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+Offset) >> FPScale8bits;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}

//-------- 128 bit uint16_t Verticals

template<bool lessthan16bit, int index>
__forceinline static void process_chunk_v_uint16_t(const uint16_t *src2_ptr, int src_pitch, __m128i &coeff01234567, __m128i &result_single_lo, __m128i &result_single_hi, const __m128i &shifttosigned)
{
  // offset table generating is what preventing us from overaddressing
  // 0-1
  __m128i src_even = _mm_load_si128(reinterpret_cast<const __m128i*>(src2_ptr + index * src_pitch)); // 4x 16bit pixels
  __m128i src_odd = _mm_load_si128(reinterpret_cast<const __m128i*>(src2_ptr + (index + 1) * src_pitch));  // 4x 16bit pixels
  __m128i src_lo = _mm_unpacklo_epi16(src_even, src_odd);
  __m128i src_hi = _mm_unpackhi_epi16(src_even, src_odd);
  if (!lessthan16bit)
  {
    src_lo = _mm_add_epi16(src_lo, shifttosigned);
    src_hi = _mm_add_epi16(src_hi, shifttosigned);
  }
  __m128i coeff = _mm_shuffle_epi32(coeff01234567, (index / 2) | ((index / 2) << 2) | ((index / 2) << 4) | ((index / 2) << 6)); // spread short pair
  result_single_lo = _mm_add_epi32(result_single_lo, _mm_madd_epi16(src_lo, coeff)); // a*b + c
  result_single_hi = _mm_add_epi32(result_single_hi, _mm_madd_epi16(src_hi, coeff)); // a*b + c
}

// program->filtersize: 1..16 special optimized, >8: normal
template<bool lessthan16bit, int _filter_size_numOfFullBlk8, int filtersizemod8>
#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
void internal_resize_v_sse41_planar_uint16_t(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size_numOfFullBlk8 = (_filter_size_numOfFullBlk8 >= 0) ? _filter_size_numOfFullBlk8 : (program->filter_size >> 3);
  const short *current_coeff = program->pixel_coefficient + program->filter_size*MinY;
  const int filter_size_numOfFullBlk8_8 = filter_size_numOfFullBlk8 << 3;

#define NON32_BYTES_ALIGNMENT
  // in AVS+ 32 bytes alignment is guaranteed
#ifdef NON32_BYTES_ALIGNMENT
  int wMod8 = (width >> 3) << 3; // uint16: 8 at a time (2x128bit)
#endif

  const __m128i zero = _mm_setzero_si128();
  const __m128i shifttosigned = _mm_set1_epi16(-32768); // for 16 bits only
  const __m128i shiftfromsigned = _mm_set1_epi32(32768 << FPScale16bits); // for 16 bits only
  const __m128i rounder = _mm_set1_epi32(1 << (FPScale16bits - 1));

  const uint16_t* src = (uint16_t *)src0;
  uint16_t* dst = (uint16_t *)dst0;
  dst_pitch >>= 1;
  src_pitch >>= 1;
  const int src_pitch8 = src_pitch << 3;

	const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
		((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));

	const int64_t Offset = 1 << (FPScale16bits - 1); // rounder

	__m128i clamp_limit_min = _mm_set1_epi16(val_min);
	__m128i clamp_limit_max = _mm_set1_epi16(val_max);

  for (int y = MinY; y < MaxY; y++) 
  {
    int offset = program->pixel_offset[y];
	const uint16_t* src_ptr = src + pitch_table[offset];

#ifdef NON32_BYTES_ALIGNMENT
    for (int x = 0; x < wMod8; x += 8)
	{ // 2x4 pixels at a time
#else
    for (int x = 0; x < width; x += 8)
	{
#endif
      __m128i result_single_lo = rounder;
      __m128i result_single_hi = rounder;

      const uint16_t* src2_ptr = src_ptr + x;

      for (int i = 0; i < filter_size_numOfFullBlk8; i++)
	  {
        __m128i coeff01234567 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(current_coeff + (i << 3))); // 4x (2x16bit) shorts for even/odd

        // offset table generating is what preventing us from overaddressing
        // 0-1
        process_chunk_v_uint16_t<lessthan16bit, 0>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        // 2-3
        process_chunk_v_uint16_t<lessthan16bit, 2>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        // 4-5
        process_chunk_v_uint16_t<lessthan16bit, 4>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        // 6-7
        process_chunk_v_uint16_t<lessthan16bit, 6>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        src2_ptr += src_pitch8;
      }

      // and the rest non-div8 chunk
      __m128i coeff01234567 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(current_coeff + filter_size_numOfFullBlk8_8)); // 4x (2x16bit) shorts for even/odd
      if (filtersizemod8 >= 2)
        process_chunk_v_uint16_t<lessthan16bit, 0>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
      if (filtersizemod8 >= 4)
        process_chunk_v_uint16_t<lessthan16bit, 2>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
      if (filtersizemod8 >= 6)
        process_chunk_v_uint16_t<lessthan16bit, 4>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
      if ((filtersizemod8 & 1)!=0)
	  { // remaining odd one
        const int index = filtersizemod8 - 1;
        __m128i src_even = _mm_load_si128(reinterpret_cast<const __m128i*>(src2_ptr + index * src_pitch)); // 8x 16bit pixels
        if (!lessthan16bit)
          src_even = _mm_add_epi16(src_even, shifttosigned);
        __m128i coeff = _mm_shuffle_epi32(coeff01234567, (index / 2) | ((index / 2) << 2) | ((index / 2) << 4) | ((index / 2) << 6));
        __m128i src_lo = _mm_unpacklo_epi16(src_even, zero); // insert zero after the unsigned->signed shift!
        __m128i src_hi = _mm_unpackhi_epi16(src_even, zero); // insert zero after the unsigned->signed shift!
        result_single_lo = _mm_add_epi32(result_single_lo, _mm_madd_epi16(src_lo, coeff)); // a*b + c
        result_single_hi = _mm_add_epi32(result_single_hi, _mm_madd_epi16(src_hi, coeff)); // a*b + c
      }

      // correct if signed, scale back, store
      __m128i result_lo = result_single_lo;
      __m128i result_hi = result_single_hi;
      if (!lessthan16bit)
	  {
        result_lo = _mm_add_epi32(result_lo, shiftfromsigned);
        result_hi = _mm_add_epi32(result_hi, shiftfromsigned);
      }
      result_lo = _mm_srai_epi32(result_lo, FPScale16bits); // shift back integer arithmetic 13 bits precision
      result_hi = _mm_srai_epi32(result_hi, FPScale16bits);

      __m128i result_8x_uint16 = _mm_packus_epi32(result_lo, result_hi);
      //if (lessthan16bit)
        result_8x_uint16 = _mm_min_epu16(result_8x_uint16,clamp_limit_max); // extra clamp for 10-14 bit
		result_8x_uint16 = _mm_max_epu16(result_8x_uint16,clamp_limit_min); // extra clamp for 10-14 bit
      _mm_store_si128(reinterpret_cast<__m128i *>(dst + x), result_8x_uint16);
    }

#ifdef NON32_BYTES_ALIGNMENT
    // Leftover, slow C
    for (int x = wMod8; x < width; x++)
	{
      int64_t result64 = Offset; // rounder
      const uint16_t* src2_ptr = src_ptr + x;
      for (int i = 0; i < program->filter_size; i++)
	  {
        //result64 += (src_ptr + pitch_table[i] / sizeof(uint16_t))[x] * (int64_t)current_coeff[i];
        result64 += (int)(*src2_ptr) * (int64_t)current_coeff[i];
        src2_ptr += src_pitch;
      }
      int result = (int)(result64 >> FPScale16bits); // scale back 13 bits
      result = result > val_max ? val_max : result < val_min ? val_min : result; // clamp 10..16 bits
      dst[x] = (uint16_t)result;
    }
#endif

    dst += dst_pitch;
    current_coeff += program->filter_size;
  }

#undef NON32_BYTES_ALIGNMENT
}


// program->filtersize: 1..16 special optimized, >8: normal
template<bool lessthan16bit, int _filter_size_numOfFullBlk8, int filtersizemod8>
void internal_resize_v_sse2_planar_uint16_t(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size_numOfFullBlk8 = (_filter_size_numOfFullBlk8 >= 0) ? _filter_size_numOfFullBlk8 : (program->filter_size >> 3);
  const short *current_coeff = program->pixel_coefficient + program->filter_size*MinY;
  const int filter_size_numOfFullBlk8_8 = filter_size_numOfFullBlk8 << 3;

#define NON32_BYTES_ALIGNMENT
  // in AVS+ 32 bytes alignment is guaranteed
#ifdef NON32_BYTES_ALIGNMENT
  int wMod8 = (width >> 3) << 3; // uint16: 8 at a time (2x128bit)
#endif

  const __m128i zero = _mm_setzero_si128();
  const __m128i shifttosigned = _mm_set1_epi16(-32768); // for 16 bits only
  const __m128i shiftfromsigned = _mm_set1_epi32(32768 << FPScale16bits); // for 16 bits only
  const __m128i rounder = _mm_set1_epi32(1 << (FPScale16bits - 1));

  const uint16_t* src = (uint16_t *)src0;
  uint16_t* dst = (uint16_t *)dst0;
  dst_pitch >>= 1;
  src_pitch >>= 1;
  const int src_pitch8 = src_pitch << 3;

	const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
		((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));

	const int64_t Offset = 1 << (FPScale16bits - 1); // rounder

	__m128i clamp_limit_min = _mm_set1_epi16(val_min);
	__m128i clamp_limit_max = _mm_set1_epi16(val_max);

  for (int y = MinY; y < MaxY; y++) 
  {
    int offset = program->pixel_offset[y];
	const uint16_t* src_ptr = src + pitch_table[offset];

#ifdef NON32_BYTES_ALIGNMENT
    for (int x = 0; x < wMod8; x += 8)
	{ // 2x4 pixels at a time
#else
    for (int x = 0; x < width; x += 8)
	{
#endif
      __m128i result_single_lo = rounder;
      __m128i result_single_hi = rounder;

      const uint16_t* src2_ptr = src_ptr + x;

      for (int i = 0; i < filter_size_numOfFullBlk8; i++)
	  {
        __m128i coeff01234567 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(current_coeff + (i << 3))); // 4x (2x16bit) shorts for even/odd

        // offset table generating is what preventing us from overaddressing
        // 0-1
        process_chunk_v_uint16_t<lessthan16bit, 0>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        // 2-3
        process_chunk_v_uint16_t<lessthan16bit, 2>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        // 4-5
        process_chunk_v_uint16_t<lessthan16bit, 4>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        // 6-7
        process_chunk_v_uint16_t<lessthan16bit, 6>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
        src2_ptr += src_pitch8;
      }

      // and the rest non-div8 chunk
      __m128i coeff01234567 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(current_coeff + filter_size_numOfFullBlk8_8)); // 4x (2x16bit) shorts for even/odd
      if (filtersizemod8 >= 2)
        process_chunk_v_uint16_t<lessthan16bit, 0>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
      if (filtersizemod8 >= 4)
        process_chunk_v_uint16_t<lessthan16bit, 2>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
      if (filtersizemod8 >= 6)
        process_chunk_v_uint16_t<lessthan16bit, 4>(src2_ptr, src_pitch, coeff01234567, result_single_lo, result_single_hi, shifttosigned);
      if ((filtersizemod8 & 1)!=0)
	  { // remaining odd one
        const int index = filtersizemod8 - 1;
        __m128i src_even = _mm_load_si128(reinterpret_cast<const __m128i*>(src2_ptr + index * src_pitch)); // 8x 16bit pixels
        if (!lessthan16bit)
          src_even = _mm_add_epi16(src_even, shifttosigned);
        __m128i coeff = _mm_shuffle_epi32(coeff01234567, (index / 2) | ((index / 2) << 2) | ((index / 2) << 4) | ((index / 2) << 6));
        __m128i src_lo = _mm_unpacklo_epi16(src_even, zero); // insert zero after the unsigned->signed shift!
        __m128i src_hi = _mm_unpackhi_epi16(src_even, zero); // insert zero after the unsigned->signed shift!
        result_single_lo = _mm_add_epi32(result_single_lo, _mm_madd_epi16(src_lo, coeff)); // a*b + c
        result_single_hi = _mm_add_epi32(result_single_hi, _mm_madd_epi16(src_hi, coeff)); // a*b + c
      }

      // correct if signed, scale back, store
      __m128i result_lo = result_single_lo;
      __m128i result_hi = result_single_hi;
      if (!lessthan16bit)
	  {
        result_lo = _mm_add_epi32(result_lo, shiftfromsigned);
        result_hi = _mm_add_epi32(result_hi, shiftfromsigned);
      }
      result_lo = _mm_srai_epi32(result_lo, FPScale16bits); // shift back integer arithmetic 13 bits precision
      result_hi = _mm_srai_epi32(result_hi, FPScale16bits);

      __m128i result_8x_uint16 = _MM_PACKUS_EPI32(result_lo, result_hi);
      //if (lessthan16bit)
        result_8x_uint16 = _MM_MIN_EPU16(result_8x_uint16,clamp_limit_max); // extra clamp for 10-14 bit
		result_8x_uint16 = _MM_MAX_EPU16(result_8x_uint16,clamp_limit_min); // extra clamp for 10-14 bit
      _mm_store_si128(reinterpret_cast<__m128i *>(dst + x), result_8x_uint16);
    }

#ifdef NON32_BYTES_ALIGNMENT
    // Leftover, slow C
    for (int x = wMod8; x < width; x++)
	{
      int64_t result64 = Offset; // rounder
      const uint16_t* src2_ptr = src_ptr + x;
      for (int i = 0; i < program->filter_size; i++)
	  {
        //result64 += (src_ptr + pitch_table[i] / sizeof(uint16_t))[x] * (int64_t)current_coeff[i];
        result64 += (int)(*src2_ptr) * (int64_t)current_coeff[i];
        src2_ptr += src_pitch;
      }
      int result = (int)(result64 >> FPScale16bits); // scale back 13 bits
      result = result > val_max ? val_max : result < val_min ? val_min : result; // clamp 10..16 bits
      dst[x] = (uint16_t)result;
    }
#endif

    dst += dst_pitch;
    current_coeff += program->filter_size;
  }

#undef NON32_BYTES_ALIGNMENT
}


//-------- uint16_t Vertical Dispatcher

template<bool lessthan16bit>
void resize_v_sse2_planar_uint16_t(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  // template<bool lessthan16bit, int _filter_size_numOfFullBlk8, int filtersizemod8>
  // filtersize 1..16: to template for optimization
  switch (program->filter_size)
  {
  case 1:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 2:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 3:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 4:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 5:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 6:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 7:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 0, 7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 8:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 0>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 9:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 10:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 11:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 12:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 13:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 14:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 15:
    internal_resize_v_sse2_planar_uint16_t<lessthan16bit, 1, 7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  default:
    switch (program->filter_size & 7)
	{
    case 0:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 0>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 1:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 2:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 3:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 4:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 5:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 6:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 7:
      internal_resize_v_sse2_planar_uint16_t<lessthan16bit, -1, 7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    }
    break;
  }
}


template<bool lessthan16bit>
#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
void resize_v_sse41_planar_uint16_t(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  // template<bool lessthan16bit, int _filter_size_numOfFullBlk8, int filtersizemod8>
  // filtersize 1..16: to template for optimization
  switch (program->filter_size)
  {
  case 1:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 2:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 3:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 4:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 5:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 6:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 7:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 0, 7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 8:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 0>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 9:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 10:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 11:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 12:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 13:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 14:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 15:
    internal_resize_v_sse41_planar_uint16_t<lessthan16bit, 1, 7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  default:
    switch (program->filter_size & 7)
	{
    case 0:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 0>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 1:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 2:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 3:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 4:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 5:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 6:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    case 7:
      internal_resize_v_sse41_planar_uint16_t<lessthan16bit, -1, 7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
      break;
    }
    break;
  }
}


//-------- 128 bit float Verticals

template<int _filtersize>
static void internal_resize_v_sse2_planar_float(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  // 1..8: special case for compiler optimization
  const int filter_size = _filtersize >= 1 ? _filtersize : program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float + filter_size*MinY;

  const float* src = (float *)src0;
  float* dst = (float *)dst0;
  dst_pitch >>= 2;
  src_pitch >>= 2;
  const int src_pitch2 = src_pitch << 1;
  const int src_pitch3 = src_pitch*3;
  const int src_pitch4 = src_pitch << 2;

  const int fsmod4 = (filter_size >> 2) << 2;
  for (int y = MinY; y < MaxY; y++)
  {
    int offset = program->pixel_offset[y];
	const float* src_ptr = src + pitch_table[offset];

    // 8 pixels/cycle (32 bytes)
    for (int x = 0; x < width; x += 8)
	{  // safe to process 8 floats, 32 bytes alignment is OK
      __m128 result_single_lo = _mm_set1_ps(0.0f);
      __m128 result_single_hi = _mm_set1_ps(0.0f);

      const float* src2_ptr = src_ptr + x;

      for (int i = 0; i < fsmod4; i += 4)
	  {
        __m128 src_single_lo;
        __m128 src_single_hi;
        __m128 coeff0123 = _mm_loadu_ps(reinterpret_cast<const float*>(current_coeff_float + i)); // loads 4 floats
        __m128 coeff;

        // unroll 4x
        // #1
        src_single_lo = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr));
        src_single_hi = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + 4));
        coeff = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(coeff0123), (0 << 0) | (0 << 2) | (0 << 4) | (0 << 6))); // spread 0th
        result_single_lo = _mm_add_ps(result_single_lo, _mm_mul_ps(coeff, src_single_lo));
        result_single_hi = _mm_add_ps(result_single_hi, _mm_mul_ps(coeff, src_single_hi));

        // #2
        src_single_lo = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + src_pitch));
        src_single_hi = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + src_pitch + 4));
        coeff = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(coeff0123), (1 << 0) | (1 << 2) | (1 << 4) | (1 << 6))); // spread 1st
        result_single_lo = _mm_add_ps(result_single_lo, _mm_mul_ps(coeff, src_single_lo));
        result_single_hi = _mm_add_ps(result_single_hi, _mm_mul_ps(coeff, src_single_hi));

        // #3
        src_single_lo = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + src_pitch2));
        src_single_hi = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + src_pitch2 + 4));
        coeff = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(coeff0123), (2 << 0) | (2 << 2) | (2 << 4) | (2 << 6))); // spread 2nd
        result_single_lo = _mm_add_ps(result_single_lo, _mm_mul_ps(coeff, src_single_lo));
        result_single_hi = _mm_add_ps(result_single_hi, _mm_mul_ps(coeff, src_single_hi));

        // #4
        src_single_lo = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + src_pitch3));
        src_single_hi = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + src_pitch3 + 4));
        coeff = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(coeff0123), (3 << 0) | (3 << 2) | (3 << 4) | (3 << 6))); // spread 3rd
        result_single_lo = _mm_add_ps(result_single_lo, _mm_mul_ps(coeff, src_single_lo));
        result_single_hi = _mm_add_ps(result_single_hi, _mm_mul_ps(coeff, src_single_hi));

        src2_ptr += src_pitch4;
      }

      // one-by-one
      for (int i = fsmod4; i < filter_size; i++)
	  {
        __m128 src_single_lo = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr));
        __m128 src_single_hi = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr + 4));
        __m128 coeff = _mm_load1_ps(reinterpret_cast<const float*>(current_coeff_float + i)); // loads 1, fills all 8 floats
        result_single_lo = _mm_add_ps(result_single_lo, _mm_mul_ps(coeff, src_single_lo)); // _mm_fmadd_ps(src_single, coeff, result_single); // a*b + c
        result_single_hi = _mm_add_ps(result_single_hi, _mm_mul_ps(coeff, src_single_hi)); // _mm_fmadd_ps(src_single, coeff, result_single); // a*b + c

        src2_ptr += src_pitch;
      }

      _mm_stream_ps(reinterpret_cast<float*>(dst + x), result_single_lo);
      _mm_stream_ps(reinterpret_cast<float*>(dst + x + 4), result_single_hi);
    }

#if 0
    // Leftover, Slow C
    for (int x = wMod4; x < width; x++)
	{
      float result = 0;
      const float* src2_ptr = src_ptr + x;
      for (int i = 0; i < filter_size; i++)
	  {
        result += (*src2_ptr) * current_coeff_float[i];
        src2_ptr += src_pitch;
      }
      dst[x] = result;
    }
#endif
    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
}

//-------- Float Vertical Dispatcher

void resize_v_sse2_planar_float(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  // 1..8: special case for compiler optimization
  switch (program->filter_size)
  {
  case 1:
    internal_resize_v_sse2_planar_float<1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 2:
    internal_resize_v_sse2_planar_float<2>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 3:
    internal_resize_v_sse2_planar_float<3>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 4:
    internal_resize_v_sse2_planar_float<4>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 5:
    internal_resize_v_sse2_planar_float<5>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 6:
    internal_resize_v_sse2_planar_float<6>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 7:
    internal_resize_v_sse2_planar_float<7>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  case 8:
    internal_resize_v_sse2_planar_float<8>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
  default:
    internal_resize_v_sse2_planar_float<-1>(dst,src,dst_pitch,src_pitch,program,width,bits_per_pixel,MinY,MaxY,pitch_table,storage,range,mode_YUY2);
    break;
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
/*
static void resize_h_pointresize(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int wMod4 = (width >> 2) << 2;

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < wMod4; x+=4)
	{
#define pixel(a) src[program->pixel_offset[x+a]]
      unsigned int data = (pixel(3) << 24) + (pixel(2) << 16) + (pixel(1) << 8) + pixel(0);
#undef pixel
      *((unsigned int *)(dst+x)) = data;
    }

    for (int x = wMod4; x < width; x++)
      dst[x] = src[program->pixel_offset[x]];

    dst += dst_pitch;
    src += src_pitch;
  }
}
*/

// make the resampling coefficient array mod8 friendly for simd, padding non-used coeffs with zeros
static void resize_h_prepare_coeff_8or16(ResamplingProgram* p,IScriptEnvironment* env,int alignFilterSize8or16)
{
  const int im0=p->target_size;
  const int im1=p->filter_size;
	const int filter_size = AlignNumber(im1,alignFilterSize8or16);
	const int target_size = AlignNumber(im0,ALIGN_RESIZER_TARGET_SIZE);
  p->filter_size_alignment = alignFilterSize8or16;

  if (p->bits_per_pixel==32)
  {
	  float *new_coeff_float = (float *)_aligned_malloc(sizeof(float)*target_size*filter_size,64);

	  if (new_coeff_float==NULL) env->ThrowError("ResizeHMT: Could not reserve memory in a resampler.");
	  std::fill_n(new_coeff_float,target_size*filter_size,0.0f);

	  float *dst_f=new_coeff_float,*src_f=p->pixel_coefficient_float;

	  for (int i=0; i<im0; i++)
	  {
		  for (int j=0; j<im1; j++)
			  dst_f[j]=src_f[j];

		  dst_f += filter_size;
		  src_f += im1;
	  }
	   myalignedfree(p->pixel_coefficient_float);
	   p->pixel_coefficient_float = new_coeff_float;
  }
  else
  {
	  short *new_coeff = (short*)_aligned_malloc(sizeof(short)*target_size*filter_size,64);

	  if (new_coeff==NULL) env->ThrowError("ResizeHMT: Could not reserve memory in a resampler.");
	  memset(new_coeff,0,sizeof(short)*target_size*filter_size);

	  short *dst=new_coeff,*src=p->pixel_coefficient;

	  for (int i=0; i<im0; i++)
	  {
		  for (int j=0; j<im1; j++)
			  dst[j]=src[j];

		  dst += filter_size;
		  src += im1;
	  }
	  myalignedfree(p->pixel_coefficient);
	  p->pixel_coefficient = new_coeff;
  }
}


static void resize_h_c_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
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
		  
				for (int i = 0; i < filter_size; i++)
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
		  
				for (int i = 0; i < filter_size; i++)
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
		  
			for (int i = 0; i < filter_size; i++)
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
		  
		  for (int i = 0; i < filter_size; i++)
			  result+=(src0+y_src_pitch)[(begin+i)]*current_coeff[i];
		  
		  (dst0 + y_dst_pitch)[x] = result;
		  current_coeff+=filter_size;
	  }
	  y_dst_pitch+=dst_pitch;
	  y_src_pitch+=src_pitch;
  }
}


//-------- 128 bit float Horizontals

__forceinline static void process_one_pixel_h_float(const float *src, int begin, int i, float *&current_coeff, __m128 &result)
{
  // 2x4 pixels
  __m128 data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin + (i << 3)));
  __m128 data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin + (i << 3) + 4));
  __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff)); // always aligned
  __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));
  __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
  __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
  result = _mm_add_ps(result, dst_l); // accumulate result.
  result = _mm_add_ps(result, dst_h);
  current_coeff += 8;
}

template<int filtersizemod8>
__forceinline static void process_one_pixel_h_float_mask(const float *src, int begin, int i, float *&current_coeff, __m128 &result, __m128 &mask)
{
  __m128 data_l_single;
  __m128 data_h_single;
  // 2x4 pixels
  if (filtersizemod8 > 4)
  { // keep low, mask high 4 pixels
    data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin + (i << 3)));
    data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin + (i << 3) + 4));
    data_h_single = _mm_and_ps(data_h_single, mask);
  }
  else if (filtersizemod8 == 4)
  { // keep low, zero high 4 pixels
    data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin + (i << 3)));
    data_h_single = _mm_setzero_ps();
  }
  else
  { // filtersizemod8 1..3
    data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin + (i << 3)));
    data_l_single = _mm_and_ps(data_l_single, mask);
    data_h_single = _mm_setzero_ps();
  }
  __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff)); // always aligned
  __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));
  __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
  __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
  result = _mm_add_ps(result, dst_l); // accumulate result.
  result = _mm_add_ps(result, dst_h);
  current_coeff += 8;
}

// filtersizealigned8: special: 1, 2. Generic: -1
template<int filtersizealigned8, int filtersizemod8>
#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void resizer_h_ssse3_generic_float(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size_numOfBlk8 = (filtersizealigned8 >= 1) ? filtersizealigned8 : (AlignNumber(program->filter_size, 8) >> 3);

  const float *src = reinterpret_cast<const float *>(src8);
  float *dst = reinterpret_cast<float *>(dst8);
  dst_pitch >>= 2;
  src_pitch >>= 2;

  // OMG! 18.01.19
  // Protection against NaN
  // When reading the last 8 consecutive pixels from right side offsets, it would access beyond-last-pixel area.
  // One SIMD cycle reads 8 bytes from (src + begin + i * 8)
  // When program->filter_size mod 8 is 1..7 then some of the last pixels should be masked because there can be NaN garbage.
  // So it's not enough to mask the coefficients by zero. Theory: let's multiply offscreen elements by 0 which works for integer samples.
  // But we are using float, so since NaN * Zero is NaN which propagates further to NaN when hadd is summing up the pixel*coeff series

  __m128 mask;
  switch (filtersizemod8 & 3)
  {
  case 3: mask = _mm_castsi128_ps(_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0)); break; // keep 0-1-2, drop #3
  case 2: mask = _mm_castsi128_ps(_mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0, 0)); break; // keep 0-1, drop #2-3
  case 1: mask = _mm_castsi128_ps(_mm_setr_epi32(0xFFFFFFFF, 0, 0, 0)); break; // keep 0, drop #1-2-3
  default: break; // for mod4 = 0 no masking needed
  }

  const int pixels_per_cycle = 8; // doing 8 is faster than 4
  const int unsafe_limit = (program->overread_possible && filtersizemod8 != 0) ? (program->source_overread_beyond_targetx / pixels_per_cycle) * pixels_per_cycle : width;

  for (int y = 0; y < height; y++)
  {
    float* current_coeff = program->pixel_coefficient_float;

    // loop for clean, non-offscreen data
    for (int x = 0; x < unsafe_limit; x += pixels_per_cycle)
	{
      __m128 result1 = _mm_set1_ps(0.0f);
      __m128 result2 = result1;
      __m128 result3 = result1;
      __m128 result4 = result1;

      // 1-4
      int begin1 = program->pixel_offset[x + 0];
      int begin2 = program->pixel_offset[x + 1];
      int begin3 = program->pixel_offset[x + 2];
      int begin4 = program->pixel_offset[x + 3];

      // begin1, result1
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin1, i, current_coeff, result1);

      // begin2, result2
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin2, i, current_coeff, result2);

      // begin3, result3
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin3, i, current_coeff, result3);

      // begin4, result4
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin4, i, current_coeff, result4);

      // this part needs ssse3
      __m128 result12 = _mm_hadd_ps(result1, result2);
      __m128 result34 = _mm_hadd_ps(result3, result4);
      __m128 result = _mm_hadd_ps(result12, result34);

      _mm_stream_ps(reinterpret_cast<float*>(dst + x), result); // 4 results at a time

      // 5-8
      result1 = _mm_set1_ps(0.0f);
      result2 = result1;
      result3 = result1;
      result4 = result1;

      begin1 = program->pixel_offset[x + 4];
      begin2 = program->pixel_offset[x + 5];
      begin3 = program->pixel_offset[x + 6];
      begin4 = program->pixel_offset[x + 7];

      // begin1, result1
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin1, i, current_coeff, result1);

      // begin2, result2
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin2, i, current_coeff, result2);

      // begin3, result3
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin3, i, current_coeff, result3);

      // begin4, result4
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_float(src, begin4, i, current_coeff, result4);

      // this part needs ssse3
     result12 = _mm_hadd_ps(result1, result2);
     result34 = _mm_hadd_ps(result3, result4);
     result = _mm_hadd_ps(result12, result34);

      _mm_stream_ps(reinterpret_cast<float*>(dst + x + 4), result); // 4 results at a time
    } // for x

      // possibly right-side offscreen
      // and the same for the rest with masking the last filtersize/8 chunk
    for (int x = unsafe_limit; x < width; x += 4)
	{
      __m128 result1 = _mm_set1_ps(0.0f);
      __m128 result2 = result1;
      __m128 result3 = result1;
      __m128 result4 = result1;

      int begin1 = program->pixel_offset[x + 0];
      int begin2 = program->pixel_offset[x + 1];
      int begin3 = program->pixel_offset[x + 2];
      int begin4 = program->pixel_offset[x + 3];

      // begin1, result1
      for (int i = 0; i < filter_size_numOfBlk8 - 1; i++)
        process_one_pixel_h_float(src, begin1, i, current_coeff, result1);
      if (begin1 < program->source_overread_offset)
        process_one_pixel_h_float(src, begin1, filter_size_numOfBlk8 - 1, current_coeff, result1);
      else
        process_one_pixel_h_float_mask<filtersizemod8>(src, begin1, filter_size_numOfBlk8 - 1, current_coeff, result1, mask);

      // begin2, result2
      for (int i = 0; i < filter_size_numOfBlk8 - 1; i++)
        process_one_pixel_h_float(src, begin2, i, current_coeff, result2);
      if (begin2 < program->source_overread_offset)
        process_one_pixel_h_float(src, begin2, filter_size_numOfBlk8 - 1, current_coeff, result2);
      else
        process_one_pixel_h_float_mask<filtersizemod8>(src, begin2, filter_size_numOfBlk8 - 1, current_coeff, result2, mask);

      // begin3, result3
      for (int i = 0; i < filter_size_numOfBlk8 - 1; i++)
        process_one_pixel_h_float(src, begin3, i, current_coeff, result3);
      if (begin3 < program->source_overread_offset)
        process_one_pixel_h_float(src, begin3, filter_size_numOfBlk8 - 1, current_coeff, result3);
      else
        process_one_pixel_h_float_mask<filtersizemod8>(src, begin3, filter_size_numOfBlk8 - 1, current_coeff, result3, mask);

      // begin4, result4
      for (int i = 0; i < filter_size_numOfBlk8 - 1; i++)
        process_one_pixel_h_float(src, begin4, i, current_coeff, result4);
      if (begin4 < program->source_overread_offset)
        process_one_pixel_h_float(src, begin4, filter_size_numOfBlk8 - 1, current_coeff, result4);
      else
        process_one_pixel_h_float_mask<filtersizemod8>(src, begin4, filter_size_numOfBlk8 - 1, current_coeff, result4, mask);

      // this part needs ssse3
      __m128 result12 = _mm_hadd_ps(result1, result2);
      __m128 result34 = _mm_hadd_ps(result3, result4);
      __m128 result = _mm_hadd_ps(result12, result34);

      _mm_stream_ps(reinterpret_cast<float*>(dst + x), result); // 4 results at a time

    } // for x

    dst += dst_pitch;
    src += src_pitch;
  }
  /*
  // check Nans
  dst -= dst_pitch * height;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x += 4) {
      if (std::isnan(dst[x]))
      {
        x = x;
      }
    }
    dst += dst_pitch;
  }
  */
}


//-------- 128 bit uint16_t Horizontals

template<bool lessthan16bit>
__forceinline static void process_one_pixel_h_uint16_t(const uint16_t *src, int begin, int i, short *&current_coeff, __m128i &result, const __m128i &shifttosigned)
{
  __m128i data_single = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin + (i << 3))); // 8 pixels
  if (!lessthan16bit)
    data_single = _mm_add_epi16(data_single, shifttosigned); // unsigned -> signed
  __m128i coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff)); // 8 coeffs
  result = _mm_add_epi32(result, _mm_madd_epi16(data_single, coeff));
  current_coeff += 8;
}

// filter_size <= 8 -> filter_size_align8 == 1 -> no loop, hope it'll be optimized
// filter_size <= 16 -> filter_size_align8 == 2 -> loop 0..1 hope it'll be optimized
// filter_size > 16 -> use parameter AlignNumber(program->filter_size_numOfFullBlk8, 8) / 8;
template<bool lessthan16bit, int filtersizealigned8>
#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void internal_resizer_h_ssse3_generic_uint16_t(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // 1 and 2: special case for compiler optimization
  const int filter_size_numOfBlk8 = (filtersizealigned8 >= 1) ? filtersizealigned8 : (AlignNumber(program->filter_size,8) >> 3);

  const __m128i zero = _mm_setzero_si128();
  const __m128i shifttosigned = _mm_set1_epi16(-32768); // for 16 bits only
  const __m128i shiftfromsigned = _mm_set1_epi32(+32768 << FPScale16bits); // for 16 bits only
  const __m128i rounder = _mm_set_epi32(0, 0, 0, 1 << (FPScale16bits - 1)); // only once

  const uint16_t *src = reinterpret_cast<const uint16_t *>(src8);
  uint16_t *dst = reinterpret_cast<uint16_t *>(dst8);
  dst_pitch >>= 1;
  src_pitch >>= 1;

	const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
		((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));

  __m128i clamp_limit_min = _mm_set1_epi16(val_min);
  __m128i clamp_limit_max = _mm_set1_epi16(val_max);

  for (int y = 0; y < height; y++)
  {
    short* current_coeff = program->pixel_coefficient;

    for (int x = 0; x < width; x += 4)
	{
      __m128i result1 = rounder;
      __m128i result2 = result1;
      __m128i result3 = result1;
      __m128i result4 = result1;

      int begin1 = program->pixel_offset[x + 0];
      int begin2 = program->pixel_offset[x + 1];
      int begin3 = program->pixel_offset[x + 2];
      int begin4 = program->pixel_offset[x + 3];

      // this part is repeated 4 times
      // begin1, result1
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin1, i, current_coeff, result1, shifttosigned);

      // begin2, result2
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin2, i, current_coeff, result2, shifttosigned);

      // begin3, result3
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin3, i, current_coeff, result3, shifttosigned);

      // begin4, result4
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin4, i, current_coeff, result4, shifttosigned);

      const __m128i sumQuad12 = _mm_hadd_epi32(result1, result2); // L1L1L1L1 + L2L2L2L2 = L1L1 L2L2
      const __m128i sumQuad34 = _mm_hadd_epi32(result3, result4); // L3L3L3L3 + L4L4L4L4 = L3L3 L4L4
      __m128i result = _mm_hadd_epi32(sumQuad12, sumQuad34); // L1L1 L2L2 + L3L3 L4L4 = L1 L2 L3 L4

      // correct if signed, scale back, store
      if (!lessthan16bit)
        result = _mm_add_epi32(result, shiftfromsigned);
      result = _mm_srai_epi32(result, FPScale16bits); // shift back integer arithmetic 13 bits precision

      __m128i result_4x_uint16 = _MM_PACKUS_EPI32(result, zero); // 4*32+zeros = lower 4*16 OK

      // extra clamp for 10-14 bit
        result_4x_uint16 = _MM_MIN_EPU16(result_4x_uint16,clamp_limit_max);
		result_4x_uint16 = _MM_MAX_EPU16(result_4x_uint16,clamp_limit_min);

      _mm_storel_epi64(reinterpret_cast<__m128i *>(dst + x), result_4x_uint16);

    }

    dst += dst_pitch;
    src += src_pitch;
  }
}


// filter_size <= 8 -> filter_size_align8 == 1 -> no loop, hope it'll be optimized
// filter_size <= 16 -> filter_size_align8 == 2 -> loop 0..1 hope it'll be optimized
// filter_size > 16 -> use parameter AlignNumber(program->filter_size_numOfFullBlk8, 8) / 8;
template<bool lessthan16bit, int filtersizealigned8>
#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
static void internal_resizer_h_sse41_generic_uint16_t(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // 1 and 2: special case for compiler optimization
  const int filter_size_numOfBlk8 = (filtersizealigned8 >= 1) ? filtersizealigned8 : (AlignNumber(program->filter_size,8) >> 3);

  const __m128i zero = _mm_setzero_si128();
  const __m128i shifttosigned = _mm_set1_epi16(-32768); // for 16 bits only
  const __m128i shiftfromsigned = _mm_set1_epi32(+32768 << FPScale16bits); // for 16 bits only
  const __m128i rounder = _mm_set_epi32(0, 0, 0, 1 << (FPScale16bits - 1)); // only once

  const uint16_t *src = reinterpret_cast<const uint16_t *>(src8);
  uint16_t *dst = reinterpret_cast<uint16_t *>(dst8);
  dst_pitch >>= 1;
  src_pitch >>= 1;

	const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
		((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));

  __m128i clamp_limit_min = _mm_set1_epi16(val_min);
  __m128i clamp_limit_max = _mm_set1_epi16(val_max);

  for (int y = 0; y < height; y++)
  {
    short* current_coeff = program->pixel_coefficient;

    for (int x = 0; x < width; x += 4)
	{
      __m128i result1 = rounder;
      __m128i result2 = result1;
      __m128i result3 = result1;
      __m128i result4 = result1;

      int begin1 = program->pixel_offset[x + 0];
      int begin2 = program->pixel_offset[x + 1];
      int begin3 = program->pixel_offset[x + 2];
      int begin4 = program->pixel_offset[x + 3];

      // this part is repeated 4 times
      // begin1, result1
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin1, i, current_coeff, result1, shifttosigned);

      // begin2, result2
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin2, i, current_coeff, result2, shifttosigned);

      // begin3, result3
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin3, i, current_coeff, result3, shifttosigned);

      // begin4, result4
      for (int i = 0; i < filter_size_numOfBlk8; i++)
        process_one_pixel_h_uint16_t<lessthan16bit>(src, begin4, i, current_coeff, result4, shifttosigned);

      const __m128i sumQuad12 = _mm_hadd_epi32(result1, result2); // L1L1L1L1 + L2L2L2L2 = L1L1 L2L2
      const __m128i sumQuad34 = _mm_hadd_epi32(result3, result4); // L3L3L3L3 + L4L4L4L4 = L3L3 L4L4
      __m128i result = _mm_hadd_epi32(sumQuad12, sumQuad34); // L1L1 L2L2 + L3L3 L4L4 = L1 L2 L3 L4

      // correct if signed, scale back, store
      if (!lessthan16bit)
        result = _mm_add_epi32(result, shiftfromsigned);
      result = _mm_srai_epi32(result, FPScale16bits); // shift back integer arithmetic 13 bits precision

      __m128i result_4x_uint16 = _mm_packus_epi32(result, zero); // 4*32+zeros = lower 4*16 OK

      // extra clamp for 10-14 bit
        result_4x_uint16 = _mm_min_epu16(result_4x_uint16,clamp_limit_max);
		result_4x_uint16 = _mm_max_epu16(result_4x_uint16,clamp_limit_min);

      _mm_storel_epi64(reinterpret_cast<__m128i *>(dst + x), result_4x_uint16);

    }

    dst += dst_pitch;
    src += src_pitch;
  }
}


//-------- 128 bit uint16_t Horizontal Dispatcher

template<bool lessthan16bit>
#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void resizer_h_ssse3_generic_uint16_t(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size_numOfBlk8 = AlignNumber(program->filter_size,8) >> 3;

  if (filter_size_numOfBlk8 == 1)
    internal_resizer_h_ssse3_generic_uint16_t<lessthan16bit, 1>(dst8,src8,dst_pitch,src_pitch,program,width,height,bits_per_pixel,range,mode_YUY2);
  else if (filter_size_numOfBlk8 == 2)
    internal_resizer_h_ssse3_generic_uint16_t<lessthan16bit, 2>(dst8,src8,dst_pitch,src_pitch,program,width,height,bits_per_pixel,range,mode_YUY2);
  else // -1: basic method, use program->filter_size
    internal_resizer_h_ssse3_generic_uint16_t<lessthan16bit, -1>(dst8,src8,dst_pitch,src_pitch,program,width,height,bits_per_pixel,range,mode_YUY2);
}


template<bool lessthan16bit>
#if defined(CLANG)
__attribute__((__target__("sse4.1")))
#endif
static void resizer_h_sse41_generic_uint16_t(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size_numOfBlk8 = AlignNumber(program->filter_size,8) >> 3;

  if (filter_size_numOfBlk8 == 1)
    internal_resizer_h_sse41_generic_uint16_t<lessthan16bit, 1>(dst8,src8,dst_pitch,src_pitch,program,width,height,bits_per_pixel,range,mode_YUY2);
  else if (filter_size_numOfBlk8 == 2)
    internal_resizer_h_sse41_generic_uint16_t<lessthan16bit, 2>(dst8,src8,dst_pitch,src_pitch,program,width,height,bits_per_pixel,range,mode_YUY2);
  else // -1: basic method, use program->filter_size
    internal_resizer_h_sse41_generic_uint16_t<lessthan16bit, -1>(dst8,src8,dst_pitch,src_pitch,program,width,height,bits_per_pixel,range,mode_YUY2);
}


#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void resizer_h_ssse3_generic(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = AlignNumber(program->filter_size,8) >> 3;
  const __m128i zero = _mm_setzero_si128();

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int Offset = 1 << (FPScale8bits-1);

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = 0; y < height; y++)
  {
    const short *current_coeff = program->pixel_coefficient;
	
    for (int x = 0; x < width; x+=4)
	{
      __m128i result1 = _mm_setr_epi32(Offset, 0, 0, 0);
      __m128i result2 = _mm_setr_epi32(Offset, 0, 0, 0);
      __m128i result3 = _mm_setr_epi32(Offset, 0, 0, 0);
      __m128i result4 = _mm_setr_epi32(Offset, 0, 0, 0);

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

	  for (int i = 0; i < filter_size; i++)
	  {
	    __m128i data, coeff, current_result;
		data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin1+(i << 3))); // 8 * 8 bit pixels
        data = _mm_unpacklo_epi8(data, zero);
	    coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff)); // 8 coeffs 14 bit scaled -> ushort OK
		current_result = _mm_madd_epi16(data, coeff);
        result1 = _mm_add_epi32(result1, current_result);
			
	    current_coeff += 8;		
	  }

      for (int i = 0; i < filter_size; i++)
	  {
		__m128i data, coeff, current_result;
        data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin2+(i << 3)));
	    data = _mm_unpacklo_epi8(data, zero);
		coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
        current_result = _mm_madd_epi16(data, coeff);
	    result2 = _mm_add_epi32(result2, current_result);

		current_coeff += 8;
      }

      for (int i = 0; i < filter_size; i++)
	  {
		__m128i data, coeff, current_result;
        data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin3+(i << 3)));
	    data = _mm_unpacklo_epi8(data, zero);
		coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
        current_result = _mm_madd_epi16(data, coeff);
	    result3 = _mm_add_epi32(result3, current_result);

		current_coeff += 8;
	  }

      for (int i = 0; i < filter_size; i++)
	  {
		__m128i data, coeff, current_result;
        data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin4+(i << 3)));
		data = _mm_unpacklo_epi8(data, zero);
	    coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
        current_result = _mm_madd_epi16(data, coeff);
	    result4 = _mm_add_epi32(result4, current_result);

        current_coeff += 8;
	  }

      __m128i result12 = _mm_hadd_epi32(result1, result2);
      __m128i result34 = _mm_hadd_epi32(result3, result4);
      __m128i result = _mm_hadd_epi32(result12, result34);

      result = _mm_srai_epi32(result, FPScale8bits);

      result = _mm_packs_epi32(result, zero);
      result = _mm_packus_epi16(result, zero);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      *((int*)(dst+x)) = _mm_cvtsi128_si32(result);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}


#if defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
static void resizer_h_ssse3_8(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const __m128i zero = _mm_setzero_si128();

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int Offset = 1 << (FPScale8bits-1);

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = 0; y < height; y++)
  {
    short *current_coeff = program->pixel_coefficient;
	
    for (int x = 0; x < width; x+=4)
	{
      __m128i result1 = _mm_setr_epi32(Offset, 0, 0, 0);
      __m128i result2 = _mm_setr_epi32(Offset, 0, 0, 0);
      __m128i result3 = _mm_setr_epi32(Offset, 0, 0, 0);
      __m128i result4 = _mm_setr_epi32(Offset, 0, 0, 0);

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      __m128i data, coeff, current_result;

      // Unroll 1
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin1));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result1 = _mm_add_epi32(result1, current_result);

      current_coeff += 8;

      // Unroll 2
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin2));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result2 = _mm_add_epi32(result2, current_result);

      current_coeff += 8;

      // Unroll 3
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin3));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result3 = _mm_add_epi32(result3, current_result);

      current_coeff += 8;

      // Unroll 4
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin4));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result4 = _mm_add_epi32(result4, current_result);

      current_coeff += 8;

      // Combine
      __m128i result12 = _mm_hadd_epi32(result1, result2);
      __m128i result34 = _mm_hadd_epi32(result3, result4);
      __m128i result = _mm_hadd_epi32(result12, result34);

      result = _mm_srai_epi32(result, FPScale8bits);

      result = _mm_packs_epi32(result, zero);
      result = _mm_packus_epi16(result, zero);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      *((int*)(dst+x)) = _mm_cvtsi128_si32(result);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}


/********************************************************************
***** Declare index of new filters for Avisynth's filter engine *****
********************************************************************/



FilteredResizeH::FilteredResizeH( PClip _child, double subrange_left, double subrange_width,int target_width,
	uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,bool _avsp, ResamplingFunction* func, IScriptEnvironment* env )
  : GenericVideoFilter(_child),
  resampling_program_luma(NULL), resampling_program_chroma(NULL),
  filter_storage_luma(NULL), filter_storage_chroma(NULL),threads(_threads),sleep(_sleep),
  avsp(_avsp)
{
  int filter_sz;
  
  src_width  = vi.width;
  src_height = vi.height;
  dst_width  = target_width;
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
  Enable_AVX2 = avsp && ((env->GetCPUFlags() & CPUF_AVX2)!=0);

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

	const int src_width = vi.IsPlanar() ? vi.width : vi.BytesFromPixels(vi.width)/pixelsize;
	const int dst_width = vi.IsPlanar() ? target_width : vi.BytesFromPixels(target_width)/pixelsize;

	if (vi.height<32) threads_number=1;
	else threads_number=threads;

  // Main resampling program
  int SizeH;

  if (desample) resampling_program_luma = func->GetDesamplingProgram(target_width, subrange_left, subrange_width, vi.width, bits_per_pixel, accuracy, 0, shift_w, SizeH, env);
  else
  {
	  resampling_program_luma = func->GetResamplingProgram(vi.width, subrange_left, subrange_width, target_width, bits_per_pixel, env);
	  SizeH=dst_width;
  }

  if (resampling_program_luma==NULL)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  if (desample) env->ThrowError("ResizeHMT: Error while GetDesamplingProgram for luma!");
	  else env->ThrowError("ResizeHMT: Error while GetResamplingProgram for luma!");
  }

  filter_sz=resampling_program_luma->filter_size;
  if (vi.width<filter_sz)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  env->ThrowError("ResizeHMT: Source luma width (%d) is too small for this resizing method, must be minimum of %d!",vi.width,filter_sz);
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

	filter_sz=resampling_program_chroma->filter_size;
	if (w_UV<filter_sz)
	{
		FreeData();
		if (threads>1) poolInterface->DeAllocateAllThreads(true);
		env->ThrowError("ResizeHMT: Source chroma width (%d) is too small for this resizing method, must be minimum of %d!",w_UV,filter_sz);
	}
	
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
	
  int8_t nPool=-1;

  if (threads_number>1)
  {
	if ((!poolInterface->RequestThreadPool(UserId,threads_number,MT_ThreadGF,nPool,false,true)) || (nPool==-1))
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
		if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);

		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=2;
			if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=3;
			if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);
		}
		else
		{
			if (isRGBPfamily)
			{
				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=4;
				if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);

				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=5;
				if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);								
			}
		}

		if (isAlphaChannel)
		{
			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=6;
			if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);												
		}

		for(uint8_t i=0; i<threads_number; i++)
			MT_ThreadGF[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep,nPool);
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
	if (pixelsize==1)
	{
		if (Enable_SSSE3)
		{
#ifdef AVX2_BUILD_POSSIBLE				
			if (Enable_AVX2)
			{
				// make the resampling coefficient array mod16 friendly for simd, padding non-used coeffs with zeros
				resize_h_prepare_coeff_8or16(program,env,16);
				return resizer_h_avx2_generic_uint8_t;
			}
			else
#endif			
			{
				// make the resampling coefficient array mod8 friendly for simd, padding non-used coeffs with zeros
				resize_h_prepare_coeff_8or16(program,env,8);
				if (program->filter_size>8) return resizer_h_ssse3_generic;
				else return resizer_h_ssse3_8;
			}
		}
		else return resize_h_c_planar;
	}
	else if (pixelsize==2)
	{ 
		if (Enable_SSSE3)
		{
			resize_h_prepare_coeff_8or16(program,env,8); // alignment of 8 is enough for AVX2 uint16_t as well
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
			resize_h_prepare_coeff_8or16(program,env,ALIGN_FLOAT_RESIZER_COEFF_SIZE); // alignment of 8 is enough for AVX2 float as well

			const int filtersizealign8 = AlignNumber(program->filter_size,8);
			const int filtersizemod8 = program->filter_size & 7;

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
	uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC,
	bool _avsp,ResamplingFunction* func, IScriptEnvironment* env )
  : GenericVideoFilter(_child),
    resampling_program_luma(NULL), resampling_program_chroma(NULL),
    src_pitch_table_luma(NULL), src_pitch_table_chromaU(NULL), src_pitch_table_chromaV(NULL),
    src_pitch_luma(-1), src_pitch_chromaU(-1), src_pitch_chromaV(-1),
    filter_storage_luma_aligned(NULL), filter_storage_luma_unaligned(NULL),
    filter_storage_chroma_aligned(NULL), filter_storage_chroma_unaligned(NULL),
	sleep(_sleep),threads(_threads),avsp(_avsp)
{
	int16_t i,filter_sz;

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
  Enable_AVX2 = avsp && ((env->GetCPUFlags() & CPUF_AVX2)!=0);

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

  if (vi.IsRGB() && !isRGBPfamily)
  {
	  if (desample)
		subrange_top = target_height - subrange_top - subrange_height; // why?
	  else
		subrange_top = vi.height - subrange_top - subrange_height; // why?
  }

  // Create resampling program and pitch table
  int SizeV;

  if (ShiftC==0) ShiftC=shift_h;

  if (desample) resampling_program_luma  = func->GetDesamplingProgram(target_height, subrange_top, subrange_height, vi.height, bits_per_pixel, accuracy, ChromaS, ShiftC, SizeV, env);
  else
  {
	  resampling_program_luma  = func->GetResamplingProgram(vi.height, subrange_top, subrange_height, target_height, bits_per_pixel, env);
	  SizeV=target_height;
  }

  if (resampling_program_luma==NULL)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  if (desample) env->ThrowError("ResizeVMT: Error while GetDesamplingProgram for luma!");
	  else env->ThrowError("ResizeVMT: Error while GetResamplingProgram for luma!");
  }

  filter_sz=resampling_program_luma->filter_size;
  if (vi.height<filter_sz)
  {
	  FreeData();
	  if (threads>1) poolInterface->DeAllocateAllThreads(true);
	  env->ThrowError("ResizeVMT: Source luma height (%d) is too small for this resizing method, must be minimum of %d!",vi.height,filter_sz);
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

    filter_sz=resampling_program_chroma->filter_size;
	if (h_UV<filter_sz)
	{
		FreeData();
		if (threads>1) poolInterface->DeAllocateAllThreads(true);
		env->ThrowError("ResizeVMT: Source chroma height (%d) is too small for this resizing method, must be minimum of %d!",h_UV,filter_sz);
	}
	
    resampler_chroma_aligned = GetResampler(true,filter_storage_chroma_aligned,resampling_program_chroma);
    resampler_chroma_unaligned = GetResampler(false,filter_storage_chroma_unaligned,resampling_program_chroma);
  }

  resampler_luma_aligned   = GetResampler(true,filter_storage_luma_aligned,resampling_program_luma);
  resampler_luma_unaligned = GetResampler(false,filter_storage_luma_unaligned,resampling_program_luma);

  threads_number=CreateMTData(threads_number,work_width,vi.height,work_width,SizeV,shift_w,shift_h);

	if (threads_number>1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			FreeData();
			poolInterface->DeAllocateAllThreads(true);
			env->ThrowError("ResizeVMT: Error with the TheadPool while getting UserId!");
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

  int8_t nPool=-1;

  if (threads_number>1)
  {
	if ((!poolInterface->RequestThreadPool(UserId,threads_number,MT_ThreadGF,nPool,false,true)) || (nPool==-1))
		env->ThrowError("ResizeHMT: Error with the TheadPool while requesting threadpool!");
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
		if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);

		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0)) f_proc=3;
			else f_proc=4;

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);

			if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0)) f_proc=5;
			else f_proc=6;

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);
		}
		else
		{
			if (isRGBPfamily)
			{
				if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0)) f_proc=7;
				else f_proc=8;

				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=f_proc;
				if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);

				if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0)) f_proc=9;
				else f_proc=10;

				for(uint8_t i=0; i<threads_number; i++)
					MT_ThreadGF[i].f_process=f_proc;
				if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);							
			}
		}
		
		if (isAlphaChannel)
		{
			if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0)) f_proc=11;
			else f_proc=12;

			for(uint8_t i=0; i<threads_number; i++)
				MT_ThreadGF[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId,nPool)) poolInterface->WaitThreadsEnd(UserId,nPool);			
		}

		for(uint8_t i=0; i<threads_number; i++)
			MT_ThreadGF[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep,nPool);
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


ResamplerV FilteredResizeV::GetResampler(bool aligned,void*& storage, ResamplingProgram* program)
{
  if (program->filter_size==1)
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


PClip FilteredResizeMT::CreateResizeV(PClip clip, double subrange_top, double subrange_height, int target_height,
                    uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC,
					bool _avsp,ResamplingFunction* func,IScriptEnvironment* env)
{
  return new FilteredResizeV(clip, subrange_top, subrange_height, target_height,_threads,_sleep,range_mode,desample,
	  accuracy,ChromaS,ShiftC,_avsp,func, env);
}


PClip FilteredResizeMT::CreateResizeH(PClip clip, double subrange_top, double subrange_height, int target_height,
                    uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,bool _avsp,ResamplingFunction* func,
					IScriptEnvironment* env)
{
  return new FilteredResizeH(clip, subrange_top, subrange_height, target_height,_threads,_sleep,range_mode,desample,
	  accuracy,_avsp,func, env);
}


PClip FilteredResizeMT::CreateResize(PClip clip, int target_width, int target_height,int _threads,
	bool _LogicalCores,bool _MaxPhysCores, bool _SetAffinity,bool _sleep,int prefetch,
	int range_mode,bool desample,int accuracy,int order,int thread_level,
	const AVSValue* args,ResamplingFunction* f,IScriptEnvironment* env)
{
  const VideoInfo& vi = clip->GetVideoInfo();
  const double subrange_left = args[0].AsFloat(0), subrange_top = args[1].AsFloat(0);
  
  if (target_height <= 0)
    env->ThrowError("ResizeMT: Height must be greater than 0.");

  if (target_width <= 0) {
    env->ThrowError("ResizeMT: Width must be greater than 0.");
  }

  if ((range_mode<0) || (range_mode>4)) env->ThrowError("ResizeMT: [range] must be between 0 and 4.");

  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("ResizeMT: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
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

  const int shift = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
  int SizeH;

  if (desample) SizeH=f->GetDesamplingData(target_width, subrange_left, subrange_width, vi.width,bits_per_pixel,shift,env);
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
		if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)) result=clip;
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
				else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
			}
			else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
		}
		if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)))
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,int(subrange_left),0,int(subrange_width),vi.height,true};
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
							v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[0],desample,accuracy,0,shift,avsp, f, env);

							VideoInfo vR = v.AsClip()->GetVideoInfo();
							int ChromaS=vR.height >> shift;

							vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
							vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
							if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[3],desample,accuracy,0,shift,avsp, f, env);
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
							result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
							result=env->Invoke(turnLeftFunction,result).AsClip();
						}
					}
					else
					{
						result=env->Invoke(turnLeftFunction,result).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
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
						v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[0],desample,accuracy,0,shift,avsp, f, env);

						VideoInfo vR = v.AsClip()->GetVideoInfo();
						int ChromaS=vR.height >> shift;

						vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
						vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
						if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[3],desample,accuracy,0,shift,avsp, f, env);
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
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
						result=env->Invoke(turnLeftFunction,result).AsClip();
					}
				}
				else
				{
					result=env->Invoke(turnLeftFunction,result).AsClip();
					result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
					result=env->Invoke(turnRightFunction,result).AsClip();
				}
			}
		}
	  }
	  else
	  {
		if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)) result=clip;
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
							v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[0],desample,accuracy,0,shift,avsp, f, env);

							VideoInfo vR = v.AsClip()->GetVideoInfo();
							int ChromaS=vR.height >> shift;

							vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
							vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
							if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[3],desample,accuracy,0,shift,avsp, f, env);
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
							result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
							result=env->Invoke(turnLeftFunction,result).AsClip();
						}
					}
					else
					{
						result=env->Invoke(turnLeftFunction,clip).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
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
						v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[0],desample,accuracy,0,shift,avsp, f, env);

						VideoInfo vR = v.AsClip()->GetVideoInfo();
						int ChromaS=vR.height >> shift;

						vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
						vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
						if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[3],desample,accuracy,0,shift,avsp, f, env);
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
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
						result=env->Invoke(turnLeftFunction,result).AsClip();
					}
				}
				else
				{
					result=env->Invoke(turnLeftFunction,clip).AsClip();
					result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
					result=env->Invoke(turnRightFunction,result).AsClip();
				}
			}
		}
		if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)))
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,0,int(subrange_top),vi.width,int(subrange_height),true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
			}
			else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
		}
	  }
  }
  else
  {	  
	  if (VFirst)
	  {
		if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)) result=clip;
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
				else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0, avsp, f, env);
			}
			else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
		}
		if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)))
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,int(subrange_left),0,int(subrange_width),vi.height,true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeH(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,avsp, f, env);
			}
			else result = CreateResizeH(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,avsp, f, env);
		}		
	  }
	  else
	  {
		if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)) result=clip;
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
				else result = CreateResizeH(clip, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,avsp, f, env);
			}
			else result = CreateResizeH(clip, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,avsp, f, env);
		}
		if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)))
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,0,int(subrange_top),vi.width,int(subrange_height),true};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
			}
			else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
		}
	  }
  }
  
  return result;
}


AVSValue __cdecl FilteredResizeMT::Create_PointResize(AVSValue args, void*, IScriptEnvironment* env)
{
  PointFilter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,args[14].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_BilinearResize(AVSValue args, void*, IScriptEnvironment* env)
{
  TriangleFilter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,args[14].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_BicubicResize(AVSValue args, void*, IScriptEnvironment* env)
{
  MitchellNetravaliFilter f(args[3].AsDblDef(1./3.), args[4].AsDblDef(1./3.));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[9].AsInt(0),
	  args[10].AsBool(true),args[11].AsBool(true),args[12].AsBool(false),args[13].AsBool(false),
	  args[14].AsInt(0),args[15].AsInt(1),false,0,0,args[16].AsInt(6),&args[5],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_LanczosResize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(args[7].AsInt(3));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,args[15].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Lanczos4Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(4);
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,args[14].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_BlackmanResize(AVSValue args, void*, IScriptEnvironment* env)
{
  BlackmanFilter f(args[7].AsInt(4));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,args[15].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Spline16Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline16Filter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,args[14].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Spline36Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline36Filter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,args[14].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_Spline64Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline64Filter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,args[14].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_GaussianResize(AVSValue args, void*, IScriptEnvironment* env)
{
  GaussianFilter f(args[7].AsFloat(30.0f));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,args[15].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_SincResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincFilter f(args[7].AsInt(4));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,args[15].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_SinPowerResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SinPowerFilter f(args[7].AsFloat(2.5f));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,args[15].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_SincLin2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincLin2Filter f(args[7].AsInt(15));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,args[15].AsInt(6),&args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_UserDefined2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
	UserDefined2Filter f(args[3].AsFloat(121.0),args[4].AsFloat(19.0));
	return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[9].AsInt(0),
		args[10].AsBool(true),args[11].AsBool(true),args[12].AsBool(false),args[13].AsBool(false),
		args[14].AsInt(0),args[15].AsInt(1),false,0,0,args[16].AsInt(6),&args[5],&f,env);
}

// Desample functions

AVSValue __cdecl FilteredResizeMT::Create_DeBilinearResize(AVSValue args, void*, IScriptEnvironment* env)
{
  TriangleFilter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),args[16].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeBicubicResize(AVSValue args, void*, IScriptEnvironment* env)
{
  MitchellNetravaliFilter f(args[3].AsDblDef(1./3.), args[4].AsDblDef(1./3.));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[9].AsInt(0),
	  args[10].AsBool(true),args[11].AsBool(true),args[12].AsBool(false),args[13].AsBool(false),
	  args[14].AsInt(0),args[15].AsInt(1),true,args[16].AsInt(0),args[17].AsInt(0),args[18].AsInt(6),
	  &args[5],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeLanczosResize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(args[7].AsInt(3));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),args[17].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeLanczos4Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(4);
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),args[16].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeBlackmanResize(AVSValue args, void*, IScriptEnvironment* env)
{
  BlackmanFilter f(args[7].AsInt(4));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),args[17].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline16Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline16Filter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),args[16].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline36Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline36Filter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),args[16].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline64Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline64Filter f;
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),args[16].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeGaussianResize(AVSValue args, void*, IScriptEnvironment* env)
{
  GaussianFilter f(args[7].AsFloat(30.0f));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),args[17].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSincResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincFilter f(args[7].AsInt(4));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),args[17].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSinPowerResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SinPowerFilter f(args[7].AsFloat(2.5f));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),args[17].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeSincLin2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincLin2Filter f(args[7].AsInt(15));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),args[17].AsInt(6),
	  &args[3],&f,env);
}

AVSValue __cdecl FilteredResizeMT::Create_DeUserDefined2Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  UserDefined2Filter f(args[3].AsFloat(121.0f),args[4].AsFloat(19.0f));
  return CreateResize(args[0].AsClip(),args[1].AsInt(),args[2].AsInt(),args[9].AsInt(0),
	  args[10].AsBool(true),args[11].AsBool(true),args[12].AsBool(false),args[13].AsBool(false),
	  args[14].AsInt(0),args[15].AsInt(1),true,args[16].AsInt(0),args[17].AsInt(0),args[18].AsInt(6),
	  &args[5],&f,env);
}
