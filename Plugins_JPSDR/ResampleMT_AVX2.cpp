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

#if _MSC_VER >= 1900

#include "./ResampleMT_AVX2.h"

extern "C" void Resize_V_AVX2_Planar_8bits_ASM(const void *src,void *dst,const void *coeff,uint32_t width16,ptrdiff_t src_pitch,
	uint32_t kernel_size_2,const uint32_t *valmin,const uint32_t *valmax,const uint32_t *rounder);
extern "C" void Resize_V_AVX2_Planar_10to14bits_ASM(const void *src,void *dst,const void *coeff,uint32_t width16,ptrdiff_t src_pitch,
	uint32_t kernel_size_2,const uint32_t *valmin,const uint32_t *valmax,const uint32_t *rounder);
extern "C" void Resize_V_AVX2_Planar_16bits_ASM(const void *src,void *dst,const void *coeff,uint32_t width16, ptrdiff_t src_pitch,
	uint32_t kernel_size_2,const uint32_t *valmin,const uint32_t *valmax,const uint32_t *rounder,
	const uint32_t *shifttosigned,const uint32_t *shiftfromsigned);
extern "C" void Resize_V_AVX2_Planar_32bits_ASM(const void *src,void *dst,const void *coeff,uint32_t width8,ptrdiff_t src_pitch,
	uint32_t kernel_size_2);
	
////////

extern "C" void Resize_H_AVX2_Planar_8bits_ASM(const void *src,void *dst,const void *coeff,ptrdiff_t src_pitch,ptrdiff_t dst_pitch,
	uint32_t kernel_size_32,uint32_t sizeh,const uint32_t *valmin,const uint32_t *valmax,const uint32_t *rounder);
extern "C" void Resize_H_AVX2_Planar_10to14bits_ASM(const void *src,void *dst,const void *coeff,ptrdiff_t src_pitch,ptrdiff_t dst_pitch,
	uint32_t kernel_size_32,uint32_t sizeh,const uint32_t *valmin,const uint32_t *valmax,const uint32_t *rounder);
extern "C" void Resize_H_AVX2_Planar_16bits_ASM(const void *src,void *dst,const void *coeff,ptrdiff_t src_pitch,ptrdiff_t dst_pitch,
	uint32_t kernel_size_32,uint32_t sizeh,const uint32_t *valmin,const uint32_t *valmax,const uint32_t *rounder,
	const uint32_t *shifttosigned,const uint32_t *shiftfromsigned);
extern "C" void Resize_H_AVX2_Planar_32bits_ASM(const void *src,void *dst,const void *coeff,ptrdiff_t src_pitch,ptrdiff_t dst_pitch,
	uint32_t kernel_size_32,uint32_t sizeh);


// Vertical Resizer

void Resize_V_AVX2_8bits(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage, const uint8_t range, const bool mode_YUY2)
{
	const ptrdiff_t src_pitch_asm = (ptrdiff_t)src_pitch;
	const uint32_t width16 = (width+15)>>4;
	
	const int filter_size = program->filter_size;
	const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
	const uint32_t kernel_size_2 = (program->filter_size_real+1)>>1;

	const uint32_t val_min0 = (range==1) ? 0 : 16;
	uint32_t val_max0 = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	val_max0 = (mode_YUY2 && ((range>=2) && (range<=3))) ? (((uint32_t)240 << 8)|235) : (val_max0 | (val_max0 << 8));
	
	const uint32_t val_min = val_min0 | (val_min0 << 8) | (val_min0 << 16) | (val_min0 << 24);
	const uint32_t val_max = val_max0 | (val_max0 << 16);
	const uint32_t rounder = 1 << (FPScale8bits - 1);

	for (int y = MinY; y < MaxY; y++)
	{
		const BYTE *src_ptr = src8 + pitch_table[program->pixel_offset[y]];
		
		Resize_V_AVX2_Planar_8bits_ASM(src_ptr,dst8,current_coeff,width16,src_pitch_asm,
			kernel_size_2,&val_min,&val_max,&rounder);
			
		dst8 += dst_pitch;
		current_coeff += filter_size;
	}
}


void Resize_V_AVX2_16bits(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage, const uint8_t range, const bool mode_YUY2)
{
	const ptrdiff_t src_pitch_asm = (ptrdiff_t)src_pitch;
	const uint32_t width16 = (width+15)>>4;

	const uint16_t *src = (const uint16_t*)src8;
	uint16_t *dst = (uint16_t*)dst8;
	dst_pitch = dst_pitch / sizeof(uint16_t);

	const uint32_t val_min0 = (range==1) ? 0 : (uint32_t)16 << (bits_per_pixel-8);
	const uint32_t val_max0 = ((range==1) || (range==4)) ? ((uint32_t)1 << bits_per_pixel)-1 : (range==2) ?
		((uint32_t)235 << (bits_per_pixel-8)) : ((uint32_t)240 << (bits_per_pixel-8));
	const uint32_t val_min = val_min0 | (val_min0 << 16);
	const uint32_t val_max = val_max0 | (val_max0 << 16);
	const uint32_t rounder = (uint32_t)1 << (FPScale16bits - 1);
	const uint32_t shifttosigned = ((uint32_t)32768) | ((uint32_t)32768 << 16);
	const uint32_t shiftfromsigned = 32768 << FPScale16bits;
	
	const int filter_size = program->filter_size;
	const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
	const uint32_t kernel_size_2 = (program->filter_size_real+1)>>1;

	if (bits_per_pixel==16)
	{
		for (int y = MinY; y < MaxY; y++)
		{
			const uint16_t *src_ptr = src + pitch_table[program->pixel_offset[y]];
			
			Resize_V_AVX2_Planar_16bits_ASM(src_ptr,dst,current_coeff,width16,src_pitch_asm,
				kernel_size_2,&val_min,&val_max,&rounder,&shifttosigned,&shiftfromsigned);
			
			dst += dst_pitch;
			current_coeff += filter_size;
		}
	}
	else
	{
		for (int y = MinY; y < MaxY; y++)
		{
			const uint16_t *src_ptr = src + pitch_table[program->pixel_offset[y]];

			Resize_V_AVX2_Planar_10to14bits_ASM(src_ptr,dst,current_coeff,width16,src_pitch_asm,
				kernel_size_2,&val_min,&val_max,&rounder);

			dst += dst_pitch;
			current_coeff += filter_size;
		}
	}
}


void Resize_V_AVX2_32bits(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage, const uint8_t range, const bool mode_YUY2)
{
	const ptrdiff_t src_pitch_asm = (ptrdiff_t)src_pitch;
	const uint32_t width8 = (width+7)>>3;
	
	const int filter_size = program->filter_size;
	const float *current_coeff = program->pixel_coefficient_float + filter_size*MinY;
	const uint32_t kernel_size_2 = (program->filter_size_real+1)>>1;

	const float *src = (const float*)src8;
	float *dst = (float*)dst8;
	dst_pitch = dst_pitch / sizeof(float);

	for (int y = MinY; y < MaxY; y++)
	{
		const float *src_ptr = src + pitch_table[program->pixel_offset[y]];

		Resize_V_AVX2_Planar_32bits_ASM(src_ptr,dst,current_coeff,width8,src_pitch_asm,kernel_size_2);
		
		dst += dst_pitch;
		current_coeff += filter_size;
	}
}


// Horizontal Resizer

void Resize_H_AVX2_8bits(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height,
	int bits_per_pixel, const uint8_t range, const bool mode_YUY2)
{
	const ptrdiff_t src_pitch_asm = (ptrdiff_t)src_pitch;
	const ptrdiff_t dst_pitch_asm = (ptrdiff_t)dst_pitch;
	const uint32_t height_asm = (uint32_t)height;

	const int filter_size = program->filter_size;
	const short *current_coeff = program->pixel_coefficient;
	const uint32_t kernel_size_32 = (program->filter_size_real + 15) >> 4;

	const uint32_t val_min0 = (range == 1) ? 0 : 16;
	uint32_t val_max0 = ((range == 1) || (range == 4)) ? 255 : (range == 2) ? 235 : 240;
	val_max0 = (mode_YUY2 && ((range >= 2) && (range <= 3))) ? (((uint32_t)240 << 8) | 235) : (val_max0 | (val_max0 << 8));

	const uint32_t val_min = val_min0 | (val_min0 << 8) | (val_min0 << 16) | (val_min0 << 24);
	const uint32_t val_max = val_max0 | (val_max0 << 16);
	const uint32_t rounder = 1 << (FPScale8bits - 1);

	for (int x = 0; x < width; x++)
	{
		const BYTE *src_ptr = src8 + program->pixel_offset[x];

		Resize_H_AVX2_Planar_8bits_ASM(src_ptr,dst8,current_coeff,src_pitch_asm,dst_pitch_asm,
			kernel_size_32,height_asm,&val_min,&val_max,&rounder);

		dst8++;
		current_coeff += filter_size;
	}
}


void Resize_H_AVX2_16bits(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height,
	int bits_per_pixel, const uint8_t range, const bool mode_YUY2)
{
	const ptrdiff_t src_pitch_asm = (ptrdiff_t)src_pitch;
	const ptrdiff_t dst_pitch_asm = (ptrdiff_t)dst_pitch;
	const uint32_t height_asm = (uint32_t)height;

	const int filter_size = program->filter_size;
	const short *current_coeff = program->pixel_coefficient;
	const uint32_t kernel_size_32 = (program->filter_size_real + 15) >> 4;

	const uint16_t *src = (const uint16_t*)src8;
	uint16_t *dst = (uint16_t*)dst8;

	const uint32_t val_min0 = (range==1) ? 0 : (uint32_t)16 << (bits_per_pixel-8);
	const uint32_t val_max0 = ((range==1) || (range==4)) ? ((uint32_t)1 << bits_per_pixel)-1 : (range==2) ?
		((uint32_t)235 << (bits_per_pixel-8)) : ((uint32_t)240 << (bits_per_pixel-8));
	const uint32_t val_min = val_min0 | (val_min0 << 16);
	const uint32_t val_max = val_max0 | (val_max0 << 16);
	const uint32_t rounder = (uint32_t)1 << (FPScale16bits - 1);
	const uint32_t shifttosigned = ((uint32_t)32768) | ((uint32_t)32768 << 16);
	const uint32_t shiftfromsigned = 32768 << FPScale16bits;

	if (bits_per_pixel==16)
	{
		for (int x = 0; x < width; x++)
		{
			const uint16_t *src_ptr = src + program->pixel_offset[x];

			Resize_H_AVX2_Planar_16bits_ASM(src_ptr,dst,current_coeff,src_pitch_asm,dst_pitch_asm,
				kernel_size_32,height_asm,&val_min,&val_max,&rounder,&shifttosigned,&shiftfromsigned);

			dst++;
			current_coeff += filter_size;
		}
	}
	else
	{
		for (int x = 0; x < width; x++)
		{
			const uint16_t *src_ptr = src + program->pixel_offset[x];

			Resize_H_AVX2_Planar_10to14bits_ASM(src_ptr,dst,current_coeff,src_pitch_asm,dst_pitch_asm,
				kernel_size_32,height_asm,&val_min,&val_max,&rounder);

			dst++;
			current_coeff += filter_size;
		}
	}
}


void Resize_H_AVX2_32bits(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height,
	int bits_per_pixel, const uint8_t range, const bool mode_YUY2)
{
	const ptrdiff_t src_pitch_asm = (ptrdiff_t)src_pitch;
	const ptrdiff_t dst_pitch_asm = (ptrdiff_t)dst_pitch;
	const uint32_t height_asm = (uint32_t)height;

	const int filter_size = program->filter_size;
	const float *current_coeff = program->pixel_coefficient_float;
	const uint32_t kernel_size_32 = (program->filter_size_real + 7) >> 3;

	const float *src = (const float*)src8;
	float *dst = (float*)dst8;

	for (int x = 0; x < width; x++)
	{
		const float *src_ptr = src + program->pixel_offset[x];

		Resize_H_AVX2_Planar_32bits_ASM(src_ptr,dst,current_coeff,src_pitch_asm,dst_pitch_asm,kernel_size_32,height_asm);

		dst++;
		current_coeff += filter_size;
	}
}

#endif
