/*
 *  AutoYUY2()
 *
 *  Adaptive YV12 upsampling. Progressive picture areas are upsampled
 *  progressively and interlaced areas are upsampled interlaced.
 *  Copyright (C) 2005 Donald A. Graft
 *	
 *  AutoYUY2 is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *   
 *  AutoYUY2 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *   
 *  You should have received a copy of the GNU General Public License
 *  along with GNU Make; see the file COPYING.  If not, write to
 *  the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA. 
 *
 */

#include <windows.h>
#include "AutoYUY2.h"

extern ThreadPoolInterface *poolInterface;

extern "C" void JPSDR_AutoYUY2_1(const uint8_t *scr_y,const uint8_t *src_u,const uint8_t *src_v,
		uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_1(const uint8_t *scr_y,const uint8_t *src_u,const uint8_t *src_v,
		uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_2(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_3(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_4(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);

extern "C" void JPSDR_AutoYUY2_SSE2_1b(const uint8_t *scr_y,const uint8_t *src_u,const uint8_t *src_v,
		uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_2b(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_3b(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_SSE2_4b(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);

extern "C" void JPSDR_AutoYUY2_AVX_1(const uint8_t *scr_y,const uint8_t *src_u,const uint8_t *src_v,
		uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_AVX_2(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_AVX_3(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_AVX_4(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);

extern "C" void JPSDR_AutoYUY2_AVX_1b(const uint8_t *scr_y,const uint8_t *src_u,const uint8_t *src_v,
		uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_AVX_2b(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_AVX_3b(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);
extern "C" void JPSDR_AutoYUY2_AVX_4b(const uint8_t *scr_y,const uint8_t *src_u1,const uint8_t *src_u2,
		const uint8_t *src_v1,const uint8_t *src_v2,uint8_t *dst,int32_t w);

extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(const void *scr_1,const void *src_2,void *dst,int w);

extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(const void *scr_1,const void *src_2,void *dst,int w);
extern "C" void JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(const void *scr_1,const void *src_2,void *dst,int w);

int __stdcall AutoYUY2::SetCacheHints(int cachehints,int frame_range)
{
  switch (cachehints)
  {
  case CACHE_GET_MTMODE :
    return MT_MULTI_INSTANCE;
  default :
    return 0;
  }
}

uint8_t AutoYUY2::CreateMTData(uint8_t max_threads,int32_t size_x,int32_t size_y)
{
	if ((max_threads<=1) || (max_threads>threads_number))
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=size_y;
		MT_Data[0].dst_Y_h_max=size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		MT_Data[0].src_UV_h_max=size_y >> 1;
		MT_Data[0].dst_UV_h_max=size_y;
		MT_Data[0].src_Y_w=size_x;
		if (output==0) MT_Data[0].dst_Y_w=size_x << 1;
		else MT_Data[0].dst_Y_w=size_x;
		MT_Data[0].src_UV_w=size_x >> 1;
		MT_Data[0].dst_UV_w=size_x >> 1;
		return(1);
	}

	int32_t dh_Y,dh_UV,h_y;
	uint8_t i,max=0;

	dh_Y=(size_y+(int32_t)max_threads-1)/(int32_t)max_threads;
	if (dh_Y<16) dh_Y=16;
	if ((dh_Y & 3)!=0) dh_Y=((dh_Y+3) >> 2) << 2;

	h_y=0;
	while (h_y<(size_y-16))
	{
		max++;
		h_y+=dh_Y;
	}

	if (max==1)
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=size_y;
		MT_Data[0].dst_Y_h_max=size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		MT_Data[0].src_UV_h_max=size_y >> 1;
		MT_Data[0].dst_UV_h_max=size_y;
		MT_Data[0].src_Y_w=size_x;
		if (output==0) MT_Data[0].dst_Y_w=size_x << 1;
		else MT_Data[0].dst_Y_w=size_x;
		MT_Data[0].src_UV_w=size_x >> 1;
		MT_Data[0].dst_UV_w=size_x >> 1;
		return(1);
	}

	dh_UV=dh_Y>>1; 

	MT_Data[0].top=true;
	MT_Data[0].bottom=false;
	MT_Data[0].src_Y_h_min=0;
	MT_Data[0].src_Y_h_max=dh_Y;
	MT_Data[0].dst_Y_h_min=0;
	MT_Data[0].dst_Y_h_max=dh_Y;
	MT_Data[0].src_UV_h_min=0;
	MT_Data[0].src_UV_h_max=dh_UV;
	MT_Data[0].dst_UV_h_min=0;
	MT_Data[0].dst_UV_h_max=dh_Y;

	i=1;
	while (i<max)
	{
		MT_Data[i].top=false;
		MT_Data[i].bottom=false;
		MT_Data[i].src_Y_h_min=MT_Data[i-1].src_Y_h_max;
		MT_Data[i].src_Y_h_max=MT_Data[i].src_Y_h_min+dh_Y;
		MT_Data[i].dst_Y_h_min=MT_Data[i-1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max=MT_Data[i].dst_Y_h_min+dh_Y;
		MT_Data[i].src_UV_h_min=MT_Data[i-1].src_UV_h_max;
		MT_Data[i].src_UV_h_max=MT_Data[i].src_UV_h_min+dh_UV;
		MT_Data[i].dst_UV_h_min=MT_Data[i-1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max=MT_Data[i].dst_UV_h_min+dh_Y;
		i++;
	}
	MT_Data[max-1].bottom=true;
	MT_Data[max-1].src_Y_h_max=size_y;
	MT_Data[max-1].dst_Y_h_max=size_y;
	MT_Data[max-1].src_UV_h_max=size_y >> 1;
	MT_Data[max-1].dst_UV_h_max=size_y;
	for (i=0; i<max; i++)
	{
		MT_Data[i].src_Y_w=size_x;
		if (output==0) MT_Data[i].dst_Y_w=size_x << 1;
		else MT_Data[i].dst_Y_w=size_x;
		MT_Data[i].src_UV_w=size_x >> 1;
		MT_Data[i].dst_UV_w=size_x >> 1;
	}
	return(max);
}


AutoYUY2::AutoYUY2(PClip _child, int _threshold, int _mode,  int _output, uint8_t _threads,bool _sleep,IScriptEnvironment* env) :
	GenericVideoFilter(_child), threshold(_threshold), mode(_mode), output(_output), threads(_threads),sleep(_sleep)
{
	bool ok;
	int16_t i,j;

	for (j=0; j<MAX_MT_THREADS; j++)
	{
		for (i=0; i<Interlaced_Tab_Size; i++)
		{
			interlaced_tab_U[j][i]=NULL;
			interlaced_tab_V[j][i]=NULL;
		}
	}
	UserId=0;

	StaticThreadpoolF=StaticThreadpool;

	for (i=0; i<MAX_MT_THREADS; i++)
	{
		MT_Thread[i].pClass=this;
		MT_Thread[i].f_process=0;
		MT_Thread[i].thread_Id=(uint8_t)i;
		MT_Thread[i].pFunc=StaticThreadpoolF;
	}

	if (vi.height<32) threads_number=1;
	else threads_number=threads;

	threads_number=CreateMTData(threads_number,vi.width,vi.height);

	if ((mode==-1) || (mode==2))
	{
		for (j=0; j<threads_number; j++)
		{
			for (i=0; i<Interlaced_Tab_Size; i++)
			{
				interlaced_tab_U[j][i]=(bool*)malloc((vi.width>>1)*sizeof(bool));
				interlaced_tab_V[j][i]=(bool*)malloc((vi.width>>1)*sizeof(bool));
			}
		}

		ok=true;
		for (j=0; j<threads_number; j++)
		{
			for (i=0; i<Interlaced_Tab_Size; i++)
			{
				ok=ok && (interlaced_tab_U[j][i]!=NULL);
				ok=ok && (interlaced_tab_V[j][i]!=NULL);
			}
		}

		if (!ok)
		{
			if (threads>1) poolInterface->DeAllocateAllThreads(true);
			FreeData();
			env->ThrowError("AutoYUY2: Memory allocation error.");
		}
	}

	switch (output)
	{
		case 0 : vi.pixel_type = VideoInfo::CS_YUY2; break;
		case 1 : vi.pixel_type = VideoInfo::CS_YV16; break;
		default : break;
	}

	for (uint16_t k=0; k<256; k++)
	{
		lookup_Upscale[k]=3*k;
		lookup_Upscale[k+256]=5*k;
		lookup_Upscale[k+512]=7*k;
	}

	SSE2_Enable=((env->GetCPUFlags()&CPUF_SSE2)!=0);
	AVX_Enable=((env->GetCPUFlags()&CPUF_AVX)!=0);

	const size_t img_size=vi.BMPSize();

	if (threads_number>1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			poolInterface->DeAllocateAllThreads(true);
			FreeData();
			env->ThrowError("AutoYUY2: Error with the TheadPool while getting UserId!");
		}
	}
}


void AutoYUY2::FreeData(void) 
{
	int16_t i,j;

	for (j=threads_number-1; j>=0; j--)
	{
		for (i=Interlaced_Tab_Size-1; i>=0; i--)
		{
			myfree(interlaced_tab_V[j][i]);
			myfree(interlaced_tab_U[j][i]);
		}
	}
}


AutoYUY2::~AutoYUY2() 
{
	if (threads_number>1) poolInterface->RemoveUserId(UserId);
	if (threads>1) poolInterface->DeAllocateAllThreads(true);
	FreeData();
}



static inline void Move_Full(const void *src_, void *dst_, const int32_t w,const int32_t h,
		int src_pitch,int dst_pitch)
{
	const uint8_t *src=(uint8_t *)src_;
	uint8_t *dst=(uint8_t *)dst_;

	if ((src_pitch==dst_pitch) && (abs(src_pitch)==w))
	{
		if (src_pitch<0)
		{
			src+=(h-1)*src_pitch;
			dst+=(h-1)*dst_pitch;
		}
		memcpy(dst,src,(size_t)h*(size_t)w);
	}
	else
	{
		for(int i=0; i<h; i++)
		{
			memcpy(dst,src,w);
			src+=src_pitch;
			dst+=dst_pitch;
		}
	}
}


void AutoYUY2::Convert_Interlaced_YV16_SSE(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Upp,*src_Un,*src_Unn,*src_Unnn;
	const uint8_t *src_V,*src_Vp,*src_Vpp,*src_Vn,*src_Vnn,*src_Vnnn;
	bool _align_U=false,_align_V=false;
	const int32_t w_U=dst_w>>1,w_V=dst_w>>1;
	const int32_t w_U4=w_U>>2,w_V4=w_V>>2;
	const int32_t w_U8=(w_U+7)>>3,w_V8=(w_V+7)>>3;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	const ptrdiff_t pitch_U_2=2*src_pitchU;
	const ptrdiff_t pitch_V_2=2*src_pitchV;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Upp=src_U-2*src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Unn=src_U+2*src_pitchU;
	src_Unnn=src_U+3*src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vpp=src_V-2*src_pitchV;
	src_Vn=src_V+src_pitchV;
	src_Vnn=src_V+2*src_pitchV;
	src_Vnnn=src_V+3*src_pitchV;

	if ((((size_t)src_U & 0x0F)==0) && ((abs(src_pitchU) & 0x0F)==0)) _align_U=true;
	if ((((size_t)src_V & 0x0F)==0) && ((abs(src_pitchV) & 0x0F)==0)) _align_V=true;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);

// U Planar
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_Un,w_U);
			dstUp+=dst_pitch_U;

			if (_align_U)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Unn,src_U,dstUp,w_U8);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_Un,src_Unnn,dstUp,w_U8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Unn,src_U,dstUp,w_U4);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_Un,src_Unnn,dstUp,w_U4);
			}
			dstUp+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		if (_align_U)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_U,src_Upp,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Up,src_Un,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Unn,src_U,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_Un,src_Unnn,dstUp,w_U8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_U,src_Upp,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Up,src_Un,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Unn,src_U,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_Un,src_Unnn,dstUp,w_U4);
		}
		dstUp+=dst_pitch_U;

		src_U+=pitch_U_2;
		src_Up+=pitch_U_2;
		src_Upp+=pitch_U_2;
		src_Un+=pitch_U_2;
		src_Unn+=pitch_U_2;
		src_Unnn+=pitch_U_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			if (_align_U)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_U,src_Upp,dstUp,w_U8);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Up,src_Un,dstUp,w_U8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_U,src_Upp,dstUp,w_U4);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Up,src_Un,dstUp,w_U4);
			}
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_Un,w_U);
			dstUp+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}


// V Planar
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_Vn,w_V);
			dstVp+=dst_pitch_V;

			if (_align_V)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Vnn,src_V,dstVp,w_V8);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_Vn,src_Vnnn,dstVp,w_V8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Vnn,src_V,dstVp,w_V4);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_Vn,src_Vnnn,dstVp,w_V4);
			}
			dstVp+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		if (_align_V)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_V,src_Vpp,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Vp,src_Vn,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Vnn,src_V,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_Vn,src_Vnnn,dstVp,w_V8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_V,src_Vpp,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Vp,src_Vn,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Vnn,src_V,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_Vn,src_Vnnn,dstVp,w_V4);
		}
		dstVp+=dst_pitch_V;

		src_V+=pitch_V_2;
		src_Vp+=pitch_V_2;
		src_Vpp+=pitch_V_2;
		src_Vn+=pitch_V_2;
		src_Vnn+=pitch_V_2;
		src_Vnnn+=pitch_V_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			if (_align_V)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b(src_V,src_Vpp,dstVp,w_V8);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b(src_Vp,src_Vn,dstVp,w_V8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3(src_V,src_Vpp,dstVp,w_V4);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2(src_Vp,src_Vn,dstVp,w_V4);
			}
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_Vn,w_V);
			dstVp+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
}


void AutoYUY2::Convert_Interlaced_YV16_AVX(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Upp,*src_Un,*src_Unn,*src_Unnn;
	const uint8_t *src_V,*src_Vp,*src_Vpp,*src_Vn,*src_Vnn,*src_Vnnn;
	bool _align_U=false,_align_V=false;
	const int32_t w_U=dst_w>>1,w_V=dst_w>>1;
	const int32_t w_U4=w_U>>2,w_V4=w_V>>2;
	const int32_t w_U8=(w_U+7)>>3,w_V8=(w_V+7)>>3;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	const ptrdiff_t pitch_U_2=2*src_pitchU;
	const ptrdiff_t pitch_V_2=2*src_pitchV;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Upp=src_U-2*src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Unn=src_U+2*src_pitchU;
	src_Unnn=src_U+3*src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vpp=src_V-2*src_pitchV;
	src_Vn=src_V+src_pitchV;
	src_Vnn=src_V+2*src_pitchV;
	src_Vnnn=src_V+3*src_pitchV;

	if ((((size_t)src_U & 0x0F)==0) && ((abs(src_pitchU) & 0x0F)==0)) _align_U=true;
	if ((((size_t)src_V & 0x0F)==0) && ((abs(src_pitchV) & 0x0F)==0)) _align_V=true;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);

// U Planar
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_Un,w_U);
			dstUp+=dst_pitch_U;

			if (_align_U)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Unn,src_U,dstUp,w_U8);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_Un,src_Unnn,dstUp,w_U8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Unn,src_U,dstUp,w_U4);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_Un,src_Unnn,dstUp,w_U4);
			}
			dstUp+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		if (_align_U)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_U,src_Upp,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Up,src_Un,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Unn,src_U,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_Un,src_Unnn,dstUp,w_U8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_U,src_Upp,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Up,src_Un,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Unn,src_U,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_Un,src_Unnn,dstUp,w_U4);
		}
		dstUp+=dst_pitch_U;

		src_U+=pitch_U_2;
		src_Up+=pitch_U_2;
		src_Upp+=pitch_U_2;
		src_Un+=pitch_U_2;
		src_Unn+=pitch_U_2;
		src_Unnn+=pitch_U_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			if (_align_U)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_U,src_Upp,dstUp,w_U8);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Up,src_Un,dstUp,w_U8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_U,src_Upp,dstUp,w_U4);
				dstUp+=dst_pitch_U;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Up,src_Un,dstUp,w_U4);
			}
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_Un,w_U);
			dstUp+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}


// V Planar
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_Vn,w_V);
			dstVp+=dst_pitch_V;

			if (_align_V)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Vnn,src_V,dstVp,w_V8);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_Vn,src_Vnnn,dstVp,w_V8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Vnn,src_V,dstVp,w_V4);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_Vn,src_Vnnn,dstVp,w_V4);
			}
			dstVp+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		if (_align_V)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_V,src_Vpp,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Vp,src_Vn,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Vnn,src_V,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_Vn,src_Vnnn,dstVp,w_V8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_V,src_Vpp,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Vp,src_Vn,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Vnn,src_V,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_Vn,src_Vnnn,dstVp,w_V4);
		}
		dstVp+=dst_pitch_V;

		src_V+=pitch_V_2;
		src_Vp+=pitch_V_2;
		src_Vpp+=pitch_V_2;
		src_Vn+=pitch_V_2;
		src_Vnn+=pitch_V_2;
		src_Vnnn+=pitch_V_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			if (_align_V)
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b(src_V,src_Vpp,dstVp,w_V8);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b(src_Vp,src_Vn,dstVp,w_V8);
			}
			else
			{
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3(src_V,src_Vpp,dstVp,w_V4);
				dstVp+=dst_pitch_V;
				JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2(src_Vp,src_Vn,dstVp,w_V4);
			}
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_Vn,w_V);
			dstVp+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
}


void AutoYUY2::Convert_Interlaced_YV16(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Upp,*src_Un,*src_Unn,*src_Unnn;
	const uint8_t *src_V,*src_Vp,*src_Vpp,*src_Vn,*src_Vnn,*src_Vnnn;
	uint8_t *dst_V,*dst_U;
	const int32_t w_U=dst_w>>1,w_V=dst_w>>1;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	const ptrdiff_t pitch_U_2=2*src_pitchU;
	const ptrdiff_t pitch_V_2=2*src_pitchV;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Upp=src_U-2*src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Unn=src_U+2*src_pitchU;
	src_Unnn=src_U+3*src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vpp=src_V-2*src_pitchV;
	src_Vn=src_V+src_pitchV;
	src_Vnn=src_V+2*src_pitchV;
	src_Vnnn=src_V+3*src_pitchV;
	dst_U=dstUp;
	dst_V=dstVp;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);


// Planar U
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dst_U,src_U,w_U);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_Un,w_U);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_U; j++)
				dst_U[j]=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_U; j++)
				dst_U[j]=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);
			dst_U+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		for(int32_t j=0; j<w_U; j++)
			dst_U[j]=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4)>>3);
		dst_U+=dst_pitch_U;

		for(int32_t j=0; j<w_U; j++)
			dst_U[j]=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4)>>3);
		dst_U+=dst_pitch_U;

		for(int32_t j=0; j<w_U; j++)
			dst_U[j]=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
		dst_U+=dst_pitch_U;

		for(int32_t j=0; j<w_U; j++)
			dst_U[j]=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);
		dst_U+=dst_pitch_U;

		src_U+=pitch_U_2;
		src_Up+=pitch_U_2;
		src_Upp+=pitch_U_2;
		src_Un+=pitch_U_2;
		src_Unn+=pitch_U_2;
		src_Unnn+=pitch_U_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			for(int32_t j=0; j<w_U; j++)
				dst_U[j]=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4)>>3);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_U; j++)
				dst_U[j]=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4)>>3);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_U,w_U);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_Un,w_U);
			dst_U+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}


// Planar V
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dst_V,src_V,w_V);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_Vn,w_V);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_V; j++)
				dst_V[j]=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_V; j++)
				dst_V[j]=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
			dst_V+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		for(int32_t j=0; j<w_V; j++)
			dst_V[j]=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4)>>3);
		dst_V+=dst_pitch_V;

		for(int32_t j=0; j<w_V; j++)
			dst_V[j]=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4)>>3);
		dst_V+=dst_pitch_V;

		for(int32_t j=0; j<w_V; j++)
			dst_V[j]=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
		dst_V+=dst_pitch_V;

		for(int32_t j=0; j<w_V; j++)
			dst_V[j]=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
		dst_V+=dst_pitch_V;

		src_V+=pitch_V_2;
		src_Vp+=pitch_V_2;
		src_Vpp+=pitch_V_2;
		src_Vn+=pitch_V_2;
		src_Vnn+=pitch_V_2;
		src_Vnnn+=pitch_V_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			for(int32_t j=0; j<w_V; j++)
				dst_V[j]=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4)>>3);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_V; j++)
				dst_V[j]=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4)>>3);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_V,w_V);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_Vn,w_V);
			dst_V+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
}


void AutoYUY2::Convert_Progressive_YV16_SSE(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Un;
	const uint8_t *src_V,*src_Vp,*src_Vn;
	bool _align_U=false,_align_V=false;
	const int32_t w_U=dst_w>>1,w_V=dst_w>>1;
	const int32_t w_U4=w_U>>2,w_V4=w_V>>2;
	const int32_t w_U8=(w_U+7)>>3,w_V8=(w_V+7)>>3;
	const int32_t h_2 = mt_data_inf.bottom ? h_Y_max-2:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 2:h_Y_min;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vn=src_V+src_pitchV;

	if ((((size_t)src_U & 0x0F)==0) && ((abs(src_pitchU) & 0x0F)==0)) _align_U=true;
	if ((((size_t)src_V & 0x0F)==0) && ((abs(src_pitchV) & 0x0F)==0)) _align_V=true;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);

// Planar U
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<2; i+=2)
		{
			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			if (_align_U) JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_U,src_Un,dstUp,w_U8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_U,src_Un,dstUp,w_U4);
			dstUp+=dst_pitch_U;

			src_U+=src_pitchU;
			src_Up+=src_pitchU;
			src_Un+=src_pitchU;
		}
	}

	for(int32_t i=h_0; i<h_2; i+=2)
	{
		if (_align_U)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_U,src_Up,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_U,src_Un,dstUp,w_U8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_U,src_Up,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_U,src_Un,dstUp,w_U4);
		}
		dstUp+=dst_pitch_U;

		src_U+=src_pitchU;
		src_Up+=src_pitchU;
		src_Un+=src_pitchU;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_2; i<h_Y_max; i+=2)
		{
			if (_align_U) JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_U,src_Up,dstUp,w_U8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_U,src_Up,dstUp,w_U4);
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			src_U+=src_pitchU;
			src_Up+=src_pitchU;
			src_Un+=src_pitchU;
		}
	}


// Planar V
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<2; i+=2)
		{
			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			if (_align_V) JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_V,src_Vn,dstVp,w_V8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_V,src_Vn,dstVp,w_V4);
			dstVp+=dst_pitch_V;

			src_V+=src_pitchV;
			src_Vp+=src_pitchV;
			src_Vn+=src_pitchV;
		}
	}

	for(int32_t i=h_0; i<h_2; i+=2)
	{
		if (_align_V)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_V,src_Vp,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_V,src_Vn,dstVp,w_V8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_V,src_Vp,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_V,src_Vn,dstVp,w_V4);
		}
		dstVp+=dst_pitch_V;

		src_V+=src_pitchV;
		src_Vp+=src_pitchV;
		src_Vn+=src_pitchV;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_2; i<h_Y_max; i+=2)
		{
			if (_align_V) JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b(src_V,src_Vp,dstVp,w_V8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4(src_V,src_Vp,dstVp,w_V4);
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			src_V+=src_pitchV;
			src_Vp+=src_pitchV;
			src_Vn+=src_pitchV;
		}
	}
}


void AutoYUY2::Convert_Progressive_YV16_AVX(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Un;
	const uint8_t *src_V,*src_Vp,*src_Vn;
	bool _align_U=false,_align_V=false;
	const int32_t w_U=dst_w>>1,w_V=dst_w>>1;
	const int32_t w_U4=w_U>>2,w_V4=w_V>>2;
	const int32_t w_U8=(w_U+7)>>3,w_V8=(w_V+7)>>3;
	const int32_t h_2 = mt_data_inf.bottom ? h_Y_max-2:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 2:h_Y_min;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vn=src_V+src_pitchV;

	if ((((size_t)src_U & 0x0F)==0) && ((abs(src_pitchU) & 0x0F)==0)) _align_U=true;
	if ((((size_t)src_V & 0x0F)==0) && ((abs(src_pitchV) & 0x0F)==0)) _align_V=true;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);

// Planar U
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<2; i+=2)
		{
			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			if (_align_U) JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_U,src_Un,dstUp,w_U8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_U,src_Un,dstUp,w_U4);
			dstUp+=dst_pitch_U;

			src_U+=src_pitchU;
			src_Up+=src_pitchU;
			src_Un+=src_pitchU;
		}
	}

	for(int32_t i=h_0; i<h_2; i+=2)
	{
		if (_align_U)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_U,src_Up,dstUp,w_U8);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_U,src_Un,dstUp,w_U8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_U,src_Up,dstUp,w_U4);
			dstUp+=dst_pitch_U;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_U,src_Un,dstUp,w_U4);
		}
		dstUp+=dst_pitch_U;

		src_U+=src_pitchU;
		src_Up+=src_pitchU;
		src_Un+=src_pitchU;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_2; i<h_Y_max; i+=2)
		{
			if (_align_U) JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_U,src_Up,dstUp,w_U8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_U,src_Up,dstUp,w_U4);
			dstUp+=dst_pitch_U;

			memcpy(dstUp,src_U,w_U);
			dstUp+=dst_pitch_U;

			src_U+=src_pitchU;
			src_Up+=src_pitchU;
			src_Un+=src_pitchU;
		}
	}


// Planar V
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<2; i+=2)
		{
			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			if (_align_V) JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_V,src_Vn,dstVp,w_V8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_V,src_Vn,dstVp,w_V4);
			dstVp+=dst_pitch_V;

			src_V+=src_pitchV;
			src_Vp+=src_pitchV;
			src_Vn+=src_pitchV;
		}
	}

	for(int32_t i=h_0; i<h_2; i+=2)
	{
		if (_align_V)
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_V,src_Vp,dstVp,w_V8);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_V,src_Vn,dstVp,w_V8);
		}
		else
		{
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_V,src_Vp,dstVp,w_V4);
			dstVp+=dst_pitch_V;
			JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_V,src_Vn,dstVp,w_V4);
		}
		dstVp+=dst_pitch_V;

		src_V+=src_pitchV;
		src_Vp+=src_pitchV;
		src_Vn+=src_pitchV;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_2; i<h_Y_max; i+=2)
		{
			if (_align_V) JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b(src_V,src_Vp,dstVp,w_V8);
			else JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4(src_V,src_Vp,dstVp,w_V4);
			dstVp+=dst_pitch_V;

			memcpy(dstVp,src_V,w_V);
			dstVp+=dst_pitch_V;

			src_V+=src_pitchV;
			src_Vp+=src_pitchV;
			src_Vn+=src_pitchV;
		}
	}
}


void AutoYUY2::Convert_Progressive_YV16(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Un;
	const uint8_t *src_V,*src_Vp,*src_Vn;
	uint8_t *dst_U,*dst_V;
	const int32_t w_U=dst_w>>1,w_V=dst_w>>1;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_2 = mt_data_inf.bottom ? h_Y_max-2:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 2:h_Y_min;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vn=src_V+src_pitchV;
	dst_U=dstUp;
	dst_V=dstVp;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);


// Planar U
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<2; i+=2)
		{
			memcpy(dst_U,src_U,w_U);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_U; j++)
				dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2)>>2);
			dst_U+=dst_pitch_U;

			src_U+=src_pitchU;
			src_Up+=src_pitchU;
			src_Un+=src_pitchU;
		}
	}

	for(int32_t i=h_0; i<h_2; i+=2)
	{
		for(int32_t j=0; j<w_U; j++)
			dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2)>>2);
		dst_U+=dst_pitch_U;

		for(int32_t j=0; j<w_U; j++)
			dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2)>>2);
		dst_U+=dst_pitch_U;

		src_U+=src_pitchU;
		src_Up+=src_pitchU;
		src_Un+=src_pitchU;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_2; i<h_Y_max; i+=2)
		{
			for(int32_t j=0; j<w_U; j++)
				dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2)>>2);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_U,w_U);
			dst_U+=dst_pitch_U;

			src_U+=src_pitchU;
			src_Up+=src_pitchU;
			src_Un+=src_pitchU;
		}
	}


// Planar V
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<2; i+=2)
		{
			memcpy(dst_V,src_V,w_V);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_V; j++)
				dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2)>>2);
			dst_V+=dst_pitch_V;

			src_V+=src_pitchV;
			src_Vp+=src_pitchV;
			src_Vn+=src_pitchV;
		}
	}

	for(int32_t i=h_0; i<h_2; i+=2)
	{
		for(int32_t j=0; j<w_V; j++)
			dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2)>>2);
		dst_V+=dst_pitch_V;

		for(int32_t j=0; j<w_V; j++)
			dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2)>>2);
		dst_V+=dst_pitch_V;

		src_V+=src_pitchV;
		src_Vp+=src_pitchV;
		src_Vn+=src_pitchV;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_2; i<h_Y_max; i+=2)
		{
			for(int32_t j=0; j<w_V; j++)
				dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2)>>2);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_V,w_V);
			dst_V+=dst_pitch_V;

			src_V+=src_pitchV;
			src_Vp+=src_pitchV;
			src_Vn+=src_pitchV;
		}
	}

}



void AutoYUY2::Convert_Automatic_YV16(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Upp,*src_Un,*src_Unn,*src_Unnn;
	const uint8_t *src_V,*src_Vp,*src_Vpp,*src_Vn,*src_Vnn,*src_Vnnn;
	uint8_t *dst_U,*dst_V;
	uint8_t index_tab_0,index_tab_1,index_tab_2;
	const int16_t threshold_=threshold;
	const int32_t w_UV=dst_w>>1;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	const ptrdiff_t pitch_U_2=2*src_pitchU;
	const ptrdiff_t pitch_V_2=2*src_pitchV;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Upp=src_U-2*src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Unn=src_U+2*src_pitchU;
	src_Unnn=src_U+3*src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vpp=src_V-2*src_pitchV;
	src_Vn=src_V+src_pitchV;
	src_Vnn=src_V+2*src_pitchV;
	src_Vnnn=src_V+3*src_pitchV;
	dst_U=dstUp;
	dst_V=dstVp;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);

// Planar U	
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dst_U,src_U,w_UV);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_Un,w_UV);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_UV; j++)
				dst_U[j]=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
			dst_U+=dst_pitch_U;

			{
				bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
				
				for(int32_t j=0; j<w_UV; j++)
				{
					if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
						(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
						(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
						((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
						itabu0[j]=true;
					else itabu0[j]=false;
					if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
						((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
						itabu1[j]=true;
					else itabu1[j]=false;
					
					dst_U[j]=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);
				}
			}
			dst_U+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}
	else
	{
		bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
				
		for(int32_t j=0; j<w_UV; j++)
		{
			if (((abs((int16_t)src_Upp[j]-(int16_t)src_Up[j])>=threshold_) &&
				(abs((int16_t)src_U[j]-(int16_t)src_Up[j])>=threshold_) &&
				(((src_Upp[j]>src_Up[j]) && (src_U[j]>src_Up[j])) ||
				((src_Upp[j]<src_Up[j]) && (src_U[j]<src_Up[j])))))
				itabu0[j]=true;
			else itabu0[j]=false;
			if (((abs((int16_t)src_Up[j]-(int16_t)src_U[j])>=threshold_) &&
				(abs((int16_t)src_Un[j]-(int16_t)src_U[j])>=threshold_) &&
				(((src_Up[j]>src_U[j]) && (src_Un[j]>src_U[j])) ||
				((src_Up[j]<src_U[j]) && (src_Un[j]<src_U[j])))))
				itabu1[j]=true;
			else itabu1[j]=false;
		}

	}

	index_tab_0=0;
	index_tab_1=1;
	index_tab_2=2;

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		{
			const bool *itabu0=interlaced_tab_U[thread_num][index_tab_0],*itabu1=interlaced_tab_U[thread_num][index_tab_1];

			// Upsample as needed.
			for(int32_t j=0; j<w_UV; j++)
			{
				if ((itabu0[j]) || (itabu1[j]))
				{
					dst_U[j]=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4)>>3);
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2) >> 2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		{
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1];
			bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
					(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
					(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
					((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
					itabu2[j]=true;
				else itabu2[j]=false;			

				// Upsample as needed.
				if ((itabu2[j]) || (itabu1[j]))
				{
					dst_U[j]=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4)>>3);
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2)>>2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		{
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1],*itabu2=interlaced_tab_U[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				// Upsample as needed.
				if ((itabu1[j]) || (itabu2[j]))
				{
					dst_U[j]=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_U[j]+2)>>2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		{
			bool *itabu0=interlaced_tab_U[thread_num][index_tab_0];
			const bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
					((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
					itabu0[j]=true;
				else itabu0[j]=false;

				// Upsample as needed.
				if ((itabu0[j]) || (itabu2[j]))
				{
					dst_U[j]=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_Unn[j]+2)>>2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		index_tab_0=(index_tab_0+2)%Interlaced_Tab_Size;
		index_tab_1=(index_tab_1+2)%Interlaced_Tab_Size;
		index_tab_2=(index_tab_2+2)%Interlaced_Tab_Size;

		src_U+=pitch_U_2;
		src_Up+=pitch_U_2;
		src_Upp+=pitch_U_2;
		src_Un+=pitch_U_2;
		src_Unn+=pitch_U_2;
		src_Unnn+=pitch_U_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			for(int32_t j=0; j<w_UV; j++)
				dst_U[j]=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4)>>3);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_UV; j++)
				dst_U[j]=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4)>>3);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_U,w_UV);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_Un,w_UV);
			dst_U+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}


// Planar V
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dst_V,src_V,w_UV);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_Vn,w_UV);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_UV; j++)
				dst_V[j]=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
			dst_V+=dst_pitch_V;

			{
				bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

				for(int32_t j=0; j<w_UV; j++)
				{
					if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
						((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
						itabv0[j]=true;
					else itabv0[j]=false;
					if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
						((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
						itabv1[j]=true;
					else itabv1[j]=false;
					
					dst_V[j]=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
				}
			}
			dst_V+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
	else
	{
		bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

		for(int32_t j=0; j<w_UV; j++)
		{
			if (((abs((int16_t)src_Vpp[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(abs((int16_t)src_V[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(((src_Vpp[j]>src_Vp[j]) && (src_V[j]>src_Vp[j])) ||
				((src_Vpp[j]<src_Vp[j]) && (src_V[j]<src_Vp[j])))))
				itabv0[j]=true;
			else itabv0[j]=false;
			if (((abs((int16_t)src_Vp[j]-(int16_t)src_V[j])>=threshold_) &&
				(abs((int16_t)src_Vn[j]-(int16_t)src_V[j])>=threshold_) &&
				(((src_Vp[j]>src_V[j]) && (src_Vn[j]>src_V[j])) ||
				((src_Vp[j]<src_V[j]) && (src_Vn[j]<src_V[j])))))
				itabv1[j]=true;
			else itabv1[j]=false;
		}
	}
	
	index_tab_0=0;
	index_tab_1=1;
	index_tab_2=2;

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		{
			const bool *itabv0=interlaced_tab_V[thread_num][index_tab_0],*itabv1=interlaced_tab_V[thread_num][index_tab_1];

			for(int32_t j=0; j<w_UV; j++)
			{
				// Upsample as needed.
				if ((itabv0[j]) || (itabv1[j]))
				{
					dst_V[j]=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4)>>3);
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2) >> 2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		{
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1];
			bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
					((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
					itabv2[j]=true;
				else itabv2[j]=false;			

				// Upsample as needed.
				if ((itabv2[j]) || (itabv1[j]))
				{
					dst_V[j]=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4)>>3);
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2)>>2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		{
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1],*itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				// Upsample as needed.
				if ((itabv1[j]) || (itabv2[j]))
				{
					dst_V[j]=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_V[j]+2)>>2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		{
			bool *itabv0=interlaced_tab_V[thread_num][index_tab_0];
			const bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
					((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
					itabv0[j]=true;
				else itabv0[j]=false;

				// Upsample as needed.
				if ((itabv0[j]) || (itabv2[j]))
				{
					dst_V[j]=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_Vnn[j]+2)>>2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		index_tab_0=(index_tab_0+2)%Interlaced_Tab_Size;
		index_tab_1=(index_tab_1+2)%Interlaced_Tab_Size;
		index_tab_2=(index_tab_2+2)%Interlaced_Tab_Size;

		src_V+=pitch_V_2;
		src_Vp+=pitch_V_2;
		src_Vpp+=pitch_V_2;
		src_Vn+=pitch_V_2;
		src_Vnn+=pitch_V_2;
		src_Vnnn+=pitch_V_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			for(int32_t j=0; j<w_UV; j++)
				dst_V[j]=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4)>>3);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_UV; j++)
				dst_V[j]=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4)>>3);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_V,w_UV);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_Vn,w_UV);
			dst_V+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
}



void AutoYUY2::Convert_Test_YV16(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstYp=(uint8_t *)mt_data_inf.dst1;
	uint8_t *dstUp=(uint8_t *)mt_data_inf.dst2;
	uint8_t *dstVp=(uint8_t *)mt_data_inf.dst3;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitch_Y=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch_Y=mt_data_inf.dst_pitch1;
	const int dst_pitch_U=mt_data_inf.dst_pitch2;
	const int dst_pitch_V=mt_data_inf.dst_pitch3;

	const uint8_t *src_U,*src_Up,*src_Upp,*src_Un,*src_Unn,*src_Unnn;
	const uint8_t *src_V,*src_Vp,*src_Vpp,*src_Vn,*src_Vnn,*src_Vnnn;
	uint8_t *dst_U,*dst_V;
	uint8_t index_tab_0,index_tab_1,index_tab_2;
	const int16_t threshold_=threshold;
	const int32_t w_UV=dst_w>>1;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	const ptrdiff_t pitch_U_2=2*src_pitchU;
	const ptrdiff_t pitch_V_2=2*src_pitchV;

	src_U=srcUp;
	src_V=srcVp;
	src_Up=src_U-src_pitchU;
	src_Upp=src_U-2*src_pitchU;
	src_Un=src_U+src_pitchU;
	src_Unn=src_U+2*src_pitchU;
	src_Unnn=src_U+3*src_pitchU;
	src_Vp=src_V-src_pitchV;
	src_Vpp=src_V-2*src_pitchV;
	src_Vn=src_V+src_pitchV;
	src_Vnn=src_V+2*src_pitchV;
	src_Vnnn=src_V+3*src_pitchV;
	dst_U=dstUp;
	dst_V=dstVp;

	Move_Full(srcYp,dstYp,dst_w,h_Y_max-h_Y_min,src_pitch_Y,dst_pitch_Y);

// Planar U	
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dst_U,src_U,w_UV);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_Un,w_UV);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_UV; j++)
				dst_U[j]=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
			dst_U+=dst_pitch_U;

			{
				bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
				
				for(int32_t j=0; j<w_UV; j++)
				{
					if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
						(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
						(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
						((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
						itabu0[j]=true;
					else itabu0[j]=false;
					if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
						((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
						itabu1[j]=true;
					else itabu1[j]=false;
					
					dst_U[j]=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);
				}
			}
			dst_U+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}
	else
	{
		bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
				
		for(int32_t j=0; j<w_UV; j++)
		{
			if (((abs((int16_t)src_Upp[j]-(int16_t)src_Up[j])>=threshold_) &&
				(abs((int16_t)src_U[j]-(int16_t)src_Up[j])>=threshold_) &&
				(((src_Upp[j]>src_Up[j]) && (src_U[j]>src_Up[j])) ||
				((src_Upp[j]<src_Up[j]) && (src_U[j]<src_Up[j])))))
				itabu0[j]=true;
			else itabu0[j]=false;
			if (((abs((int16_t)src_Up[j]-(int16_t)src_U[j])>=threshold_) &&
				(abs((int16_t)src_Un[j]-(int16_t)src_U[j])>=threshold_) &&
				(((src_Up[j]>src_U[j]) && (src_Un[j]>src_U[j])) ||
				((src_Up[j]<src_U[j]) && (src_Un[j]<src_U[j])))))
				itabu1[j]=true;
			else itabu1[j]=false;
		}

	}

	index_tab_0=0;
	index_tab_1=1;
	index_tab_2=2;

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		{
			const bool *itabu0=interlaced_tab_U[thread_num][index_tab_0],*itabu1=interlaced_tab_U[thread_num][index_tab_1];

			// Upsample as needed.
			for(int32_t j=0; j<w_UV; j++)
			{
				if ((itabu0[j]) || (itabu1[j]))
				{
					dst_U[j]=239;
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2) >> 2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		{
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1];
			bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
					(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
					(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
					((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
					itabu2[j]=true;
				else itabu2[j]=false;			

				// Upsample as needed.
				if ((itabu2[j]) || (itabu1[j]))
				{
					dst_U[j]=239;
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2)>>2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		{
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1],*itabu2=interlaced_tab_U[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				// Upsample as needed.
				if ((itabu1[j]) || (itabu2[j]))
				{
					dst_U[j]=239;
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_U[j]+2)>>2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		{
			bool *itabu0=interlaced_tab_U[thread_num][index_tab_0];
			const bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
					((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
					itabu0[j]=true;
				else itabu0[j]=false;

				// Upsample as needed.
				if ((itabu0[j]) || (itabu2[j]))
				{
					dst_U[j]=239;
				}
				else
				{
					dst_U[j]=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_Unn[j]+2)>>2);
				}
			}
		}
		dst_U+=dst_pitch_U;

		index_tab_0=(index_tab_0+2)%Interlaced_Tab_Size;
		index_tab_1=(index_tab_1+2)%Interlaced_Tab_Size;
		index_tab_2=(index_tab_2+2)%Interlaced_Tab_Size;

		src_U+=pitch_U_2;
		src_Up+=pitch_U_2;
		src_Upp+=pitch_U_2;
		src_Un+=pitch_U_2;
		src_Unn+=pitch_U_2;
		src_Unnn+=pitch_U_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			for(int32_t j=0; j<w_UV; j++)
				dst_U[j]=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4)>>3);
			dst_U+=dst_pitch_U;

			for(int32_t j=0; j<w_UV; j++)
				dst_U[j]=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4)>>3);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_U,w_UV);
			dst_U+=dst_pitch_U;

			memcpy(dst_U,src_Un,w_UV);
			dst_U+=dst_pitch_U;

			src_U+=pitch_U_2;
			src_Up+=pitch_U_2;
			src_Upp+=pitch_U_2;
			src_Un+=pitch_U_2;
			src_Unn+=pitch_U_2;
			src_Unnn+=pitch_U_2;
		}
	}


// Planar V
	if (mt_data_inf.top)
	{
		for(int32_t i=0; i<4; i+=4)
		{
			memcpy(dst_V,src_V,w_UV);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_Vn,w_UV);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_UV; j++)
				dst_V[j]=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
			dst_V+=dst_pitch_V;

			{
				bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

				for(int32_t j=0; j<w_UV; j++)
				{
					if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
						((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
						itabv0[j]=true;
					else itabv0[j]=false;
					if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
						((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
						itabv1[j]=true;
					else itabv1[j]=false;
					
					dst_V[j]=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
				}
			}
			dst_V+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
	else
	{
		bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

		for(int32_t j=0; j<w_UV; j++)
		{
			if (((abs((int16_t)src_Vpp[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(abs((int16_t)src_V[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(((src_Vpp[j]>src_Vp[j]) && (src_V[j]>src_Vp[j])) ||
				((src_Vpp[j]<src_Vp[j]) && (src_V[j]<src_Vp[j])))))
				itabv0[j]=true;
			else itabv0[j]=false;
			if (((abs((int16_t)src_Vp[j]-(int16_t)src_V[j])>=threshold_) &&
				(abs((int16_t)src_Vn[j]-(int16_t)src_V[j])>=threshold_) &&
				(((src_Vp[j]>src_V[j]) && (src_Vn[j]>src_V[j])) ||
				((src_Vp[j]<src_V[j]) && (src_Vn[j]<src_V[j])))))
				itabv1[j]=true;
			else itabv1[j]=false;
		}
	}
	
	index_tab_0=0;
	index_tab_1=1;
	index_tab_2=2;

	for(int32_t i=h_0; i<h_4; i+=4)
	{
		{
			const bool *itabv0=interlaced_tab_V[thread_num][index_tab_0],*itabv1=interlaced_tab_V[thread_num][index_tab_1];

			for(int32_t j=0; j<w_UV; j++)
			{
				// Upsample as needed.
				if ((itabv0[j]) || (itabv1[j]))
				{
					dst_V[j]=239;
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2) >> 2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		{
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1];
			bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
					((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
					itabv2[j]=true;
				else itabv2[j]=false;			

				// Upsample as needed.
				if ((itabv2[j]) || (itabv1[j]))
				{
					dst_V[j]=239;
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2)>>2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		{
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1],*itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				// Upsample as needed.
				if ((itabv1[j]) || (itabv2[j]))
				{
					dst_V[j]=239;
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_V[j]+2)>>2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		{
			bool *itabv0=interlaced_tab_V[thread_num][index_tab_0];
			const bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for(int32_t j=0; j<w_UV; j++)
			{
				if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
					((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
					itabv0[j]=true;
				else itabv0[j]=false;

				// Upsample as needed.
				if ((itabv0[j]) || (itabv2[j]))
				{
					dst_V[j]=239;
				}
				else
				{
					dst_V[j]=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_Vnn[j]+2)>>2);
				}
			}
		}
		dst_V+=dst_pitch_V;

		index_tab_0=(index_tab_0+2)%Interlaced_Tab_Size;
		index_tab_1=(index_tab_1+2)%Interlaced_Tab_Size;
		index_tab_2=(index_tab_2+2)%Interlaced_Tab_Size;

		src_V+=pitch_V_2;
		src_Vp+=pitch_V_2;
		src_Vpp+=pitch_V_2;
		src_Vn+=pitch_V_2;
		src_Vnn+=pitch_V_2;
		src_Vnnn+=pitch_V_2;
	}

	if (mt_data_inf.bottom)
	{
		for(int32_t i=h_4; i<h_Y_max; i+=4)
		{
			for(int32_t j=0; j<w_UV; j++)
				dst_V[j]=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4)>>3);
			dst_V+=dst_pitch_V;

			for(int32_t j=0; j<w_UV; j++)
				dst_V[j]=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4)>>3);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_V,w_UV);
			dst_V+=dst_pitch_V;

			memcpy(dst_V,src_Vn,w_UV);
			dst_V+=dst_pitch_V;

			src_V+=pitch_V_2;
			src_Vp+=pitch_V_2;
			src_Vpp+=pitch_V_2;
			src_Vn+=pitch_V_2;
			src_Vnn+=pitch_V_2;
			src_Vnnn+=pitch_V_2;
		}
	}
}



void AutoYUY2::Convert_Progressive_YUY2(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	YUYV *dst=(YUYV *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up, *src_Un;
	const uint8_t *src_V,*src_Vp, *src_Vn;
	const int32_t w_UV=dst_w>>2;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_2 = mt_data_inf.bottom ? h_Y_max-2:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 2:h_Y_min;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Un = srcUp + src_pitchU;
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vn = srcVp + src_pitchV;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<2; y+=2)
		{

			JPSDR_AutoYUY2_1(src_Y,src_U,src_V,(uint8_t *)dst,dst_w>>2);
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2) >> 2);
					dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2) >> 2);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU;
			src_Up += src_pitchU;
			src_Un += src_pitchU;
			src_V += src_pitchV;
			src_Vp += src_pitchV;
			src_Vn += src_pitchV;
		}
	}

	for (int32_t y = h_0; y < h_2; y+=2)
	{
		{
			int32_t i=0;

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;
				dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2) >> 2);
				dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2) >> 2);
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;
				dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2) >> 2);
				dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2) >> 2);
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;
		
		src_U += src_pitchU;
		src_Up += src_pitchU;
		src_Un += src_pitchU;
		src_V += src_pitchV;
		src_Vp += src_pitchV;
		src_Vn += src_pitchV;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_2; y < h_Y_max; y+=2)
		{
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2) >> 2);
					dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2) >> 2);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			JPSDR_AutoYUY2_1(src_Y,src_U,src_V,(uint8_t *)dst,dst_w>>2);
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;
		
			src_U += src_pitchU;
			src_Up += src_pitchU;
			src_Un += src_pitchU;
			src_V += src_pitchV;
			src_Vp += src_pitchV;
			src_Vn += src_pitchV;	
		}
	}
}



void AutoYUY2::Convert_Progressive_YUY2_SSE(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstp=(uint8_t *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up, *src_Un;
	const uint8_t *src_V,*src_Vp, *src_Vn;
	const int32_t w_8=(dst_w+15)>>4;
	const int32_t h_2 = mt_data_inf.bottom ? h_Y_max-2:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 2:h_Y_min;

	bool _align=false;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Un = srcUp + src_pitchU;
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vn = srcVp + src_pitchV;

	if ((((size_t)srcYp & 0x0F)==0) && ((abs(src_pitchY) & 0x0F)==0)) _align=true;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<2; y+=2)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_SSE2_1b(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_4b(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_SSE2_1(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_4(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;

			src_U += src_pitchU;
			src_Up += src_pitchU;
			src_Un += src_pitchU;
			src_V += src_pitchV;
			src_Vp += src_pitchV;
			src_Vn += src_pitchV;
		}
	}

	for (int32_t y = h_0; y < h_2; y+=2)
	{
		if (_align)
		{
			JPSDR_AutoYUY2_SSE2_4b(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_4b(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
		}
		else
		{
			JPSDR_AutoYUY2_SSE2_4(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_4(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
		}
		dstp += dst_pitch;
		src_Y += src_pitchY;
		
		src_U += src_pitchU;
		src_Up += src_pitchU;
		src_Un += src_pitchU;
		src_V += src_pitchV;
		src_Vp += src_pitchV;
		src_Vn += src_pitchV;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_2; y < h_Y_max; y+=2)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_SSE2_4b(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1b(src_Y,src_U,src_V,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_SSE2_4(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1(src_Y,src_U,src_V,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;
		
			src_U += src_pitchU;
			src_Up += src_pitchU;
			src_Un += src_pitchU;
			src_V += src_pitchV;
			src_Vp += src_pitchV;
			src_Vn += src_pitchV;	
		}
	}
}


void AutoYUY2::Convert_Progressive_YUY2_AVX(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstp=(uint8_t *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up, *src_Un;
	const uint8_t *src_V,*src_Vp, *src_Vn;
	const int32_t w_8=(dst_w+15)>>4;
	const int32_t h_2 = mt_data_inf.bottom ? h_Y_max-2:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 2:h_Y_min;

	bool _align=false;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Un = srcUp + src_pitchU;
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vn = srcVp + src_pitchV;

	if ((((size_t)srcYp & 0x0F)==0) && ((abs(src_pitchY) & 0x0F)==0)) _align=true;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<2; y+=2)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_AVX_1b(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_4b(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_AVX_1(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_4(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;

			src_U += src_pitchU;
			src_Up += src_pitchU;
			src_Un += src_pitchU;
			src_V += src_pitchV;
			src_Vp += src_pitchV;
			src_Vn += src_pitchV;
		}
	}

	for (int32_t y = h_0; y < h_2; y+=2)
	{
		if (_align)
		{
			JPSDR_AutoYUY2_AVX_4b(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_4b(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
		}
		else
		{
			JPSDR_AutoYUY2_AVX_4(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_4(src_Y,src_U,src_Un,src_V,src_Vn,dstp,w_8);
		}
		dstp += dst_pitch;
		src_Y += src_pitchY;
		
		src_U += src_pitchU;
		src_Up += src_pitchU;
		src_Un += src_pitchU;
		src_V += src_pitchV;
		src_Vp += src_pitchV;
		src_Vn += src_pitchV;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_2; y < h_Y_max; y+=2)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_AVX_4b(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1b(src_Y,src_U,src_V,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_AVX_4(src_Y,src_U,src_Up,src_V,src_Vp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1(src_Y,src_U,src_V,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;
		
			src_U += src_pitchU;
			src_Up += src_pitchU;
			src_Un += src_pitchU;
			src_V += src_pitchV;
			src_Vp += src_pitchV;
			src_Vn += src_pitchV;	
		}
	}
}


void AutoYUY2::Convert_Interlaced_YUY2(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	YUYV *dst=(YUYV *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up, *src_Un, *src_Unn, *src_Unnn,*src_Upp;
	const uint8_t *src_V,*src_Vp, *src_Vn, *src_Vnn, *src_Vnnn,*src_Vpp;
	const int32_t w_UV=dst_w>>2;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Upp = srcUp - (2*src_pitchU);
	src_Un = srcUp + src_pitchU;
	src_Unn = srcUp + (2 * src_pitchU);
	src_Unnn = srcUp + (3 * src_pitchU);
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vpp = srcVp - (2*src_pitchV);
	src_Vn = srcVp + src_pitchV;
	src_Vnn = srcVp + (2 * src_pitchV);
	src_Vnnn = srcVp + (3 * src_pitchV);
	
	const ptrdiff_t src_pitchU_2=src_pitchU << 1;
	const ptrdiff_t src_pitchV_2=src_pitchV << 1;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<4; y+=4)
		{

			JPSDR_AutoYUY2_1(src_Y,src_U,src_V,(uint8_t *)dst,dst_w>>2);
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			JPSDR_AutoYUY2_1(src_Y,src_Un,src_Vn,(uint8_t *)dst,dst_w>>2);
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
					dst[j].v=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}

	for (int32_t y = h_0; y < h_4; y+=4)
	{
		{
			int32_t i=0;

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;
				dst[j].u=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4) >> 3);
				dst[j].v=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4) >> 3);
			}
		}

		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;
				dst[j].u=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4) >> 3);
				dst[j].v=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4) >> 3);
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;
		
		{
			int32_t i=0;

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;
				dst[j].u=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4) >> 3);
				dst[j].v=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4) >> 3);
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;
				dst[j].u=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4) >> 3);
				dst[j].v=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4) >> 3);
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		src_U += src_pitchU_2;
		src_Up += src_pitchU_2;
		src_Upp += src_pitchU_2;
		src_Un += src_pitchU_2;
		src_Unn += src_pitchU_2;
		src_Unnn += src_pitchU_2;
		src_V += src_pitchV_2;
		src_Vp += src_pitchV_2;
		src_Vpp += src_pitchV_2;
		src_Vn += src_pitchV_2;
		src_Vnn += src_pitchV_2;
		src_Vnnn += src_pitchV_2;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_4; y < h_Y_max; y+=4)
		{
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4) >> 3);
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4) >> 3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;
			
				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4) >> 3);
					dst[j].v=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4) >> 3);
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;
		
			JPSDR_AutoYUY2_1(src_Y,src_U,src_V,(uint8_t *)dst,dst_w>>2);
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			JPSDR_AutoYUY2_1(src_Y,src_Un,src_Vn,(uint8_t *)dst,dst_w>>2);
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;	
		}
	}
}



void AutoYUY2::Convert_Interlaced_YUY2_SSE(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstp=(uint8_t *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up, *src_Un, *src_Unn, *src_Unnn,*src_Upp;
	const uint8_t *src_V,*src_Vp, *src_Vn, *src_Vnn, *src_Vnnn,*src_Vpp;
	const int32_t w_8=(dst_w+15)>>4;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	bool _align=false;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Upp = srcUp - (2*src_pitchU);
	src_Un = srcUp + src_pitchU;
	src_Unn = srcUp + (2 * src_pitchU);
	src_Unnn = srcUp + (3 * src_pitchU);
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vpp = srcVp - (2*src_pitchV);
	src_Vn = srcVp + src_pitchV;
	src_Vnn = srcVp + (2 * src_pitchV);
	src_Vnnn = srcVp + (3 * src_pitchV);

	const ptrdiff_t src_pitchU_2=src_pitchU << 1;
	const ptrdiff_t src_pitchV_2=src_pitchV << 1;

	if ((((size_t)srcYp & 0x0F)==0) && ((abs(src_pitchY) & 0x0F)==0)) _align=true;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<4; y+=4)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_SSE2_1b(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1b(src_Y,src_Un,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_2b(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_3b(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_SSE2_1(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1(src_Y,src_Un,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_2(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_3(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}

	for (int32_t y = h_0; y < h_4; y+=4)
	{
		if (_align)
		{
			JPSDR_AutoYUY2_SSE2_3b(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_2b(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_2b(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_3b(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
		}
		else
		{
			JPSDR_AutoYUY2_SSE2_3(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_2(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_2(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_SSE2_3(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
		}
		dstp += dst_pitch;
		src_Y += src_pitchY;

		src_U += src_pitchU_2;
		src_Up += src_pitchU_2;
		src_Upp += src_pitchU_2;
		src_Un += src_pitchU_2;
		src_Unn += src_pitchU_2;
		src_Unnn += src_pitchU_2;
		src_V += src_pitchV_2;
		src_Vp += src_pitchV_2;
		src_Vpp += src_pitchV_2;
		src_Vn += src_pitchV_2;
		src_Vnn += src_pitchV_2;
		src_Vnnn += src_pitchV_2;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_4; y < h_Y_max; y+=4)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_SSE2_3b(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_2b(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1b(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1b(src_Y,src_Un,src_Vn,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_SSE2_3(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_2(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_SSE2_1(src_Y,src_Un,src_Vn,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;	
		}
	}
}



void AutoYUY2::Convert_Interlaced_YUY2_AVX(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	uint8_t *dstp=(uint8_t *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up, *src_Un, *src_Unn, *src_Unnn,*src_Upp;
	const uint8_t *src_V,*src_Vp, *src_Vn, *src_Vnn, *src_Vnnn,*src_Vpp;
	const int32_t w_8=(dst_w+15)>>4;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	bool _align=false;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Upp = srcUp - (2*src_pitchU);
	src_Un = srcUp + src_pitchU;
	src_Unn = srcUp + (2 * src_pitchU);
	src_Unnn = srcUp + (3 * src_pitchU);
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vpp = srcVp - (2*src_pitchV);
	src_Vn = srcVp + src_pitchV;
	src_Vnn = srcVp + (2 * src_pitchV);
	src_Vnnn = srcVp + (3 * src_pitchV);

	const ptrdiff_t src_pitchU_2=src_pitchU << 1;
	const ptrdiff_t src_pitchV_2=src_pitchV << 1;

	if ((((size_t)srcYp & 0x0F)==0) && ((abs(src_pitchY) & 0x0F)==0)) _align=true;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<4; y+=4)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_AVX_1b(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1b(src_Y,src_Un,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_2b(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_3b(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_AVX_1(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1(src_Y,src_Un,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_2(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_3(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}

	for (int32_t y = h_0; y < h_4; y+=4)
	{
		if (_align)
		{
			JPSDR_AutoYUY2_AVX_3b(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_2b(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_2b(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_3b(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
		}
		else
		{
			JPSDR_AutoYUY2_AVX_3(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_2(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_2(src_Y,src_Unn,srcUp,src_Vnn,srcVp,dstp,w_8);
			dstp += dst_pitch;
			src_Y += src_pitchY;
			JPSDR_AutoYUY2_AVX_3(src_Y,src_Un,src_Unnn,src_Vn,src_Vnnn,dstp,w_8);
		}
		dstp += dst_pitch;
		src_Y += src_pitchY;

		src_U += src_pitchU_2;
		src_Up += src_pitchU_2;
		src_Upp += src_pitchU_2;
		src_Un += src_pitchU_2;
		src_Unn += src_pitchU_2;
		src_Unnn += src_pitchU_2;
		src_V += src_pitchV_2;
		src_Vp += src_pitchV_2;
		src_Vpp += src_pitchV_2;
		src_Vn += src_pitchV_2;
		src_Vnn += src_pitchV_2;
		src_Vnnn += src_pitchV_2;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_4; y < h_Y_max; y+=4)
		{
			if (_align)
			{
				JPSDR_AutoYUY2_AVX_3b(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_2b(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1b(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1b(src_Y,src_Un,src_Vn,dstp,w_8);
			}
			else
			{
				JPSDR_AutoYUY2_AVX_3(src_Y,src_U,src_Upp,src_V,src_Vpp,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_2(src_Y,src_Up,src_Un,src_Vp,src_Vn,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1(src_Y,src_U,src_V,dstp,w_8);
				dstp += dst_pitch;
				src_Y += src_pitchY;
				JPSDR_AutoYUY2_AVX_1(src_Y,src_Un,src_Vn,dstp,w_8);
			}
			dstp += dst_pitch;
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;	
		}
	}
}



void AutoYUY2::Convert_Automatic_YUY2(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	YUYV *dst=(YUYV *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up,*src_Un,*src_Unn,*src_Unnn,*src_Upp;
	const uint8_t *src_V,*src_Vp,*src_Vn,*src_Vnn,*src_Vnnn,*src_Vpp;
	uint8_t index_tab_0,index_tab_1,index_tab_2;
	const int32_t w_UV=dst_w>>2;
	const int16_t threshold_=threshold;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Upp = srcUp - (2*src_pitchU);
	src_Un = srcUp + src_pitchU;
	src_Unn = srcUp + (2 * src_pitchU);
	src_Unnn = srcUp + (3 * src_pitchU);
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vpp = srcVp - (2*src_pitchV);
	src_Vn = srcVp + src_pitchV;
	src_Vnn = srcVp + (2 * src_pitchV);
	src_Vnnn = srcVp + (3 * src_pitchV);
	
	const ptrdiff_t src_pitchU_2=src_pitchU << 1;
	const ptrdiff_t src_pitchV_2=src_pitchV << 1;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<4; y+=4)
		{
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_U[j];
					dst[j].v=src_V[j];
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_Un[j];
					dst[j].v=src_Vn[j];
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
					dst[j].v=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;
				bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
				bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

				for (int32_t j=0; j<w_UV; j++)
				{	
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;

					if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
						(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
						(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
						((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
						itabu0[j]=true;
					else itabu0[j]=false;
					if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
						((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
						itabu1[j]=true;
					else itabu1[j]=false;

					dst[j].u=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);

					if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
						((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
						itabv0[j]=true;
					else itabv0[j]=false;
					if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
						((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
						itabv1[j]=true;
					else itabv1[j]=false;
				
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}
	else
	{
		bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
		bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

		for (int32_t j=0; j<w_UV; j++)
		{	
			if (((abs((int16_t)src_Upp[j]-(int16_t)src_Up[j])>=threshold_) &&
				(abs((int16_t)src_U[j]-(int16_t)src_Up[j])>=threshold_) &&
				(((src_Upp[j]>src_Up[j]) && (src_U[j]>src_Up[j])) ||
				((src_Upp[j]<src_Up[j]) && (src_U[j]<src_Up[j])))))
				itabu0[j]=true;
			else itabu0[j]=false;
			if (((abs((int16_t)src_Up[j]-(int16_t)src_U[j])>=threshold_) &&
				(abs((int16_t)src_Un[j]-(int16_t)src_U[j])>=threshold_) &&
				(((src_Up[j]>src_U[j]) && (src_Un[j]>src_U[j])) ||
				((src_Up[j]<src_U[j]) && (src_Un[j]<src_U[j])))))
				itabu1[j]=true;
			else itabu1[j]=false;

			if (((abs((int16_t)src_Vpp[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(abs((int16_t)src_V[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(((src_Vpp[j]>src_Vp[j]) && (src_V[j]>src_Vp[j])) ||
				((src_Vpp[j]<src_Vp[j]) && (src_V[j]<src_Vp[j])))))
				itabv0[j]=true;
			else itabv0[j]=false;
			if (((abs((int16_t)src_Vp[j]-(int16_t)src_V[j])>=threshold_) &&
				(abs((int16_t)src_Vn[j]-(int16_t)src_V[j])>=threshold_) &&
				(((src_Vp[j]>src_V[j]) && (src_Vn[j]>src_V[j])) ||
				((src_Vp[j]<src_V[j]) && (src_Vn[j]<src_V[j])))))
				itabv1[j]=true;
			else itabv1[j]=false;
		}
	}

	
	index_tab_0=0;
	index_tab_1=1;
	index_tab_2=2;

	for (int32_t y = h_0; y < h_4; y+=4)
	{
		{
			int32_t i=0;
			const bool *itabu0=interlaced_tab_U[thread_num][index_tab_0],*itabu1=interlaced_tab_U[thread_num][index_tab_1];
			const bool *itabv0=interlaced_tab_V[thread_num][index_tab_0],*itabv1=interlaced_tab_V[thread_num][index_tab_1];

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				// Upsample as needed.
				if ((itabu0[j]) || (itabu1[j]))
				{
					dst[j].u=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4) >> 3);
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2) >> 2);
				}

				// Upsample as needed.
				if ((itabv0[j]) || (itabv1[j]))
				{
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4) >> 3);
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2) >> 2);
				}
			}
		}

		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1];
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1];
			bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];
			bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for (int32_t j=0; j<w_UV; j++)
			{			
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
					(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
					(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
					((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
					itabu2[j]=true;
				else itabu2[j]=false;			

				// Upsample as needed.
				if ((itabu2[j]) || (itabu1[j]))
				{
					dst[j].u=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4) >> 3);
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2) >> 2);
				}

				if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
					((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
					itabv2[j]=true;
				else itabv2[j]=false;			

				// Upsample as needed.
				if ((itabv2[j]) || (itabv1[j]))
				{
					dst[j].v=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4) >> 3);
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2) >> 2);
				}
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;
		
		{
			int32_t i=0;
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1],*itabu2=interlaced_tab_U[thread_num][index_tab_2];
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1],*itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				// Upsample as needed.
				if ((itabu1[j]) || (itabu2[j]))
				{
					dst[j].u=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4) >> 3);
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_U[j]+2) >> 2);
				}

				// Upsample as needed.
				if ((itabv1[j]) || (itabv2[j]))
				{
					dst[j].v=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4) >> 3);
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_V[j]+2) >> 2);
				}
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;
			const bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];
			const bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];
			bool *itabu0=interlaced_tab_U[thread_num][index_tab_0];
			bool *itabv0=interlaced_tab_V[thread_num][index_tab_0];

			for (int32_t j=0; j<w_UV; j++)
			{			
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
					((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
					itabu0[j]=true;
				else itabu0[j]=false;

				// Upsample as needed.
				if ((itabu0[j]) || (itabu2[j]))
				{
					dst[j].u=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4) >> 3);
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_Unn[j]+2) >> 2);
				}

				if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
					((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
					itabv0[j]=true;
				else itabv0[j]=false;

				// Upsample as needed.
				if ((itabv0[j]) || (itabv2[j]))
				{
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4) >> 3);
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_Vnn[j]+2) >> 2);
				}
			}
		}
		
		index_tab_0=(index_tab_0+2)%Interlaced_Tab_Size;
		index_tab_1=(index_tab_1+2)%Interlaced_Tab_Size;
		index_tab_2=(index_tab_2+2)%Interlaced_Tab_Size;
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		src_U += src_pitchU_2;
		src_Up += src_pitchU_2;
		src_Upp += src_pitchU_2;
		src_Un += src_pitchU_2;
		src_Unn += src_pitchU_2;
		src_Unnn += src_pitchU_2;
		src_V += src_pitchV_2;
		src_Vp += src_pitchV_2;
		src_Vpp += src_pitchV_2;
		src_Vn += src_pitchV_2;
		src_Vnn += src_pitchV_2;
		src_Vnnn += src_pitchV_2;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_4; y < h_Y_max; y+=4)
		{
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4) >> 3);
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4) >> 3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4) >> 3);
					dst[j].v=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4) >> 3);
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;
		
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_U[j];
					dst[j].v=src_V[j];
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_Un[j];
					dst[j].v=src_Vn[j];
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}
}



void AutoYUY2::Convert_Test_YUY2(uint8_t thread_num)
{
	const MT_Data_Info_AutoYUY2 mt_data_inf=MT_Data[thread_num];

	const uint8_t *srcYp=(const uint8_t *)mt_data_inf.src1;
	const uint8_t *srcUp=(const uint8_t *)mt_data_inf.src2;
	const uint8_t *srcVp=(const uint8_t *)mt_data_inf.src3;
	YUYV *dst=(YUYV *)mt_data_inf.dst1;
	const int32_t dst_w=mt_data_inf.dst_Y_w;
	const int32_t h_Y_min=mt_data_inf.src_Y_h_min;
	const int32_t h_Y_max=mt_data_inf.src_Y_h_max;
	const int src_pitchY=mt_data_inf.src_pitch1;
	const int src_pitchU=mt_data_inf.src_pitch2;
	const int src_pitchV=mt_data_inf.src_pitch3;
	const int dst_pitch=mt_data_inf.dst_pitch1;

	const uint8_t *src_Y;
	const uint8_t *src_U,*src_Up,*src_Un,*src_Unn,*src_Unnn,*src_Upp;
	const uint8_t *src_V,*src_Vp,*src_Vn,*src_Vnn,*src_Vnnn,*src_Vpp;
	uint8_t index_tab_0,index_tab_1,index_tab_2;
	const int32_t w_UV=dst_w>>2;
	const int16_t threshold_=threshold;
	const uint16_t *lookup=lookup_Upscale;
	const int32_t h_4 = mt_data_inf.bottom ? h_Y_max-4:h_Y_max;
	const int32_t h_0 = mt_data_inf.top ? 4:h_Y_min;

	src_Y = srcYp;
	src_U = srcUp;
	src_Up = srcUp - src_pitchU;
	src_Upp = srcUp - (2*src_pitchU);
	src_Un = srcUp + src_pitchU;
	src_Unn = srcUp + (2 * src_pitchU);
	src_Unnn = srcUp + (3 * src_pitchU);
	src_V = srcVp;
	src_Vp = srcVp - src_pitchV;
	src_Vpp = srcVp - (2*src_pitchV);
	src_Vn = srcVp + src_pitchV;
	src_Vnn = srcVp + (2 * src_pitchV);
	src_Vnnn = srcVp + (3 * src_pitchV);
	
	const ptrdiff_t src_pitchU_2=src_pitchU << 1;
	const ptrdiff_t src_pitchV_2=src_pitchV << 1;

	if (mt_data_inf.top)
	{
		for (int32_t y=0; y<4; y+=4)
		{
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_U[j];
					dst[j].v=src_V[j];
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_Un[j];
					dst[j].v=src_Vn[j];
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_Unn[j]]+lookup[(uint16_t)src_U[j]+256]+4)>>3);
					dst[j].v=(uint8_t)((lookup[src_Vnn[j]]+lookup[(uint16_t)src_V[j]+256]+4)>>3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;
				bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
				bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

				for (int32_t j=0; j<w_UV; j++)
				{	
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;

					if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
						(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
						(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
						((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
						itabu0[j]=true;
					else itabu0[j]=false;
					if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
						(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
						((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
						itabu1[j]=true;
					else itabu1[j]=false;

					dst[j].u=(uint8_t)((lookup[(uint16_t)src_Un[j]+512]+(uint16_t)src_Unnn[j]+4)>>3);

					if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
						(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
						((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
						itabv0[j]=true;
					else itabv0[j]=false;
					if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
						(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
						((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
						itabv1[j]=true;
					else itabv1[j]=false;
				
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_Vn[j]+512]+(uint16_t)src_Vnnn[j]+4)>>3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}
	else
	{
		bool *itabu0=interlaced_tab_U[thread_num][0],*itabu1=interlaced_tab_U[thread_num][1];
		bool *itabv0=interlaced_tab_V[thread_num][0],*itabv1=interlaced_tab_V[thread_num][1];

		for (int32_t j=0; j<w_UV; j++)
		{	
			if (((abs((int16_t)src_Upp[j]-(int16_t)src_Up[j])>=threshold_) &&
				(abs((int16_t)src_U[j]-(int16_t)src_Up[j])>=threshold_) &&
				(((src_Upp[j]>src_Up[j]) && (src_U[j]>src_Up[j])) ||
				((src_Upp[j]<src_Up[j]) && (src_U[j]<src_Up[j])))))
				itabu0[j]=true;
			else itabu0[j]=false;
			if (((abs((int16_t)src_Up[j]-(int16_t)src_U[j])>=threshold_) &&
				(abs((int16_t)src_Un[j]-(int16_t)src_U[j])>=threshold_) &&
				(((src_Up[j]>src_U[j]) && (src_Un[j]>src_U[j])) ||
				((src_Up[j]<src_U[j]) && (src_Un[j]<src_U[j])))))
				itabu1[j]=true;
			else itabu1[j]=false;

			if (((abs((int16_t)src_Vpp[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(abs((int16_t)src_V[j]-(int16_t)src_Vp[j])>=threshold_) &&
				(((src_Vpp[j]>src_Vp[j]) && (src_V[j]>src_Vp[j])) ||
				((src_Vpp[j]<src_Vp[j]) && (src_V[j]<src_Vp[j])))))
				itabv0[j]=true;
			else itabv0[j]=false;
			if (((abs((int16_t)src_Vp[j]-(int16_t)src_V[j])>=threshold_) &&
				(abs((int16_t)src_Vn[j]-(int16_t)src_V[j])>=threshold_) &&
				(((src_Vp[j]>src_V[j]) && (src_Vn[j]>src_V[j])) ||
				((src_Vp[j]<src_V[j]) && (src_Vn[j]<src_V[j])))))
				itabv1[j]=true;
			else itabv1[j]=false;
		}
	}

	
	index_tab_0=0;
	index_tab_1=1;
	index_tab_2=2;

	for (int32_t y = h_0; y < h_4; y+=4)
	{
		{
			int32_t i=0;
			const bool *itabu0=interlaced_tab_U[thread_num][index_tab_0],*itabu1=interlaced_tab_U[thread_num][index_tab_1];
			const bool *itabv0=interlaced_tab_V[thread_num][index_tab_0],*itabv1=interlaced_tab_V[thread_num][index_tab_1];

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				// Upsample as needed.
				if ((itabu0[j]) || (itabu1[j]))
				{
					dst[j].u=239;
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Up[j]+2) >> 2);
				}

				// Upsample as needed.
				if ((itabv0[j]) || (itabv1[j]))
				{
					dst[j].v=239;
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vp[j]+2) >> 2);
				}
			}
		}

		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1];
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1];
			bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];
			bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for (int32_t j=0; j<w_UV; j++)
			{			
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				if (((abs((int16_t)src_U[j]-(int16_t)src_Un[j])>=threshold_) &&
					(abs((int16_t)src_Unn[j]-(int16_t)src_Un[j])>=threshold_) &&
					(((src_U[j]>src_Un[j]) && (src_Unn[j]>src_Un[j])) ||
					((src_U[j]<src_Un[j]) && (src_Unn[j]<src_Un[j])))))
					itabu2[j]=true;
				else itabu2[j]=false;			

				// Upsample as needed.
				if ((itabu2[j]) || (itabu1[j]))
				{
					dst[j].u=239;
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_U[j]]+(uint16_t)src_Un[j]+2) >> 2);
				}

				if (((abs((int16_t)src_V[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(abs((int16_t)src_Vnn[j]-(int16_t)src_Vn[j])>=threshold_) &&
					(((src_V[j]>src_Vn[j]) && (src_Vnn[j]>src_Vn[j])) ||
					((src_V[j]<src_Vn[j]) && (src_Vnn[j]<src_Vn[j])))))
					itabv2[j]=true;
				else itabv2[j]=false;			

				// Upsample as needed.
				if ((itabv2[j]) || (itabv1[j]))
				{
					dst[j].v=239;
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_V[j]]+(uint16_t)src_Vn[j]+2) >> 2);
				}
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;
		
		{
			int32_t i=0;
			const bool *itabu1=interlaced_tab_U[thread_num][index_tab_1],*itabu2=interlaced_tab_U[thread_num][index_tab_2];
			const bool *itabv1=interlaced_tab_V[thread_num][index_tab_1],*itabv2=interlaced_tab_V[thread_num][index_tab_2];

			for (int32_t j=0; j<w_UV; j++)
			{
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				// Upsample as needed.
				if ((itabu1[j]) || (itabu2[j]))
				{
					dst[j].u=239;
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_U[j]+2) >> 2);
				}

				// Upsample as needed.
				if ((itabv1[j]) || (itabv2[j]))
				{
					dst[j].v=239;
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_V[j]+2) >> 2);
				}
			}
		}
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		{
			int32_t i=0;
			const bool *itabu2=interlaced_tab_U[thread_num][index_tab_2];
			const bool *itabv2=interlaced_tab_V[thread_num][index_tab_2];
			bool *itabu0=interlaced_tab_U[thread_num][index_tab_0];
			bool *itabv0=interlaced_tab_V[thread_num][index_tab_0];

			for (int32_t j=0; j<w_UV; j++)
			{			
				dst[j].y1=src_Y[i];
				dst[j].y2=src_Y[i+1];
				i+=2;

				if (((abs((int16_t)src_Un[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(abs((int16_t)src_Unnn[j]-(int16_t)src_Unn[j])>=threshold_) &&
					(((src_Un[j]>src_Unn[j]) && (src_Unnn[j]>src_Unn[j])) ||
					((src_Un[j]<src_Unn[j]) && (src_Unnn[j]<src_Unn[j])))))
					itabu0[j]=true;
				else itabu0[j]=false;

				// Upsample as needed.
				if ((itabu0[j]) || (itabu2[j]))
				{
					dst[j].u=239;
				}
				else
				{
					dst[j].u=(uint8_t)((lookup[src_Un[j]]+(uint16_t)src_Unn[j]+2) >> 2);
				}

				if (((abs((int16_t)src_Vn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(abs((int16_t)src_Vnnn[j]-(int16_t)src_Vnn[j])>=threshold_) &&
					(((src_Vn[j]>src_Vnn[j]) && (src_Vnnn[j]>src_Vnn[j])) ||
					((src_Vn[j]<src_Vnn[j]) && (src_Vnnn[j]<src_Vnn[j])))))
					itabv0[j]=true;
				else itabv0[j]=false;

				// Upsample as needed.
				if ((itabv0[j]) || (itabv2[j]))
				{
					dst[j].v=239;
				}
				else
				{
					dst[j].v=(uint8_t)((lookup[src_Vn[j]]+(uint16_t)src_Vnn[j]+2) >> 2);
				}
			}
		}
		
		index_tab_0=(index_tab_0+2)%Interlaced_Tab_Size;
		index_tab_1=(index_tab_1+2)%Interlaced_Tab_Size;
		index_tab_2=(index_tab_2+2)%Interlaced_Tab_Size;
		
		dst=(YUYV *)((uint8_t *)dst+dst_pitch);
		src_Y += src_pitchY;

		src_U += src_pitchU_2;
		src_Up += src_pitchU_2;
		src_Upp += src_pitchU_2;
		src_Un += src_pitchU_2;
		src_Unn += src_pitchU_2;
		src_Unnn += src_pitchU_2;
		src_V += src_pitchV_2;
		src_Vp += src_pitchV_2;
		src_Vpp += src_pitchV_2;
		src_Vn += src_pitchV_2;
		src_Vnn += src_pitchV_2;
		src_Vnnn += src_pitchV_2;	
	}

	if (mt_data_inf.bottom)
	{
		for (int32_t y = h_4; y < h_Y_max; y+=4)
		{
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[(uint16_t)src_U[j]+512]+(uint16_t)src_Upp[j]+4) >> 3);
					dst[j].v=(uint8_t)((lookup[(uint16_t)src_V[j]+512]+(uint16_t)src_Vpp[j]+4) >> 3);
				}
			}

			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=(uint8_t)((lookup[src_Up[j]]+lookup[(uint16_t)src_Un[j]+256]+4) >> 3);
					dst[j].v=(uint8_t)((lookup[src_Vp[j]]+lookup[(uint16_t)src_Vn[j]+256]+4) >> 3);
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;
		
			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_U[j];
					dst[j].v=src_V[j];
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			{
				int32_t i=0;

				for (int32_t j=0; j<w_UV; j++)
				{
					dst[j].y1=src_Y[i];
					dst[j].y2=src_Y[i+1];
					i+=2;
					dst[j].u=src_Un[j];
					dst[j].v=src_Vn[j];
				}
			}
		
			dst=(YUYV *)((uint8_t *)dst+dst_pitch);
			src_Y += src_pitchY;

			src_U += src_pitchU_2;
			src_Up += src_pitchU_2;
			src_Upp += src_pitchU_2;
			src_Un += src_pitchU_2;
			src_Unn += src_pitchU_2;
			src_Unnn += src_pitchU_2;
			src_V += src_pitchV_2;
			src_Vp += src_pitchV_2;
			src_Vpp += src_pitchV_2;
			src_Vn += src_pitchV_2;
			src_Vnn += src_pitchV_2;
			src_Vnnn += src_pitchV_2;
		}
	}
}



void AutoYUY2::StaticThreadpool(void *ptr)
{
	const Public_MT_Data_Thread *data=(const Public_MT_Data_Thread *)ptr;
	AutoYUY2 *ptrClass=(AutoYUY2 *)data->pClass;
	
	switch(data->f_process)
	{
		case 1 : ptrClass->Convert_Progressive_YUY2(data->thread_Id);
			break;
		case 2 : ptrClass->Convert_Progressive_YUY2_SSE(data->thread_Id);
			break;
		case 3 : ptrClass->Convert_Interlaced_YUY2(data->thread_Id);
			break;
		case 4 : ptrClass->Convert_Interlaced_YUY2_SSE(data->thread_Id);
			break;
		case 5 : ptrClass->Convert_Automatic_YUY2(data->thread_Id);
			break;
		case 6 : ptrClass->Convert_Test_YUY2(data->thread_Id);
			break;
		case 7 : ptrClass->Convert_Progressive_YV16(data->thread_Id);
			break;
		case 8 : ptrClass->Convert_Progressive_YV16_SSE(data->thread_Id);
			break;
		case 9 : ptrClass->Convert_Interlaced_YV16(data->thread_Id);
			break;
		case 10 : ptrClass->Convert_Interlaced_YV16_SSE(data->thread_Id);
			break;
		case 11 : ptrClass->Convert_Automatic_YV16(data->thread_Id);
			break;
		case 12 : ptrClass->Convert_Test_YV16(data->thread_Id);
			break;
		case 13 : ptrClass->Convert_Progressive_YUY2_AVX(data->thread_Id);
			break;
		case 14 : ptrClass->Convert_Interlaced_YUY2_AVX(data->thread_Id);
			break;
		case 15 : ptrClass->Convert_Progressive_YV16_AVX(data->thread_Id);
			break;
		case 16 : ptrClass->Convert_Interlaced_YV16_AVX(data->thread_Id);
			break;
		default : ;
	}
}


PVideoFrame __stdcall AutoYUY2::GetFrame(int n, IScriptEnvironment* env) 
{
	PVideoFrame src = child->GetFrame(n,env);
	PVideoFrame dst = env->NewVideoFrame(vi,64);
	uint8_t *dst_Y,*dst_U,*dst_V;
	const uint8_t *srcYp = src->GetReadPtr(PLANAR_Y);
	const uint8_t *srcUp = src->GetReadPtr(PLANAR_U);
	const uint8_t *srcVp = src->GetReadPtr(PLANAR_V);
	int dst_h,dst_w;
	int dst_pitch_Y,dst_pitchU,dst_pitchV;
	int src_pitch_Y,src_pitchU,src_pitchV;
	uint8_t f_proc;

	switch(output)
	{
		case 0 : 
			dst_Y = dst->GetWritePtr();
			dst_h = dst->GetHeight();
			dst_w = dst->GetRowSize();
			dst_pitch_Y = dst->GetPitch();
			dst_U=NULL;
			dst_V=NULL;
			dst_pitchU=0;
			dst_pitchV=0;
			break;
		case 1 :
			dst_h = dst->GetHeight(PLANAR_Y);
			dst_w = dst->GetRowSize(PLANAR_Y);
			dst_Y = dst->GetWritePtr(PLANAR_Y);
			dst_U=dst->GetWritePtr(PLANAR_U);
			dst_V=dst->GetWritePtr(PLANAR_V);
			dst_pitch_Y = dst->GetPitch(PLANAR_Y);
			dst_pitchU=dst->GetPitch(PLANAR_U);
			dst_pitchV=dst->GetPitch(PLANAR_V);
			break;
		default : break;
	}
		
	src_pitch_Y = src->GetPitch(PLANAR_Y);
	src_pitchU = src->GetPitch(PLANAR_U);
	src_pitchV = src->GetPitch(PLANAR_V);

	if (threads_number>1)
	{
		if (!poolInterface->RequestThreadPool(UserId,threads_number,MT_Thread,-1,false))
			env->ThrowError("AutoYUY2: Error with the TheadPool while requesting threadpool !");
	}

	for(uint8_t i=0; i<threads_number; i++)
	{
		MT_Data[i].src1=(void *)(srcYp+(MT_Data[i].src_Y_h_min*src_pitch_Y));
		MT_Data[i].src2=(void *)(srcUp+(MT_Data[i].src_UV_h_min*src_pitchU));
		MT_Data[i].src3=(void *)(srcVp+(MT_Data[i].src_UV_h_min*src_pitchV));
		MT_Data[i].src_pitch1=src_pitch_Y;
		MT_Data[i].src_pitch2=src_pitchU;
		MT_Data[i].src_pitch3=src_pitchV;
		MT_Data[i].dst1=(void *)(dst_Y+(MT_Data[i].dst_Y_h_min*dst_pitch_Y));
		MT_Data[i].dst2=(void *)(dst_U+(MT_Data[i].dst_UV_h_min*dst_pitchU));
		MT_Data[i].dst3=(void *)(dst_V+(MT_Data[i].dst_UV_h_min*dst_pitchV));
		MT_Data[i].dst_pitch1=dst_pitch_Y;
		MT_Data[i].dst_pitch2=dst_pitchU;
		MT_Data[i].dst_pitch3=dst_pitchV;
	}

	switch(output)
	{
		case 0 :
			switch(mode)
			{
				case -1 : f_proc=5; break;
				case 0 :
					if ((dst_w & 0x7)==0)
					{
						if (AVX_Enable) f_proc=13;
						else
						{
							if (SSE2_Enable) f_proc=2;
							else f_proc=1;
						}
					}
					else f_proc=1;
					break;
				case 1 :
					if ((dst_w & 0x7)==0)
					{
						if (AVX_Enable) f_proc=14;
						else
						{
							if (SSE2_Enable) f_proc=4;
							else f_proc=3;
						}
					}
					else f_proc=3;
					break;
				case 2 : f_proc=6; break;
				default : f_proc=0; break;
			}
			break;
		case 1 :
			switch(mode)
			{
				case -1 : f_proc=11; break;
				case 0 :
					if ((dst_w & 0x7)==0)
					{
						if (AVX_Enable) f_proc=15;
						else
						{
							if (SSE2_Enable) f_proc=8;
							else f_proc=7;
						}
					}
					else f_proc=7;
					break;
				case 1 :
					if ((dst_w & 0x7)==0)
					{
						if (AVX_Enable) f_proc=16;
						else
						{
							if (SSE2_Enable) f_proc=10;
							else f_proc=9;
						}
					}
					else f_proc=9;
					break;
				case 2 : f_proc=12; break;
				default : f_proc=0; break;
			}
			break;
		default : f_proc=0; break;
	}

	if (threads_number>1)
	{
		for(uint8_t i=0; i<threads_number; i++)
			MT_Thread[i].f_process=f_proc;
		if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

		for(uint8_t i=0; i<threads_number; i++)
			MT_Thread[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep);
	}
	else
	{
		switch(f_proc)
		{
			case 1 : Convert_Progressive_YUY2(0); break;
			case 2 : Convert_Progressive_YUY2_SSE(0); break;
			case 3 : Convert_Interlaced_YUY2(0); break;
			case 4 : Convert_Interlaced_YUY2_SSE(0); break;
			case 5 : Convert_Automatic_YUY2(0); break;
			case 6 : Convert_Test_YUY2(0); break;
			case 7 : Convert_Progressive_YV16(0); break;
			case 8 : Convert_Progressive_YV16_SSE(0); break;
			case 9 : Convert_Interlaced_YV16(0); break;
			case 10 : Convert_Interlaced_YV16_SSE(0); break;
			case 11 : Convert_Automatic_YV16(0); break;
			case 12 : Convert_Test_YV16(0); break;
			case 13 : Convert_Progressive_YUY2_AVX(0); break;
			case 14 : Convert_Interlaced_YUY2_AVX(0); break;
			case 15 : Convert_Progressive_YV16_AVX(0); break;
			case 16 : Convert_Interlaced_YV16_AVX(0); break;
			default : break;
		}
	}

	return dst;
}
