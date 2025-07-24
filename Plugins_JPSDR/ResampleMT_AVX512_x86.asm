;
;                    AVX512 ResampleMT for Avs+/Avisynth 2.6.x
;
;   Copyright (C) 2025 JPSDR from PINTERF AVX2
;
;   This program is free software; you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation; either version 2 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;
;   You should have received a copy of the GNU General Public License
;   along with this program; if not, write to the Free Software
;   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
;
;

.xmm
.model flat,c

.code

Resize_V_AVX512_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,width32:dword,src_pitch:dword,
	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_V_AVX512_Planar_8bits_ASM

	push ebx
	push edi
	push esi		

	mov esi,valmin
	vbroadcastss ymm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss ymm6,dword ptr[esi]
	mov esi,rounder
	vbroadcastss zmm7,dword ptr[esi]

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX512_Planar_8bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vmovdqa64 zmm0,zmm7
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_8bits_loop_2:
	vpmovzxbw zmm2,YMMWORD ptr [esi]
	vpmovzxbw zmm4,YMMWORD ptr [esi+ebx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[edx+eax]

	vpmaddwd zmm3,zmm3,zmm4
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add esi,ebx
	add eax,4
	add esi,ebx
	loop Resize_V_AVX512_Planar_8bits_loop_2

	vpsrad zmm0,zmm0,14 ;FPScale8bits = 14
	vpsrad zmm1,zmm1,14

	vpackusdw zmm0,zmm0,zmm1

	vextracti32x8 ymm3,zmm0,1
	vextracti128 xmm1,ymm0,1
	vextracti128 xmm4,ymm3,1
	vpackuswb xmm0,xmm0,xmm1
	vpackuswb xmm3,xmm3,xmm4
	vpmaxub xmm0,xmm0,xmm5
	vpmaxub xmm3,xmm3,xmm5
	vpminub xmm0,xmm0,xmm6
	vpminub xmm3,xmm3,xmm6
	vmovdqa XMMWORD ptr [edi],xmm0
	vmovdqa XMMWORD ptr [edi+16],xmm3
		
	add edi,32
	add src,32
	dec width32
	jnz Resize_V_AVX512_Planar_8bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX512_Planar_8bits_ASM endp


Resize_V_AVX512_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,width32:dword,
	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_V_AVX512_Planar_10to14bits_ASM

	push ebx
	push edi
	push esi		

	mov esi,valmin
	vbroadcastss zmm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss zmm6,dword ptr[esi]
	mov esi,rounder
	vbroadcastss zmm7,dword ptr[esi]

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX512_Planar_10to14bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vmovdqa64 zmm0,zmm7
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_10to14bits_loop_2:
	vmovdqa64 zmm2,ZMMWORD ptr [esi]
	vmovdqa64 zmm4,ZMMWORD ptr [esi+ebx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[edx+eax]

	vpmaddwd zmm3,zmm3,zmm4
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add esi,ebx
	add eax,4
	add esi,ebx
	loop Resize_V_AVX512_Planar_10to14bits_loop_2

	vpsrad zmm0,zmm0,13 ;FPScale16bits = 13
	vpsrad zmm1,zmm1,13

	vpackusdw zmm0,zmm0,zmm1

	vpmaxuw zmm0,zmm0,zmm5
	vpminuw zmm0,zmm0,zmm6
	
	vmovdqa64 ZMMWORD ptr [edi],zmm0
		
	add edi,64
	add src,64
	dec width32
	jnz Resize_V_AVX512_Planar_10to14bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX512_Planar_10to14bits_ASM endp


Resize_V_AVX512_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,width32:dword,
	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword,
	shifttosigned:dword,shiftfromsigned:dword

    public Resize_V_AVX512_Planar_16bits_ASM
	
	local ZValMin : ZMMWORD
	local ZValMax : ZMMWORD
	local ZShiftfromSigned : ZMMWORD
	
	push ebx
	push edi
	push esi		

	mov esi,shifttosigned
	vbroadcastss zmm5,dword ptr[esi]
	mov esi,valmin
	vbroadcastss zmm6,dword ptr[esi]
	vmovdqu64 ZValMin,zmm6
	mov esi,valmax
	vbroadcastss zmm6,dword ptr[esi]
	vmovdqu64 ZValMax,zmm6
	mov esi,shiftfromsigned
	vbroadcastss zmm6,dword ptr[esi]
	vmovdqu64 ZShiftfromSigned,zmm6
	mov esi,rounder
	vbroadcastss zmm7,dword ptr[esi]

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX512_Planar_16bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vmovdqa64 zmm0,zmm7
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_16bits_loop_2:
	vmovdqa64 zmm2,ZMMWORD ptr [esi]
	vmovdqa64 zmm4,ZMMWORD ptr [esi+ebx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[edx+eax]

	vpaddw zmm3,zmm3,zmm5
	vpaddw zmm2,zmm2,zmm5
	
	vpmaddwd zmm3,zmm3,zmm4
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add esi,ebx
	add eax,4
	add esi,ebx
	loop Resize_V_AVX512_Planar_16bits_loop_2

	vpaddw zmm0,zmm0,zmm6
	vpaddw zmm1,zmm1,zmm6

	vpsrad zmm0,zmm0,13 ;FPScale16bits = 13
	vpsrad zmm1,zmm1,13

	vmovdqu64 zmm6,ZValMin

	vpackusdw zmm0,zmm0,zmm1

	vpmaxuw zmm0,zmm0,zmm6

	vmovdqu64 zmm6,ZValMax
	vpminuw zmm0,zmm0,zmm6

	vmovdqa64 ZMMWORD ptr [edi],zmm0

	vmovdqu64 zmm6,ZShiftfromSigned
		
	add edi,64
	add src,64
	dec width32
	jnz Resize_V_AVX512_Planar_16bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX512_Planar_16bits_ASM endp


Resize_V_AVX512_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,
	src_pitch:dword,kernel_size_2:dword

    public Resize_V_AVX512_Planar_32bits_ASM

	push ebx
	push edi
	push esi		

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX512_Planar_32bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1

Resize_V_AVX512_Planar_32bits_loop_2:
	vmovaps zmm2,ZMMWORD ptr [esi]
	vmovaps zmm3,ZMMWORD ptr [esi+ebx]

	vbroadcastss zmm4,dword ptr[edx+eax]
	vbroadcastss zmm5,dword ptr[edx+eax+4]
	
	vfmadd231ps zmm0,zmm2,zmm4
	vfmadd231ps zmm1,zmm3,zmm5

	add esi,ebx
	add eax,8
	add esi,ebx
	loop Resize_V_AVX512_Planar_32bits_loop_2
	
	vaddps zmm0,zmm0,zmm1
	
	vmovaps ZMMWORD ptr [edi],zmm0
		
	add edi,64
	add src,64
	dec width16
	jnz short Resize_V_AVX512_Planar_32bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX512_Planar_32bits_ASM endp

end
