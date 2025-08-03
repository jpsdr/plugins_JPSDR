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
; AVX512F,AVX512BW,AVX512DQ,AVX512VL

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
	vbroadcastss xmm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss xmm6,dword ptr[esi]
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
	vpmovzxbw zmm2,YMMWORD ptr[esi]
	vpmovzxbw zmm4,YMMWORD ptr[esi+ebx]

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

	vmovdqa XMMWORD ptr[edi],xmm0
	vmovdqa XMMWORD ptr[edi+16],xmm3
		
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
	vmovdqa64 zmm2,ZMMWORD ptr[esi]
	vmovdqa64 zmm4,ZMMWORD ptr[esi+ebx]

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
	
	vmovdqa64 ZMMWORD ptr[edi],zmm0
		
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
	vmovdqa64 zmm2,ZMMWORD ptr[esi]
	vmovdqa64 zmm4,ZMMWORD ptr[esi+ebx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[edx+eax]

	vpaddw zmm2,zmm2,zmm5
	vpaddw zmm3,zmm3,zmm5
	
	vpmaddwd zmm2,zmm2,zmm4
	vpmaddwd zmm3,zmm3,zmm4

	vpaddd zmm0,zmm0,zmm2
	vpaddd zmm1,zmm1,zmm3

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

	vmovdqa64 ZMMWORD ptr[edi],zmm0

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
	vbroadcastss zmm4,dword ptr[edx+eax]
	vbroadcastss zmm5,dword ptr[edx+eax+4]
	
	vfmadd231ps zmm0,zmm4,ZMMWORD ptr[esi]
	vfmadd231ps zmm1,zmm5,ZMMWORD ptr[esi+ebx]

	add esi,ebx
	add eax,8
	add esi,ebx
	loop Resize_V_AVX512_Planar_32bits_loop_2
	
	vaddps zmm0,zmm0,zmm1
	
	vmovaps ZMMWORD ptr[edi],zmm0
		
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


Resize_H_AVX512_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_64:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_H_AVX512_Planar_8bits_ASM
	
	local sizehx2 : dword

	push ebx
	push edi
	push esi		

	mov esi,valmin
	vbroadcastss xmm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss xmm6,dword ptr[esi]
	mov esi,rounder
	vbroadcastss xmm7,dword ptr[esi]
	
	mov eax,sizeh
	shr eax,1
	jz Resize_H_AVX512_Planar_8bits_1
	mov sizehx2,eax

Resize_H_AVX512_Planar_8bits_loop_1:
	mov eax,32
	mov edx,64
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_64 ;kernel_size_64 = (kernel_size + 31) >> 5
	mov ebx,src_pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_8bits_loop_2:
	vmovdqa64 zmm4,ZMMWORD ptr[edi]			;coeff

	vpmovzxbw zmm2,YMMWORD ptr[esi]			;src
	vpmovzxbw zmm3,YMMWORD ptr[esi+ebx]		;src+src_pitch
	
	vpmaddwd zmm2,zmm2,zmm4
	vpmaddwd zmm3,zmm3,zmm4

	vpaddd zmm0,zmm0,zmm2
	vpaddd zmm1,zmm1,zmm3

	add esi,eax
	add edi,edx
	loop Resize_H_AVX512_Planar_8bits_loop_2

	vextracti32x8 ymm2,zmm0,1
	vextracti32x8 ymm3,zmm1,1
	
	vphaddd ymm0,ymm0,ymm2
	vphaddd ymm1,ymm1,ymm3
	
	vextracti128 xmm2,ymm0,1
	vextracti128 xmm3,ymm1,1

	vphaddd xmm0,xmm0,xmm2
	vphaddd xmm1,xmm1,xmm3
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	sal ebx,1

	vphaddd xmm0,xmm0,xmm1

	add src,ebx

	vpaddd xmm0,xmm0,xmm7

	mov ebx,dst_pitch

	vpsrad xmm0,xmm0,14

	mov edi,dst

	vpackusdw xmm0,xmm0,xmm0

	mov edx,ebx

	vpackuswb xmm0,xmm0,xmm0
	
	sal edx,1
	
	vpmaxub xmm0,xmm0,xmm5

	add dst,edx

	vpminub xmm0,xmm0,xmm6
	
	vpextrb eax,xmm0,0
	vpextrb edx,xmm0,2
	
	mov byte ptr[edi],al
	mov byte ptr[edi+ebx],dl
	
	dec sizehx2
	jnz Resize_H_AVX512_Planar_8bits_loop_1
	
Resize_H_AVX512_Planar_8bits_1:	
	test sizeh,1
	jz short Resize_H_AVX512_Planar_8bits_end

	mov eax,32
	mov edx,64
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_64

	vpxord zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_8bits_loop_3:
	vmovdqa64 zmm4,ZMMWORD ptr[edi]	;coeff

	vpmovzxbw zmm2,YMMWORD ptr[esi]	;src
	
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm0,zmm0,zmm2

	add esi,eax
	add edi,edx
	loop Resize_H_AVX512_Planar_8bits_loop_3

	vextracti32x8 ymm2,zmm0,1
	
	vphaddd ymm0,ymm0,ymm2
	
	vextracti128 xmm2,ymm0,1

	vphaddd xmm0,xmm0,xmm2
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0

	vpaddd xmm0,xmm0,xmm7 ;rounder

	vpsrad xmm0,xmm0,14 ;FPScale16bits = 14

	mov edi,dst

	vpackusdw xmm0,xmm0,xmm0
	vpackuswb xmm0,xmm0,xmm0

	vpmaxub xmm0,xmm0,xmm5
	vpminub xmm0,xmm0,xmm6
	
	vpextrb eax,xmm0,0
	
	mov byte ptr[edi],al

Resize_H_AVX512_Planar_8bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX512_Planar_8bits_ASM endp


Resize_H_AVX512_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_64:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_H_AVX512_Planar_10to14bits_ASM
	
	local sizehx2 : dword

	push ebx
	push edi
	push esi		

	mov esi,valmin
	vbroadcastss xmm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss xmm6,dword ptr[esi]
	mov esi,rounder
	vbroadcastss xmm7,dword ptr[esi]
	
	mov eax,sizeh
	shr eax,1
	jz Resize_H_AVX512_Planar_10to14bits_1
	mov sizehx2,eax

Resize_H_AVX512_Planar_10to14bits_loop_1:
	mov eax,64
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_64 ;kernel_size_64 = (kernel_size + 31) >> 5
	mov ebx,src_pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_10to14bits_loop_2:
	vmovdqa64 zmm4,ZMMWORD ptr[edi] ;coeff

	vpmaddwd zmm2,zmm4,ZMMWORD ptr[esi]		;src
	vpmaddwd zmm3,zmm4,ZMMWORD ptr[esi+ebx]	;src+src_pitch

	vpaddd zmm0,zmm0,zmm2
	vpaddd zmm1,zmm1,zmm3

	add esi,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_10to14bits_loop_2

	vextracti32x8 ymm2,zmm0,1
	vextracti32x8 ymm3,zmm1,1
	
	vphaddd ymm0,ymm0,ymm2
	vphaddd ymm1,ymm1,ymm3
	
	vextracti128 xmm2,ymm0,1
	vextracti128 xmm3,ymm1,1

	vphaddd xmm0,xmm0,xmm2
	vphaddd xmm1,xmm1,xmm3
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	sal ebx,1

	vphaddd xmm0,xmm0,xmm1

	add src,ebx

	vpaddd xmm0,xmm0,xmm7 ;rounder

	mov ebx,dst_pitch

	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov edi,dst

	vpackusdw xmm0,xmm0,xmm0

	mov edx,ebx

	vpmaxuw xmm0,xmm0,xmm5
	
	sal edx,1
	
	vpminuw xmm0,xmm0,xmm6

	add dst,edx
	
	vpextrw eax,xmm0,0
	vpextrw edx,xmm0,2
	
	mov word ptr[edi],ax
	mov word ptr[edi+ebx],dx
	
	dec sizehx2
	jnz Resize_H_AVX512_Planar_10to14bits_loop_1

Resize_H_AVX512_Planar_10to14bits_1:
	test sizeh,1
	jz short Resize_H_AVX512_Planar_10to14bits_end

	mov eax,64
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_64

	vpxord zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_10to14bits_loop_3:
	vmovdqa64 zmm4,ZMMWORD ptr[edi]
	vpmaddwd zmm2,zmm4,ZMMWORD ptr[esi]

	vpaddd zmm0,zmm0,zmm2

	add esi,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_10to14bits_loop_3

	vextracti32x8 ymm2,zmm0,1
	
	vphaddd ymm0,ymm0,ymm2
	
	vextracti128 xmm2,ymm0,1

	vphaddd xmm0,xmm0,xmm2
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0

	vpaddd xmm0,xmm0,xmm7

	vpsrad xmm0,xmm0,13

	mov edi,dst

	vpackusdw xmm0,xmm0,xmm0

	vpmaxuw xmm0,xmm0,xmm5
	vpminuw xmm0,xmm0,xmm6
	
	vpextrw eax,xmm0,0
	
	mov word ptr[edi],ax

Resize_H_AVX512_Planar_10to14bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX512_Planar_10to14bits_ASM endp


Resize_H_AVX512_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_64:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword,
	shifttosigned:dword,shiftfromsigned:dword

    public Resize_H_AVX512_Planar_16bits_ASM
	
	local sizehx2 : dword
	local XValMin : XMMWORD
	local XValMax : XMMWORD
	local XShiftfromSigned : XMMWORD

	push ebx
	push edi
	push esi		

	mov esi,shifttosigned
	vbroadcastss zmm5,dword ptr[esi]
	mov esi,valmin
	vbroadcastss xmm6,dword ptr[esi]
	vmovdqu XValMin,xmm6
	mov esi,valmax
	vbroadcastss xmm6,dword ptr[esi]
	vmovdqu XValMax,xmm6
	mov esi,shiftfromsigned
	vbroadcastss xmm6,dword ptr[esi]
	vmovdqu XShiftfromSigned,xmm6
	mov esi,rounder
	vbroadcastss xmm7,dword ptr[esi]
	
	mov eax,sizeh
	shr eax,1
	jz Resize_H_AVX512_Planar_16bits_1
	mov sizehx2,eax

Resize_H_AVX512_Planar_16bits_loop_1:
	mov eax,64
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_64 ;kernel_size_64 = (kernel_size + 31) >> 5
	mov ebx,src_pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_16bits_loop_2:
	vmovdqa64 zmm4,ZMMWORD ptr[edi]			;coeff

	vpaddw zmm2,zmm5,ZMMWORD ptr[esi]		;src
	vpaddw zmm3,zmm5,ZMMWORD ptr[esi+ebx]	;src+src_pitch

	vpmaddwd zmm2,zmm2,zmm4
	vpmaddwd zmm3,zmm3,zmm4

	vpaddd zmm0,zmm0,zmm2
	vpaddd zmm1,zmm1,zmm3

	add esi,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_16bits_loop_2

	vextracti32x8 ymm2,zmm0,1
	vextracti32x8 ymm3,zmm1,1
	
	vphaddd ymm0,ymm0,ymm2
	vphaddd ymm1,ymm1,ymm3
	
	vextracti128 xmm2,ymm0,1
	vextracti128 xmm3,ymm1,1

	vphaddd xmm0,xmm0,xmm2
	vphaddd xmm1,xmm1,xmm3
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	sal ebx,1

	vphaddd xmm0,xmm0,xmm1

	add src,ebx

	vpaddd xmm0,xmm0,xmm6 ;ShiftfromSigned

	mov ebx,dst_pitch

	vpaddd xmm0,xmm0,xmm7 ;rounder

	vmovdqu xmm6,XValMin

	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov edi,dst

	vpackusdw xmm0,xmm0,xmm0

	mov edx,ebx

	vpmaxuw xmm0,xmm0,xmm6
	
	sal edx,1
	
	vmovdqu xmm6,XValMax

	add dst,edx
	
	vpminuw xmm0,xmm0,xmm6

	vmovdqu xmm6,XShiftfromSigned

	vpextrw eax,xmm0,0
	vpextrw edx,xmm0,2
	
	mov word ptr[edi],ax
	mov word ptr[edi+ebx],dx
	
	dec sizehx2
	jnz Resize_H_AVX512_Planar_16bits_loop_1

Resize_H_AVX512_Planar_16bits_1:
	test sizeh,1
	jz Resize_H_AVX512_Planar_16bits_end

	mov eax,64
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_64

	vpxord zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_16bits_loop_3:
	vmovdqa64 zmm4,ZMMWORD ptr[edi]	;coeff

	vpaddw zmm2,zmm5,ZMMWORD ptr[esi]

	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm0,zmm0,zmm2

	add esi,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_16bits_loop_3

	vextracti32x8 ymm2,zmm0,1
	
	vphaddd ymm0,ymm0,ymm2
	
	vextracti128 xmm2,ymm0,1

	vphaddd xmm0,xmm0,xmm2
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0

	vpaddd xmm0,xmm0,xmm6 ;ShiftfromSigned

	mov edi,dst

	vpaddd xmm0,xmm0,xmm7 ;rounder

	vmovdqu xmm6,XValMin

	vpsrad xmm0,xmm0,13

	vpackusdw xmm0,xmm0,xmm0

	vpmaxuw xmm0,xmm0,xmm6
	vmovdqu xmm6,XValMax
	vpminuw xmm0,xmm0,xmm6

	vpextrw eax,xmm0,0
	
	mov word ptr[edi],ax

Resize_H_AVX512_Planar_16bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX512_Planar_16bits_ASM endp


Resize_H_AVX512_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_64:dword,sizeh:dword

    public Resize_H_AVX512_Planar_32bits_ASM
	
	local sizehx4 : dword

	push ebx
	push edi
	push esi		
	
	mov eax,sizeh
	shr eax,2
	jz Resize_H_AVX512_Planar_32bits_1
	mov sizehx4,eax

Resize_H_AVX512_Planar_32bits_loop_1:
	mov ebx,src_pitch
	mov eax,64
	mov ecx,kernel_size_64 ;kernel_size_64 = (kernel_size + 15) >> 4
	mov edi,coeff
	mov esi,src
	mov edx,esi
	add edx,ebx

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1
	vxorps zmm2,zmm2,zmm2
	vxorps zmm3,zmm3,zmm3

Resize_H_AVX512_Planar_32bits_loop_2:
	vmovaps zmm4,ZMMWORD ptr[edi] ;coeff
	
	vfmadd231ps zmm0,zmm4,ZMMWORD ptr[esi]			;src
	vfmadd231ps zmm1,zmm4,ZMMWORD ptr[edx]			;src+src_pitch
	vfmadd231ps zmm2,zmm4,ZMMWORD ptr[esi+2*ebx]	;src+2*src_pitch
	vfmadd231ps zmm3,zmm4,ZMMWORD ptr[edx+2*ebx]	;src+3*src_pitch

	add esi,eax
	add edx,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_32bits_loop_2

	vextractf32x8 ymm5,zmm0,1
	vextractf32x8 ymm6,zmm1,1
	vhaddps ymm0,ymm0,ymm5
	vhaddps ymm1,ymm1,ymm6
	vextractf32x8 ymm5,zmm2,1
	vextractf32x8 ymm6,zmm3,1
	vhaddps ymm2,ymm2,ymm5
	vhaddps ymm3,ymm3,ymm6

	vextractf128 xmm5,ymm0,1
	vextractf128 xmm6,ymm1,1
	vhaddps xmm0,xmm0,xmm5
	vhaddps xmm1,xmm1,xmm6
	vextractf128 xmm5,ymm2,1
	vextractf128 xmm6,ymm3,1
	vhaddps xmm2,xmm2,xmm5
	vhaddps xmm3,xmm3,xmm6

	vhaddps xmm0,xmm0,xmm0
	vhaddps xmm1,xmm1,xmm1
	vhaddps xmm2,xmm2,xmm2
	vhaddps xmm3,xmm3,xmm3

	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	mov edi,dst

	vhaddps xmm0,xmm0,xmm2

	mov edx,dst_pitch

	vpextrd eax,xmm0,0
	mov dword ptr[edi],eax
	sal ebx,2
	vpextrd eax,xmm0,1
	add edi,edx
	add src,ebx
	mov dword ptr[edi],eax
	mov ebx,edx
	add edi,edx
	vpextrd eax,xmm0,2
	sal ebx,2
	mov dword ptr[edi],eax
	add edi,edx
	vpextrd eax,xmm0,3
	add dst,ebx
	mov dword ptr[edi],eax
	
	dec sizehx4
	jnz Resize_H_AVX512_Planar_32bits_loop_1

Resize_H_AVX512_Planar_32bits_1:
	test sizeh,2
	jz Resize_H_AVX512_Planar_32bits_2

	mov ebx,src_pitch
	mov eax,64
	mov ecx,kernel_size_64 ;kernel_size_64 = (kernel_size + 15) >> 4
	mov edi,coeff
	mov esi,src

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_32bits_loop_3:
	vmovaps zmm4,ZMMWORD ptr[edi] ;coeff
	
	vfmadd231ps zmm0,zmm4,ZMMWORD ptr[esi]			;src
	vfmadd231ps zmm1,zmm4,ZMMWORD ptr[esi+ebx]			;src+src_pitch

	add esi,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_32bits_loop_3

	vextractf32x8 ymm5,zmm0,1
	vextractf32x8 ymm6,zmm1,1
	vhaddps ymm0,ymm0,ymm5
	vhaddps ymm1,ymm1,ymm6

	vextractf128 xmm5,ymm0,1
	vextractf128 xmm6,ymm1,1
	vhaddps xmm0,xmm0,xmm5
	vhaddps xmm1,xmm1,xmm6

	vhaddps xmm0,xmm0,xmm0
	vhaddps xmm1,xmm1,xmm1

	vhaddps xmm0,xmm0,xmm1

	mov edi,dst
	vpextrd eax,xmm0,0
	mov edx,dst_pitch
	mov dword ptr[edi],eax
	sal ebx,1
	vpextrd eax,xmm0,2
	add edi,edx
	add src,ebx
	mov dword ptr[edi],eax
	mov ebx,edx
	sal ebx,1
	add dst,ebx

Resize_H_AVX512_Planar_32bits_2:
	test sizeh,1
	jz short Resize_H_AVX512_Planar_32bits_end

	mov ebx,src_pitch
	mov eax,64
	mov ecx,kernel_size_64 ;kernel_size_64 = (kernel_size + 15) >> 4
	mov edi,coeff
	mov esi,src

	vxorps zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_32bits_loop_4:
	vmovaps zmm4,ZMMWORD ptr[edi] 			;coeff
	
	vfmadd231ps zmm0,zmm4,ZMMWORD ptr[esi]	;src

	add esi,eax
	add edi,eax
	loop Resize_H_AVX512_Planar_32bits_loop_4

	vextractf32x8 ymm5,zmm0,1
	vhaddps ymm0,ymm0,ymm5

	vextractf128 xmm5,ymm0,1
	vhaddps xmm0,xmm0,xmm5

	vhaddps xmm0,xmm0,xmm0

	vhaddps xmm0,xmm0,xmm0

	mov edi,dst
	vpextrd eax,xmm0,0
	mov dword ptr[edi],eax

Resize_H_AVX512_Planar_32bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX512_Planar_32bits_ASM endp

end
