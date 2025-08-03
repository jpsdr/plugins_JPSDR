;
;                    AVX2 ResampleMT for Avs+/Avisynth 2.6.x
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

Resize_V_AVX_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,src_pitch:dword,
	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_V_AVX_Planar_8bits_ASM

	push ebx
	push edi
	push esi		

	mov esi,valmin
	vbroadcastss xmm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss xmm6,dword ptr[esi]
	mov esi,rounder
	vbroadcastss xmm7,dword ptr[esi]

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX_Planar_8bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vmovdqa xmm0,xmm7
	vmovdqa xmm1,xmm7

Resize_V_AVX_Planar_8bits_loop_2:
	vpmovzxbw xmm2,qword ptr[esi]
	vpmovzxbw xmm4,qword ptr[esi+ebx]

	vpunpckhwd xmm3,xmm2,xmm4
	vpunpcklwd xmm2,xmm2,xmm4

	vbroadcastss xmm4,dword ptr[edx+eax]

	vpmaddwd xmm3,xmm3,xmm4
	vpmaddwd xmm2,xmm2,xmm4

	vpaddd xmm1,xmm1,xmm3
	vpaddd xmm0,xmm0,xmm2

	add esi,ebx
	add eax,4
	add esi,ebx
	loop Resize_V_AVX_Planar_8bits_loop_2

	vpsrad xmm0,xmm0,14 ;FPScale8bits = 14
	vpsrad xmm1,xmm1,14

	vpackusdw xmm0,xmm0,xmm1

	vpackuswb xmm0,xmm0,xmm0
	vpmaxub xmm0,xmm0,xmm5
	vpminub xmm0,xmm0,xmm6

	vmovq qword ptr[edi],xmm0
		
	add edi,8
	add src,8
	dec width8
	jnz Resize_V_AVX_Planar_8bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX_Planar_8bits_ASM endp


Resize_V_AVX_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,
	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_V_AVX_Planar_10to14bits_ASM

	push ebx
	push edi
	push esi		

	mov esi,valmin
	vbroadcastss xmm5,dword ptr[esi]
	mov esi,valmax
	vbroadcastss xmm6,dword ptr[esi]
	mov esi,rounder
	vbroadcastss xmm7,dword ptr[esi]

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX_Planar_10to14bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vmovdqa xmm0,xmm7
	vmovdqa xmm1,xmm7

Resize_V_AVX_Planar_10to14bits_loop_2:
	vmovdqa xmm2,XMMWORD ptr[esi]
	vmovdqa xmm4,XMMWORD ptr[esi+ebx]

	vpunpckhwd xmm3,xmm2,xmm4
	vpunpcklwd xmm2,xmm2,xmm4

	vbroadcastss xmm4,dword ptr[edx+eax]

	vpmaddwd xmm3,xmm3,xmm4
	vpmaddwd xmm2,xmm2,xmm4

	vpaddd xmm1,xmm1,xmm3
	vpaddd xmm0,xmm0,xmm2

	add esi,ebx
	add eax,4
	add esi,ebx
	loop Resize_V_AVX_Planar_10to14bits_loop_2

	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	vpsrad xmm1,xmm1,13

	vpackusdw xmm0,xmm0,xmm1

	vpmaxuw xmm0,xmm0,xmm5
	vpminuw xmm0,xmm0,xmm6
	
	vmovdqa XMMWORD ptr[edi],xmm0
		
	add edi,16
	add src,16
	dec width8
	jnz Resize_V_AVX_Planar_10to14bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX_Planar_10to14bits_ASM endp


Resize_V_AVX_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,
	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword,
	shifttosigned:dword,shiftfromsigned:dword

    public Resize_V_AVX_Planar_16bits_ASM
	
	local XValMin : XMMWORD
	local XValMax : XMMWORD
	local XShiftfromSigned : XMMWORD
	
	push ebx
	push edi
	push esi		

	mov esi,shifttosigned
	vbroadcastss xmm5,dword ptr[esi]
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

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX_Planar_16bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vmovdqa xmm0,xmm7
	vmovdqa xmm1,xmm7

Resize_V_AVX_Planar_16bits_loop_2:
	vmovdqa xmm2,XMMWORD ptr[esi]
	vmovdqa xmm4,XMMWORD ptr[esi+ebx]

	vpunpckhwd xmm3,xmm2,xmm4
	vpunpcklwd xmm2,xmm2,xmm4

	vbroadcastss xmm4,dword ptr[edx+eax]

	vpaddw xmm2,xmm2,xmm5
	vpaddw xmm3,xmm3,xmm5
	
	vpmaddwd xmm2,xmm2,xmm4
	vpmaddwd xmm3,xmm3,xmm4

	vpaddd xmm0,xmm0,xmm2
	vpaddd xmm1,xmm1,xmm3

	add esi,ebx
	add eax,4
	add esi,ebx
	loop Resize_V_AVX_Planar_16bits_loop_2

	vpaddw xmm0,xmm0,xmm6
	vpaddw xmm1,xmm1,xmm6

	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	vpsrad xmm1,xmm1,13

	vmovdqu xmm6,XValMin

	vpackusdw xmm0,xmm0,xmm1

	vpmaxuw xmm0,xmm0,xmm6

	vmovdqu xmm6,XValMax
	vpminuw xmm0,xmm0,xmm6

	vmovdqa XMMWORD ptr[edi],xmm0

	vmovdqu xmm6,XShiftfromSigned
		
	add edi,16
	add src,16
	dec width8
	jnz Resize_V_AVX_Planar_16bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX_Planar_16bits_ASM endp


Resize_V_AVX_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,
	src_pitch:dword,kernel_size_2:dword

    public Resize_V_AVX_Planar_32bits_ASM

	push ebx
	push edi
	push esi		

	mov edx,coeff
	mov ebx,src_pitch

	mov edi,dst

Resize_V_AVX_Planar_32bits_loop_1:
	mov ecx,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov esi,src
	xor eax,eax

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1

Resize_V_AVX_Planar_32bits_loop_2:
	vbroadcastss ymm4,dword ptr[edx+eax]
	vbroadcastss ymm5,dword ptr[edx+eax+4]
	
	vmulps ymm2,ymm4,YMMWORD ptr[esi]
	vmulps ymm3,ymm5,YMMWORD ptr[esi+ebx]
	
	vaddps ymm0,ymm0,ymm2
	vaddps ymm1,ymm1,ymm3

	add esi,ebx
	add eax,8
	add esi,ebx
	loop Resize_V_AVX_Planar_32bits_loop_2
	
	vaddps ymm0,ymm0,ymm1
	
	vmovaps YMMWORD ptr[edi],ymm0
		
	add edi,32
	add src,32
	dec width8
	jnz short Resize_V_AVX_Planar_32bits_loop_1
	
	vzeroupper

	pop esi
	pop edi
	pop ebx

	ret

Resize_V_AVX_Planar_32bits_ASM endp


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


Resize_H_AVX_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_16:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_H_AVX_Planar_8bits_ASM
	
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
	jz Resize_H_AVX_Planar_8bits_1
	mov sizehx2,eax

Resize_H_AVX_Planar_8bits_loop_1:
	mov eax,8
	mov edx,16
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_16 ;kernel_size_16 = (kernel_size + 7) >> 3
	mov ebx,src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1

Resize_H_AVX_Planar_8bits_loop_2:
	vmovdqa xmm4,XMMWORD ptr[edi]			;coeff

	vpmovzxbw xmm2,qword ptr[esi]			;src
	vpmovzxbw xmm3,qword ptr[esi+ebx]		;src+src_pitch
	
	vpmaddwd xmm2,xmm2,xmm4
	vpmaddwd xmm3,xmm3,xmm4

	vpaddd xmm0,xmm0,xmm2
	vpaddd xmm1,xmm1,xmm3

	add esi,eax
	add edi,edx
	loop Resize_H_AVX_Planar_8bits_loop_2
	
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
	jnz Resize_H_AVX_Planar_8bits_loop_1
	
Resize_H_AVX_Planar_8bits_1:	
	test sizeh,1
	jz short Resize_H_AVX_Planar_8bits_end

	mov eax,8
	mov edx,16
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_16

	vpxor xmm0,xmm0,xmm0

Resize_H_AVX_Planar_8bits_loop_3:
	vmovdqa xmm4,XMMWORD ptr[edi]	;coeff

	vpmovzxbw xmm2,qword ptr[esi]	;src
	
	vpmaddwd xmm2,xmm2,xmm4

	vpaddd xmm0,xmm0,xmm2

	add esi,eax
	add edi,edx
	loop Resize_H_AVX_Planar_8bits_loop_3
	
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

Resize_H_AVX_Planar_8bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX_Planar_8bits_ASM endp


Resize_H_AVX_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_16:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

    public Resize_H_AVX_Planar_10to14bits_ASM
	
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
	jz Resize_H_AVX_Planar_10to14bits_1
	mov sizehx2,eax

Resize_H_AVX_Planar_10to14bits_loop_1:
	mov eax,16
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_16 ;kernel_size_16 = (kernel_size + 7) >> 3
	mov ebx,src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1

Resize_H_AVX_Planar_10to14bits_loop_2:
	vmovdqa xmm4,XMMWORD ptr[edi] ;coeff

	vpmaddwd xmm2,xmm4,XMMWORD ptr[esi]		;src
	vpmaddwd xmm3,xmm4,XMMWORD ptr[esi+ebx]	;src+src_pitch

	vpaddd xmm0,xmm0,xmm2
	vpaddd xmm1,xmm1,xmm3

	add esi,eax
	add edi,eax
	loop Resize_H_AVX_Planar_10to14bits_loop_2
	
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
	jnz Resize_H_AVX_Planar_10to14bits_loop_1

Resize_H_AVX_Planar_10to14bits_1:
	test sizeh,1
	jz short Resize_H_AVX_Planar_10to14bits_end

	mov eax,16
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_16

	vpxor xmm0,xmm0,xmm0

Resize_H_AVX_Planar_10to14bits_loop_3:
	vmovdqa xmm4,XMMWORD ptr[edi]
	vpmaddwd xmm2,xmm4,XMMWORD ptr[esi]

	vpaddd xmm0,xmm0,xmm2

	add esi,eax
	add edi,eax
	loop Resize_H_AVX_Planar_10to14bits_loop_3
	
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

Resize_H_AVX_Planar_10to14bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX_Planar_10to14bits_ASM endp


Resize_H_AVX_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_16:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword,
	shifttosigned:dword,shiftfromsigned:dword

    public Resize_H_AVX_Planar_16bits_ASM
	
	local sizehx2 : dword
	local XValMin : XMMWORD
	local XValMax : XMMWORD
	local XShiftfromSigned : XMMWORD

	push ebx
	push edi
	push esi		

	mov esi,shifttosigned
	vbroadcastss xmm5,dword ptr[esi]
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
	jz Resize_H_AVX_Planar_16bits_1
	mov sizehx2,eax

Resize_H_AVX_Planar_16bits_loop_1:
	mov eax,16
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_16 ;kernel_size_16 = (kernel_size + 7) >> 3
	mov ebx,src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1

Resize_H_AVX_Planar_16bits_loop_2:
	vmovdqa xmm4,XMMWORD ptr[edi]			;coeff

	vpaddw xmm2,xmm5,XMMWORD ptr[esi]		;src
	vpaddw xmm3,xmm5,XMMWORD ptr[esi+ebx]	;src+src_pitch

	vpmaddwd xmm2,xmm2,xmm4
	vpmaddwd xmm3,xmm3,xmm4

	vpaddd xmm0,xmm0,xmm2
	vpaddd xmm1,xmm1,xmm3

	add esi,eax
	add edi,eax
	loop Resize_H_AVX_Planar_16bits_loop_2
	
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
	jnz Resize_H_AVX_Planar_16bits_loop_1

Resize_H_AVX_Planar_16bits_1:
	test sizeh,1
	jz Resize_H_AVX_Planar_16bits_end

	mov eax,16
	mov edi,coeff
	mov esi,src
	mov ecx,kernel_size_16

	vpxor xmm0,xmm0,xmm0

Resize_H_AVX_Planar_16bits_loop_3:
	vmovdqa xmm4,XMMWORD ptr[edi]	;coeff

	vpaddw xmm2,xmm5,XMMWORD ptr[esi]

	vpmaddwd xmm2,xmm2,xmm4

	vpaddd xmm0,xmm0,xmm2

	add esi,eax
	add edi,eax
	loop Resize_H_AVX_Planar_16bits_loop_3
	
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

Resize_H_AVX_Planar_16bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX_Planar_16bits_ASM endp


Resize_H_AVX_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
	kernel_size_32:dword,sizeh:dword

    public Resize_H_AVX_Planar_32bits_ASM
	
	local sizehx4 : dword

	push ebx
	push edi
	push esi		
	
	mov eax,sizeh
	shr eax,2
	jz Resize_H_AVX_Planar_32bits_1
	mov sizehx4,eax

Resize_H_AVX_Planar_32bits_loop_1:
	mov ebx,src_pitch
	mov eax,32
	mov ecx,kernel_size_32 ;kernel_size_32 = (kernel_size + 7) >> 3
	mov edi,coeff
	mov esi,src
	mov edx,esi
	add edx,ebx

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1
	vxorps ymm2,ymm2,ymm2
	vxorps ymm3,ymm3,ymm3

Resize_H_AVX_Planar_32bits_loop_2:
	vmovaps ymm4,YMMWORD ptr[edi] ;coeff
	
	vmulps ymm5,ymm4,YMMWORD ptr[esi]		;src
	vmulps ymm6,ymm4,YMMWORD ptr[edx]		;src + pitch
	vaddps ymm0,ymm0,ymm5
	vaddps ymm1,ymm1,ymm6

	vmulps ymm5,ymm4,YMMWORD ptr[esi+2*ebx]	;src+2*src_pitch
	vmulps ymm6,ymm4,YMMWORD ptr[edx+2*ebx]	;src+3*src_pitch
	vaddps ymm2,ymm2,ymm5
	vaddps ymm3,ymm3,ymm6

	add esi,eax
	add edx,eax
	add edi,eax
	loop Resize_H_AVX_Planar_32bits_loop_2

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
	jnz Resize_H_AVX_Planar_32bits_loop_1

Resize_H_AVX_Planar_32bits_1:
	test sizeh,2
	jz Resize_H_AVX_Planar_32bits_2

	mov ebx,src_pitch
	mov eax,32
	mov ecx,kernel_size_32 ;kernel_size_32 = (kernel_size + 7) >> 3
	mov edi,coeff
	mov esi,src

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1

Resize_H_AVX_Planar_32bits_loop_3:
	vmovaps ymm4,YMMWORD ptr[edi] ;coeff

	vmulps ymm5,ymm4,YMMWORD ptr[esi]		;src
	vmulps ymm6,ymm4,YMMWORD ptr[esi+ebx]	;src + pitch
	vaddps ymm0,ymm0,ymm5
	vaddps ymm1,ymm1,ymm6

	add esi,eax
	add edi,eax
	loop Resize_H_AVX_Planar_32bits_loop_3

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

Resize_H_AVX_Planar_32bits_2:
	test sizeh,1
	jz short Resize_H_AVX_Planar_32bits_end

	mov ebx,src_pitch
	mov eax,32
	mov ecx,kernel_size_32 ;kernel_size_32 = (kernel_size + 7) >> 3
	mov edi,coeff
	mov esi,src

	vxorps ymm0,ymm0,ymm0

Resize_H_AVX_Planar_32bits_loop_4:
	vmovaps ymm4,YMMWORD ptr[edi] 			;coeff

	vmulps ymm5,ymm4,YMMWORD ptr[esi]		;src
	vaddps ymm0,ymm0,ymm5

	add esi,eax
	add edi,eax
	loop Resize_H_AVX_Planar_32bits_loop_4

	vextractf128 xmm5,ymm0,1
	vhaddps xmm0,xmm0,xmm5

	vhaddps xmm0,xmm0,xmm0

	vhaddps xmm0,xmm0,xmm0

	mov edi,dst
	vpextrd eax,xmm0,0
	mov dword ptr[edi],eax

Resize_H_AVX_Planar_32bits_end:
	vzeroupper

	pop esi
	pop edi
	pop ebx
	
	ret

Resize_H_AVX_Planar_32bits_ASM endp

end
