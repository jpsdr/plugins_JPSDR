;
;  AutoYUY2()
;
;  Adaptive YV12 upsampling. Progressive picture areas are upsampled
;  progressively and interlaced areas are upsampled interlaced.
;  Copyright (C) 2005 Donald A. Graft
;  ASM part made by JPSDR
;	
;  AutoYUY2 is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2, or (at your option)
;  any later version.
;   
;  AutoYUY2 is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
;  GNU General Public License for more details.
;   
;  You should have received a copy of the GNU General Public License
;  along with GNU Make; see the file COPYING.  If not, write to
;  the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA. 
;

.586
.xmm
.model flat,c

.code


JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_8_AVX2 proc src1:dword,src2:dword,dst:dword,w32:dword

	public JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_8_AVX2
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	mov ecx,w32	
	mov ebx,32
	
Convert_Planar420_to_Planar422_x3x1_8_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[esi+eax]
	vmovdqa ymm1,YMMWORD ptr[edx+eax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgb ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[edi+eax],ymm2
	add eax,ebx
	loop Convert_Planar420_to_Planar422_x3x1_8_AVX2_1
	
	vzeroupper
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_8_AVX2 endp


JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_16_AVX2 proc src1:dword,src2:dword,dst:dword,w16:dword

	public JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_16_AVX2
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax	
	mov ecx,w16
	mov ebx,32
	
Convert_Planar420_to_Planar422_x3x1_16_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[esi+eax]
	vmovdqa ymm1,YMMWORD ptr[edx+eax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgw ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[edi+eax],ymm2
	add eax,ebx
	loop Convert_Planar420_to_Planar422_x3x1_16_AVX2_1
	
	vzeroupper
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_16_AVX2 endp


JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_8_AVX2 proc src1:dword,src2:dword,dst:dword,w32:dword

	public JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_8_AVX2
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	mov ecx,w32	
	mov ebx,32
	
Convert_Planar420_to_Planar422_x3x5_8_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[esi+eax]
	vmovdqa ymm1,YMMWORD ptr[edx+eax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm2,ymm2,ymm1
	vpavgb ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgb ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[edi+eax],ymm2
	add eax,ebx
	loop Convert_Planar420_to_Planar422_x3x5_8_AVX2_1
	
	vzeroupper
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_8_AVX2 endp


JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_16_AVX2 proc src1:dword,src2:dword,dst:dword,w16:dword

	public JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_16_AVX2
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	mov ecx,w16
	mov ebx,32
	
Convert_Planar420_to_Planar422_x3x5_16_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[esi+eax]
	vmovdqa ymm1,YMMWORD ptr[edx+eax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm2,ymm2,ymm1
	vpavgw ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgw ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[edi+eax],ymm2
	add eax,ebx
	loop Convert_Planar420_to_Planar422_x3x5_16_AVX2_1
	
	vzeroupper
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_16_AVX2 endp


JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_8_AVX2 proc src1:dword,src2:dword,dst:dword,w32:dword

	public JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_8_AVX2
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	mov ecx,w32	
	mov ebx,32
	
Convert_Planar420_to_Planar422_x7x1_8_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[esi+eax]
	vmovdqa ymm1,YMMWORD ptr[edx+eax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm1,ymm1,ymm2
	vpavgb ymm1,ymm1,ymm2
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm1,ymm1,ymm0
	vmovdqa YMMWORD ptr[edi+eax],ymm1
	add eax,ebx
	loop Convert_Planar420_to_Planar422_x7x1_8_AVX2_1
	
	vzeroupper
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_8_AVX2 endp


JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_16_AVX2 proc src1:dword,src2:dword,dst:dword,w16:dword

	public JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_16_AVX2
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	mov ecx,w16
	mov ebx,32
	
Convert_Planar420_to_Planar422_x7x1_16_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[esi+eax]
	vmovdqa ymm1,YMMWORD ptr[edx+eax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm1,ymm1,ymm2
	vpavgw ymm1,ymm1,ymm2
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm1,ymm1,ymm0
	vmovdqa YMMWORD ptr[edi+eax],ymm1
	add eax,ebx
	loop Convert_Planar420_to_Planar422_x7x1_16_AVX2_1
	
	vzeroupper
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_16_AVX2 endp


end





