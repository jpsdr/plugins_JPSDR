;
;                    nnedi3 for Avs+/Avisynth 2.6.x
;
;   Copyright (C) 2010-2011 Kevin Stone
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
;   Modified by JPSDR
;

.data

FLT_EPSILON  equ   1.192092896e-07

align 16

sign_bits_f_zero_l qword 7FFFFFFF00000000h,7FFFFFFF7FFFFFFFh
sign_bits_f qword 2 dup(7FFFFFFF7FFFFFFFh)
ones_f real4 4 dup(1.0)

flt_epsilon_sse real4 4 dup(FLT_EPSILON)

min_weight_sum real4 4 dup(1.0e-10)
five_f real4 4 dup(5.0)

sse_half real4 4 dup(0.5)

data segment align(64)

exp_hi real4 16 dup(80.0)
exp_lo real4 16 dup(-80.0)

; exp from:  A Fast, Compact Approximation of the Exponential Function (1998)
;            Nicol N. Schraudolph

e0_mult real4 16 dup(12102203.161561486)   ; (1.0/ln(2))*(2^23)
e0_bias real4 16 dup(1064866805.0)         ; (2^23)*127.0-486411.0

; exp from Loren Merritt

e1_scale real4 16 dup(1.4426950409)        ; 1/ln(2)
e1_bias real4 16 dup(12582912.0)           ; 3<<22
e1_c1 real4 16 dup(0.701277797)
e1_c2 real4 16 dup(0.237348593)
e1_c0 real4 16 dup(1.00035)

; exp from Intel Approximate Math (AM) Library

exp_rln2 real4 16 dup(1.442695041)
am_0p5 real4 16 dup(0.5)
epi32_1 sdword 16 dup(1)
exp_c2 real4 16 dup(1.428606820e-6)
exp_c1 real4 16 dup(6.931457520e-1)
exp_q0 real4 16 dup(3.001985051e-6)
exp_p0 real4 16 dup(1.261771931e-4)
epi32_0x7f sdword 16 dup(7Fh)
exp_q1 real4 16 dup(2.524483403e-3)
exp_p1 real4 16 dup(3.029944077e-2)
exp_q2 real4 16 dup(2.272655482e-1)
am_1 real4 16 dup(1.0)
exp_q3 real4 16 dup(2.0)

w_19 sword 32 dup(19)
w_3 sword 32 dup(3)
uw_16 word 32 dup(16)
ub_1 byte 64 dup(1)

d_19 sdword 16 dup(19)
d_3 sdword 16 dup(3)
ud_16 dword 16 dup(16)
uw_1 word 32 dup(1)

f_19 real4 16 dup(0.59375)
f_3 real4 16 dup(0.09375)

sign_bits_f_32 qword 8 dup(7FFFFFFF7FFFFFFFh)
ones_f_32 real4 16 dup(1.0)

.code


;computeNetwork0_i16_AVX512 proc inputf:dword,weightsf:dword,ptr_d:dword
; inputf = rcx
; weightsf = rdx
; ptr_d = r8

computeNetwork0_i16_AVX512 proc public frame

	.endprolog

		mov rax,1

		vmovdqa64 zmm17,ZMMWORD ptr [rcx]
		vpmaddwd zmm0,zmm17,ZMMWORD ptr [rdx]
		vpmaddwd zmm1,zmm17,ZMMWORD ptr [rdx+64]
		vpmaddwd zmm2,zmm17,ZMMWORD ptr [rdx+128]
		vpmaddwd zmm3,zmm17,ZMMWORD ptr [rdx+192]
		
		vmovdqa64 zmm17,ZMMWORD ptr [rcx+64]
		vpmaddwd zmm4,zmm17,ZMMWORD ptr [rdx+256]
		vpmaddwd zmm5,zmm17,ZMMWORD ptr [rdx+320]
		vpmaddwd zmm16,zmm17,ZMMWORD ptr [rdx+384]
		vpmaddwd zmm17,zmm17,ZMMWORD ptr [rdx+448]
		vpaddd zmm0,zmm0,zmm4
		vpaddd zmm1,zmm1,zmm5
		vpaddd zmm2,zmm2,zmm16
		vpaddd zmm3,zmm3,zmm17

		vmovdqa32 ymm18,ymm0
		vmovdqa32 ymm19,ymm1
		
		vextracti32x8 ymm16,zmm0,1
		vextracti32x8 ymm5,zmm1,1
		vextracti32x8 ymm17,zmm2,1
		vextracti32x8 ymm4,zmm3,1

		vpunpckhqdq ymm0,ymm16,ymm5
		vpunpckhqdq ymm1,ymm17,ymm4
		vpunpcklqdq ymm16,ymm16,ymm5
		vpunpcklqdq ymm17,ymm17,ymm4
		vpaddd ymm16,ymm16,ymm0
		vpaddd ymm17,ymm17,ymm1

		vmovdqa32 ymm0,ymm18
		vmovdqa32 ymm1,ymm19
		
		vpunpckhqdq ymm4,ymm0,ymm1
		vpunpckhqdq ymm5,ymm2,ymm3
		vpunpcklqdq ymm0,ymm0,ymm1
		vpunpcklqdq ymm2,ymm2,ymm3
		vpaddd ymm0,ymm0,ymm4
		vpaddd ymm2,ymm2,ymm5

		vpaddd ymm0,ymm0,ymm16
		vpaddd ymm2,ymm2,ymm17

		vextracti128 xmm4,ymm0,1
		vextracti128 xmm5,ymm2,1
		vpaddd xmm0,xmm0,xmm4
		vpaddd xmm2,xmm2,xmm5
		
		vshufps xmm16,xmm0,xmm2,221
		vshufps xmm0,xmm0,xmm2,136
		vpaddd xmm0,xmm0,xmm16
		vcvtdq2ps xmm0,xmm0
		vmulps xmm0,xmm0,XMMWORD ptr [rdx+512]
		vaddps xmm0,xmm0,XMMWORD ptr [rdx+528]
		vmovaps xmm1,xmm0
		vandps xmm0,xmm0,XMMWORD ptr sign_bits_f_zero_l
		vaddps xmm0,xmm0,XMMWORD ptr ones_f
		vrcpps xmm0,xmm0
		vmulps xmm0,xmm0,xmm1
		vpshufd xmm1,xmm0,0
		vpshufd xmm2,xmm0,85
		vpshufd xmm3,xmm0,170
		vpshufd xmm4,xmm0,255
		vmulps xmm1,xmm1,XMMWORD ptr [rdx+544]
		vmulps xmm2,xmm2,XMMWORD ptr [rdx+544+16]
		vmulps xmm3,xmm3,XMMWORD ptr [rdx+544+32]
		vmulps xmm4,xmm4,XMMWORD ptr [rdx+544+48]
		vaddps xmm1,xmm1,xmm2
		vaddps xmm3,xmm3,xmm4
		vaddps xmm1,xmm1,xmm3
		vaddps xmm1,xmm1,XMMWORD ptr [rdx+544+64]
		vmovaps xmm17,xmm1
		vandps xmm1,xmm1,XMMWORD ptr sign_bits_f
		vmovaps xmm3,xmm0
		vaddps xmm1,xmm1,XMMWORD ptr ones_f
		vrcpps xmm1,xmm1
		vmulps xmm17,xmm17,xmm1
		vpshufd xmm0,xmm0,0
		vpshufd xmm1,xmm3,85
		vpshufd xmm2,xmm3,170
		vpshufd xmm3,xmm3,255
		vmulps xmm0,xmm0,XMMWORD ptr [rdx+624]
		vmulps xmm1,xmm1,XMMWORD ptr [rdx+624+16]
		vmulps xmm2,xmm2,XMMWORD ptr [rdx+624+32]
		vmulps xmm3,xmm3,XMMWORD ptr [rdx+624+48]
		vpshufd xmm4,xmm17,0
		vpshufd xmm5,xmm17,85
		vpshufd xmm16,xmm17,170
		vpshufd xmm17,xmm17,255
		vmulps xmm4,xmm4,XMMWORD ptr [rdx+624+64]
		vmulps xmm5,xmm5,XMMWORD ptr [rdx+624+80]
		vmulps xmm16,xmm16,XMMWORD ptr [rdx+624+96]
		vmulps xmm17,xmm17,XMMWORD ptr [rdx+624+112]
		vaddps xmm0,xmm0,xmm1
		vaddps xmm2,xmm2,xmm3
		vaddps xmm4,xmm4,xmm5
		vaddps xmm16,xmm16,xmm17
		vaddps xmm0,xmm0,xmm2
		vaddps xmm4,xmm4,xmm16
		vaddps xmm0,xmm0,xmm4
		vaddps xmm0,xmm0,XMMWORD ptr [rdx+624+128]
		vmovhlps xmm1,xmm1,xmm0
		vmaxps xmm0,xmm0,xmm1
		vpshuflw xmm1,xmm0,14
		vcomiss xmm1,xmm0
		jbe short finish_2
		xor rax,rax
finish_2:
		mov BYTE PTR[r8],al

	vzeroupper

	ret

computeNetwork0_i16_AVX512 endp


;computeNetwork0new_AVX512 proc datai:dword,weights:dword,ptr_d:dword
; datai = rcx
; weights = rdx
; ptr_d = r8

computeNetwork0new_AVX512 proc public frame

	.endprolog

		vmovdqa64 zmm17,ZMMWORD ptr [rcx]
		vpmaddwd zmm0,zmm17,ZMMWORD ptr [rdx]
		vpmaddwd zmm1,zmm17,ZMMWORD ptr [rdx+64]
		vpmaddwd zmm2,zmm17,ZMMWORD ptr [rdx+128]
		vpmaddwd zmm3,zmm17,ZMMWORD ptr [rdx+192]
		
		vmovdqa64 zmm17,ZMMWORD ptr [rcx+64]
		vpmaddwd zmm4,zmm17,ZMMWORD ptr [rdx+256]
		vpmaddwd zmm5,zmm17,ZMMWORD ptr [rdx+320]
		vpmaddwd zmm16,zmm17,ZMMWORD ptr [rdx+384]
		vpmaddwd zmm17,zmm17,ZMMWORD ptr [rdx+448]
		vpaddd zmm0,zmm0,zmm4
		vpaddd zmm1,zmm1,zmm5
		vpaddd zmm2,zmm2,zmm16
		vpaddd zmm3,zmm3,zmm17
		
		vextracti32x8 ymm16,zmm0,1
		vextracti32x8 ymm5,zmm1,1
		vextracti32x8 ymm17,zmm2,1
		vextracti32x8 ymm4,zmm3,1

		vpunpckhqdq ymm18,ymm16,ymm5
		vpunpckhqdq ymm19,ymm17,ymm4
		vpunpcklqdq ymm16,ymm16,ymm5
		vpunpcklqdq ymm17,ymm17,ymm4
		vpaddd ymm16,ymm16,ymm18
		vpaddd ymm17,ymm17,ymm19

		vpunpckhqdq ymm4,ymm0,ymm1
		vpunpckhqdq ymm5,ymm2,ymm3
		vpunpcklqdq ymm0,ymm0,ymm1
		vpunpcklqdq ymm2,ymm2,ymm3
		vpaddd ymm0,ymm0,ymm4
		vpaddd ymm2,ymm2,ymm5

		vpaddd ymm0,ymm0,ymm16
		vpaddd ymm2,ymm2,ymm17

		vextracti128 xmm4,ymm0,1
		vextracti128 xmm5,ymm2,1
		vpaddd xmm0,xmm0,xmm4
		vpaddd xmm2,xmm2,xmm5
		
		vshufps xmm16,xmm0,xmm2,221
		vshufps xmm0,xmm0,xmm2,136	
		
		vpaddd xmm0,xmm0,xmm16
		vcvtdq2ps xmm0,xmm0		
		vmulps xmm0,xmm0,XMMWORD ptr [rdx+512]
		vaddps xmm0,xmm0,XMMWORD ptr [rdx+528]
		vmovaps xmm1,xmm0
		vandps xmm0,xmm0,XMMWORD ptr sign_bits_f
		vaddps xmm0,xmm0,XMMWORD ptr ones_f
		vrcpps xmm0,xmm0
		vmulps xmm0,xmm0,xmm1
		vpshufd xmm1,xmm0,0
		vpshufd xmm2,xmm0,85
		vpshufd xmm3,xmm0,170
		vpshufd xmm4,xmm0,255
		vmulps xmm1,xmm1,XMMWORD ptr [rdx+544]
		vmulps xmm2,xmm2,XMMWORD ptr [rdx+560]
		vmulps xmm3,xmm3,XMMWORD ptr [rdx+576]
		vmulps xmm4,xmm4,XMMWORD ptr [rdx+592]
		vpxor xmm0,xmm0,xmm0
		vaddps xmm1,xmm1,xmm2
		vaddps xmm3,xmm3,xmm4
		vaddps xmm1,xmm1,xmm3
		vaddps xmm1,xmm1,XMMWORD ptr [rdx+608]
		vcmpps xmm1,xmm1,xmm0,1
		vpackssdw xmm1,xmm1,xmm0
		vpacksswb xmm1,xmm1,xmm0
				
		vmovd eax,xmm1
		xor eax,0FFFFFFFFh
		and eax,001010101h
		mov [r8],eax

	vzeroupper

	ret

computeNetwork0new_AVX512 endp


; From FMA3
;computeNetwork0_AVX512 proc input:dword,weights:dword,ptr_d:dword
; input = rcx
; weights = rdx
; ptr_d = r8

computeNetwork0_AVX512 proc public frame

	sub rsp,40
	.allocstack 40
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
		mov rax,1

		vmovaps zmm4,ZMMWORD ptr [rcx]
		vmulps zmm0,zmm4,ZMMWORD ptr [rdx]
		vmulps zmm1,zmm4,ZMMWORD ptr [rdx+64]
		vmulps zmm2,zmm4,ZMMWORD ptr [rdx+128]
		vmulps zmm3,zmm4,ZMMWORD ptr [rdx+192]
		
		vmovaps zmm4,ZMMWORD ptr [rcx+64]
		vfmadd231ps zmm0,zmm4,ZMMWORD ptr [rdx+256]
		vfmadd231ps zmm1,zmm4,ZMMWORD ptr [rdx+320]
		vfmadd231ps zmm2,zmm4,ZMMWORD ptr [rdx+384]
		vfmadd231ps zmm3,zmm4,ZMMWORD ptr [rdx+448]
	
		vmovaps zmm4,ZMMWORD ptr [rcx+128]
		vfmadd231ps zmm0,zmm4,ZMMWORD ptr [rdx+512]
		vfmadd231ps zmm1,zmm4,ZMMWORD ptr [rdx+576]
		vfmadd231ps zmm2,zmm4,ZMMWORD ptr [rdx+640]
		vfmadd231ps zmm3,zmm4,ZMMWORD ptr [rdx+704]
		
		vextractf32x8 ymm6,zmm0,1
		vextractf32x8 ymm5,zmm1,1
		vextractf32x8 ymm7,zmm2,1
		vextractf32x8 ymm4,zmm3,1
		
		vhaddps ymm0,ymm0,ymm1
		vhaddps ymm2,ymm2,ymm3
		vhaddps ymm0,ymm0,ymm2

		vhaddps ymm6,ymm6,ymm5
		vhaddps ymm7,ymm7,ymm4
		vhaddps ymm6,ymm6,ymm7
		
		vaddps ymm0,ymm0,ymm6

		vextractf128 xmm4,ymm0,1
		vaddps xmm0,xmm0,xmm4
		
		vaddps xmm0,xmm0,XMMWORD ptr [rdx+768]
		
		vmovaps xmm1,xmm0
		vandps xmm0,xmm0,XMMWORD ptr sign_bits_f_zero_l
		vaddps xmm0,xmm0,XMMWORD ptr ones_f
		vrcpps xmm0,xmm0
		vmulps xmm0,xmm0,xmm1
		
		vpshufd xmm1,xmm0,0
		vpshufd xmm2,xmm0,85
		vpshufd xmm3,xmm0,170
		vpshufd xmm4,xmm0,255
		
		vmulps xmm1,xmm1,XMMWORD ptr [rdx+784]
		vfmadd231ps xmm1,xmm2,XMMWORD ptr [rdx+784+16]
		vmulps xmm3,xmm3,XMMWORD ptr [rdx+784+32]
		vfmadd231ps xmm3,xmm4,XMMWORD ptr [rdx+784+48]
		vaddps xmm1,xmm1,xmm3
		vaddps xmm1,xmm1,XMMWORD ptr [rdx+784+64]
		
		vmovaps xmm7,xmm1
		vandps xmm1,xmm1, XMMWORD ptr sign_bits_f
		vmovaps xmm3,xmm0
		vaddps xmm1,xmm1,XMMWORD ptr ones_f
		vrcpps xmm1,xmm1
		vmulps xmm7,xmm7,xmm1

		vpshufd xmm0,xmm0,0
		vpshufd xmm1,xmm3,85
		vpshufd xmm2,xmm3,170
		vpshufd xmm3,xmm3,255
		vmulps xmm0,xmm0,XMMWORD ptr [rdx+864]
		vfmadd231ps xmm0,xmm1,XMMWORD ptr [rdx+864+16]
		vmulps xmm2,xmm2,XMMWORD ptr [rdx+864+32]
		vfmadd231ps xmm2,xmm3,XMMWORD ptr [rdx+864+48]
		
		vpshufd xmm4,xmm7,0
		vpshufd xmm5,xmm7,85
		vpshufd xmm6,xmm7,170
		vpshufd xmm7,xmm7,255
		
		vmulps xmm4,xmm4,XMMWORD ptr [rdx+864+64]
		vfmadd231ps xmm4,xmm5,XMMWORD ptr [rdx+864+80]
		vmulps xmm6,xmm6,XMMWORD ptr [rdx+864+96]
		vfmadd231ps xmm6,xmm7,XMMWORD ptr [rdx+864+112]		
		
		vaddps xmm0,xmm0,xmm2
		vaddps xmm4,xmm4,xmm6
		vaddps xmm0,xmm0,xmm4
		vaddps xmm0,xmm0,XMMWORD ptr [rdx+864+128]
		vmovhlps xmm1,xmm1,xmm0
		vmaxps xmm0,xmm0,xmm1
		vpshuflw xmm1,xmm0,14
		vcomiss xmm1,xmm0
		jbe short finish_1a
		xor rax,rax
finish_1a:
		mov BYTE PTR[r8],al

	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,40

	vzeroupper

	ret

computeNetwork0_AVX512 endp


;processLine0_AVX512_ASM proc tempu:dword,width_:dword,dstp:dword,src3p:dword,src_pitch:dword,val_min_max:dword
; tempu = rcx
; width_ = edx
; dstp = r8
; src3p = r9

processLine0_AVX512_ASM proc public frame

src_pitch equ dword ptr[rbp+48]
val_min_max equ qword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	.endprolog

		mov rax,rcx
		mov ecx,edx
		movsxd rdx,src_pitch
		mov r10,val_min_max
		
		lea r11,[r9+rdx*4]
		
		vmovdqa64 zmm18,ZMMWORD ptr w_19
		vmovdqa64 zmm19,ZMMWORD ptr w_3
		vmovdqa64 zmm20,ZMMWORD ptr ub_1
		vmovdqa64 zmm21,ZMMWORD ptr uw_16		
		vmovdqa64 zmm22,ZMMWORD ptr[r10]
		vmovdqa64 zmm23,ZMMWORD ptr[r10+64]
		vpxord zmm16,zmm16,zmm16
		vpxord zmm17,zmm17,zmm17
		
xloop:
		vmovdqu64 zmm0,ZMMWORD ptr[r9+rdx*2]
		vmovdqu64 zmm1,ZMMWORD ptr[r11]
		vpunpckhbw zmm2,zmm0,zmm17
		vpunpckhbw zmm3,zmm1,zmm17
		vpunpcklbw zmm0,zmm0,zmm17
		vpunpcklbw zmm1,zmm1,zmm17
		vpaddw zmm0,zmm0,zmm1
		vpaddw zmm2,zmm2,zmm3
		vpmullw zmm0,zmm0,zmm18
		vpmullw zmm2,zmm2,zmm18

		vmovdqu64 zmm1,ZMMWORD ptr[r9]
		vmovdqu64 zmm3,ZMMWORD ptr[r11+rdx*2]
		vpunpckhbw zmm4,zmm1,zmm17
		vpunpckhbw zmm5,zmm3,zmm17
		vpunpcklbw zmm1,zmm1,zmm17
		vpunpcklbw zmm3,zmm3,zmm17
		vpaddw zmm1,zmm1,zmm3
		vpaddw zmm4,zmm4,zmm5
		vpmullw zmm1,zmm1,zmm19
		vpmullw zmm4,zmm4,zmm19

		vmovdqu64 zmm5,ZMMWORD ptr[rax]
		vpsubusw zmm0,zmm0,zmm1
		vpsubusw zmm2,zmm2,zmm4
		vpxord zmm5,zmm5,zmm20
		vpaddusw zmm0,zmm0,zmm21
		vpaddusw zmm2,zmm2,zmm21
		vpsadbw zmm5,zmm5,zmm17		
		vpsraw zmm0,zmm0,5
		vpsraw zmm2,zmm2,5

		vmovdqa64 zmm3,zmm5
		vpminuw zmm0,zmm0,zmm23
		vpsrldq zmm5,zmm5,8
		vpminuw zmm2,zmm2,zmm23
		vpaddusw zmm5,zmm5,zmm3
		vpmaxuw zmm0,zmm0,zmm22
		vpmaxuw zmm2,zmm2,zmm22

		vextracti32x8 ymm1,zmm5,1

		vextracti128 xmm3,ymm5,1
		vextracti128 xmm4,ymm1,1

		vpackuswb zmm0,zmm0,zmm2

		vpaddusw xmm5,xmm5,xmm3
		vpaddusw xmm1,xmm1,xmm4

		vmovdqu64 ZMMWORD ptr[r8],zmm0

		vpaddusw xmm16,xmm16,xmm5		
		vpaddusw xmm16,xmm16,xmm1
		
		add r9,64
		add r11,64
		add rax,64
		add r8,64
		sub ecx,64
		jnz xloop
					
		xor  rax,rax
		vmovd eax,xmm16		

	pop rbp

	vzeroupper

	ret

processLine0_AVX512_ASM endp


;processLine0_AVX512_ASM_16 proc tempu:dword,width_:dword,dstp:dword,src3p:dword,src_pitch:dword,val_min_max:dword
; tempu = rcx
; width_ = edx
; dstp = r8
; src3p = r9

processLine0_AVX512_ASM_16 proc public frame

src_pitch equ dword ptr[rbp+48]
val_min_max equ qword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	.endprolog
		
		mov rax,rcx
		mov ecx,edx
		movsxd rdx,src_pitch
		mov r10,val_min_max
			
		lea r11,[r9+rdx*4]
		
		vmovdqa64 zmm18,ZMMWORD ptr d_19
		vmovdqa64 zmm19,ZMMWORD ptr d_3
		vmovdqa32 ymm20,YMMWORD ptr ub_1
		vmovdqa64 zmm21,ZMMWORD ptr ud_16
		vmovdqa64 zmm22,ZMMWORD ptr[r10]
		vmovdqa64 zmm23,ZMMWORD ptr[r10+64]

		vpxord zmm16,zmm16,zmm16
		vpxord zmm17,zmm17,zmm17
		
xloop_16:
		vmovdqa64 zmm0,ZMMWORD ptr[r9+rdx*2]
		vmovdqa64 zmm1,ZMMWORD ptr[r11]
		vpunpckhwd zmm2,zmm0,zmm17
		vpunpckhwd zmm3,zmm1,zmm17
		vpunpcklwd zmm0,zmm0,zmm17
		vpunpcklwd zmm1,zmm1,zmm17
		vpaddd zmm0,zmm0,zmm1
		vpaddd zmm2,zmm2,zmm3
		vpmulld zmm0,zmm0,zmm18
		vpmulld zmm2,zmm2,zmm18

		vmovdqa64 zmm1,ZMMWORD ptr[r9]
		vmovdqa64 zmm3,ZMMWORD ptr[r11+rdx*2]
		vpunpckhwd zmm4,zmm1,zmm17
		vpunpckhwd zmm5,zmm3,zmm17
		vpunpcklwd zmm1,zmm1,zmm17
		vpunpcklwd zmm3,zmm3,zmm17
		vpaddd zmm1,zmm1,zmm3
		vpaddd zmm4,zmm4,zmm5
		vpmulld zmm1,zmm1,zmm19
		vpmulld zmm4,zmm4,zmm19

		vmovdqa ymm5,YMMWORD ptr[rax]
		vpsubd zmm0,zmm0,zmm1
		vpsubd zmm2,zmm2,zmm4
		vpxord ymm5,ymm5,ymm20
		vpaddd zmm0,zmm0,zmm21
		vpaddd zmm2,zmm2,zmm21
		vpsadbw ymm5,ymm5,ymm17
		vpsrad zmm0,zmm0,5		
		vpsrad zmm2,zmm2,5
		vmovdqa ymm3,ymm5

		vpackusdw zmm0,zmm0,zmm2
		vpsrldq ymm5,ymm5,8
		vpminuw zmm0,zmm0,zmm23
		vpaddusw ymm5,ymm5,ymm3
		vpmaxuw zmm0,zmm0,zmm22

		vextracti128 xmm3,ymm5,1
		vmovdqa64 ZMMWORD ptr [r8],zmm0

		vpaddusw xmm5,xmm5,xmm3
		vpaddusw xmm16,xmm16,xmm5
		
		add r9,64
		add r11,64
		add rax,64
		add r8,64
		sub ecx,32
		jnz xloop_16
					
		xor  rax,rax
		vmovd eax,xmm16		
		
	pop rbp
	
	vzeroupper
	
	ret
		
processLine0_AVX512_ASM_16 endp


;processLine0_AVX512_ASM_32 proc tempu:dword,width_:dword,dstp:dword,src3p:dword,src_pitch:dword
; tempu = rcx
; width_ = edx
; dstp = r8
; src3p = r9

processLine0_AVX512_ASM_32 proc public frame

src_pitch equ dword ptr[rbp+48]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	.endprolog

		mov rax,rcx
		mov ecx,edx
		movsxd rdx,src_pitch
		
		lea r10,[r9+rdx*4]

		vpxord zmm5,zmm5,zmm5
		vpxord zmm16,zmm16,zmm16
		vmovaps zmm17,ZMMWORD ptr f_19
		vmovaps zmm18,ZMMWORD ptr f_3
		vmovdqa32 xmm19,XMMWORD ptr uw_1
		
xloop_32:
		vmovdqa xmm4,XMMWORD ptr [rax]
		vmovaps zmm2,ZMMWORD ptr[r9]		
		vmovaps zmm0,ZMMWORD ptr[r9+rdx*2]
		vpunpckhbw xmm20,xmm4,xmm16
		vpunpcklbw xmm4,xmm4,xmm16

		vmovaps zmm1,ZMMWORD ptr[r10]
		vmovaps zmm3,ZMMWORD ptr[r10+rdx*2]		
		vaddps zmm0,zmm0,zmm1
		vpxord xmm4,xmm4,xmm19
		vpxord xmm20,xmm20,xmm19
		vaddps zmm2,zmm2,zmm3		
		vpsadbw xmm4,xmm4,xmm16
		vpsadbw xmm20,xmm20,xmm16

		vmulps zmm0,zmm0,zmm17
		vmovdqa xmm3,xmm4
		vmovdqa32 xmm1,xmm20
		vmulps zmm2,zmm2,zmm18
		vpsrldq xmm4,xmm4,8
		vpsrldq xmm20,xmm20,8
		vpaddusw xmm4,xmm4,xmm3
		vpaddusw xmm20,xmm20,xmm1
		vsubps zmm0,zmm0,zmm2

		vpaddusw xmm4,xmm4,xmm20
		vmovaps ZMMWORD ptr[r8],zmm0
		vpaddusw xmm5,xmm5,xmm4

		add r9,64
		add r10,64
		add rax,16
		add r8,64
		sub ecx,16
		jnz xloop_32
				
		xor  rax,rax
		vmovd eax,xmm5

	pop rbp

	vzeroupper		

	ret

processLine0_AVX512_ASM_32 endp


; From FMA3
;weightedAvgElliottMul5_m16_AVX512 proc ptr_w:dword,n:dword,mstd:dword
; ptr_w = rcx
; n = edx
; mstd = r8

weightedAvgElliottMul5_m16_AVX512 proc public frame

	sub rsp,24
	.allocstack 24
	vmovdqa XMMWORD ptr[rsp],xmm8
	.savexmm128 xmm8,0
	.endprolog

		mov edx,edx  ; RAZ high

		vmovdqa64 zmm16,ZMMWORD ptr sign_bits_f_32
		vmovaps zmm17,ZMMWORD ptr ones_f_32
	
		lea r10,[rcx+rdx*4]
		xor r9,r9
		vxorps zmm0,zmm0,zmm0
		vxorps zmm1,zmm1,zmm1
		
nloop_52:
		vmovaps zmm2,ZMMWORD ptr [rcx+r9*4]
		vmovaps zmm4,ZMMWORD ptr [r10+r9*4]
		vaddps zmm0,zmm0,zmm2
		vandps zmm5,zmm4,zmm16
		vaddps zmm5,zmm5,zmm17

		vextractf32x8 ymm8,zmm5,1
		vrcpps ymm5,ymm5
		vrcpps ymm8,ymm8
		vinsertf32x8 zmm5,zmm5,ymm8,1

		vmulps zmm4,zmm4,zmm5
		vfmadd231ps zmm1,zmm2,zmm4
		
		add r9,16
		sub edx,16
		jnz short nloop_52

		vextractf32x8 ymm2,zmm0,1
		vextractf32x8 ymm3,zmm1,1
		
		vaddps ymm0,ymm0,ymm2
		vaddps ymm1,ymm1,ymm3

		vextractf128 xmm2,ymm0,1
		vextractf128 xmm3,ymm1,1
		vaddps xmm0,xmm0,xmm2
		vaddps xmm1,xmm1,xmm3
		
		vmovhlps xmm2,xmm2,xmm0
		vmovhlps xmm3,xmm3,xmm1
		vaddps xmm0,xmm0,xmm2
		vaddps xmm1,xmm1,xmm3
		vpshuflw xmm2,xmm0,14
		vpshuflw xmm3,xmm1,14
		vaddss xmm0,xmm0,xmm2
		vaddss xmm1,xmm1,xmm3
		vcomiss xmm0,dword ptr min_weight_sum
		jbe short nodiv2
		vmulss xmm1,xmm1,dword ptr five_f
		vrcpss xmm0,xmm0,xmm0
		vmulss xmm1,xmm1,xmm0
		jmp short finish_52
nodiv2:
		vxorps xmm1,xmm1,xmm1
finish_52:
		vmulss xmm1,xmm1,dword ptr[r8+4]
		vaddss xmm1,xmm1,dword ptr[r8]
		vaddss xmm1,xmm1,dword ptr[r8+12]
		vmovss dword ptr[r8+12],xmm1

	vmovdqa xmm8,XMMWORD ptr[rsp]	
	add rsp,24

	vzeroupper

	ret

weightedAvgElliottMul5_m16_AVX512 endp


; From FMA3
;dotProd_m32_m16_AVX512 proc data_:dword,weights:dword,vals:dword,n:dword,len:dword,istd:dword
; data_ = rcx
; weights = rdx
; vals = r8
; n = r9d

dotProd_m32_m16_AVX512 proc public frame

len equ dword ptr[rbp+48]
istd equ qword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push r12
	.pushreg r12
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	.endprolog
		
		mov r11,rdx
		mov rax,r8
		mov ebx,r9d
		mov r12d,len
		mov r10,rcx
		
		mov r13,128
		mov r14,512
nloop_2:
		mov rcx,r10
		vxorps zmm0,zmm0,zmm0
		vxorps zmm1,zmm1,zmm1
		vxorps zmm2,zmm2,zmm2
		vxorps zmm3,zmm3,zmm3
		mov edx,r12d
lloop_2:
		vmovaps zmm17,ZMMWORD ptr[rcx]
		vfmadd231ps zmm0,zmm17,ZMMWORD ptr[r11]
		vfmadd231ps zmm1,zmm17,ZMMWORD ptr[r11+64]
		vfmadd231ps zmm2,zmm17,ZMMWORD ptr[r11+128]
		vfmadd231ps zmm3,zmm17,ZMMWORD ptr[r11+192]
		
		vmovaps zmm17,ZMMWORD ptr[rcx+64]
		vfmadd231ps zmm0,zmm17,ZMMWORD ptr[r11+256]
		vfmadd231ps zmm1,zmm17,ZMMWORD ptr[r11+320]
		vfmadd231ps zmm2,zmm17,ZMMWORD ptr[r11+384]
		vfmadd231ps zmm3,zmm17,ZMMWORD ptr[r11+448]
				
		add rcx,r13
		add r11,r14
		sub edx,32
		jnz short lloop_2
		
		vextractf32x8 ymm4,zmm0,1
		vextractf32x8 ymm5,zmm1,1
		vextractf32x8 ymm16,zmm2,1
		vextractf32x8 ymm17,zmm3,1
		
		vaddps ymm0,ymm0,ymm4
		vaddps ymm1,ymm1,ymm5
		vaddps ymm2,ymm2,ymm16
		vaddps ymm3,ymm3,ymm17
				
		vextractf128 xmm4,ymm0,1
		vextractf128 xmm5,ymm1,1
		vextractf32x4 xmm16,ymm2,1
		vextractf32x4 xmm17,ymm3,1
		vaddps xmm0,xmm0,xmm4
		vaddps xmm1,xmm1,xmm5
		vaddps xmm2,xmm2,xmm16
		vaddps xmm3,xmm3,xmm17								
				
		haddps xmm0,xmm1
		haddps xmm2,xmm3
		haddps xmm0,xmm2		
		
		vmovaps XMMWORD ptr[rax],xmm0
		add rax,16
		sub ebx,4
		jnz nloop_2
		
		mov rcx,istd
		mov rax,r8
		vmovss xmm17,dword ptr[rcx]
		mov edx,r9d
		vshufps xmm17,xmm17,xmm17,0
		xor rcx,rcx
		vinsertf32x4 ymm17,ymm17,xmm17,1
aloop_2:
		vmulps ymm0,ymm17,YMMWORD ptr[rax+rcx*4]
		vmulps ymm2,ymm17,YMMWORD ptr[rax+rcx*4+32]
		vaddps ymm0,ymm0,YMMWORD ptr[r11+rcx*4]
		vaddps ymm2,ymm2,YMMWORD ptr[r11+rcx*4+32]
		vmovaps YMMWORD ptr[rax+rcx*4],ymm0
		vmovaps YMMWORD ptr[rax+rcx*4+32],ymm2
		add rcx,16
		sub edx,16
		jnz short aloop_2

	vzeroupper
	
	pop r14
	pop r13
	pop r12
	pop rbx
	pop rbp

	ret

dotProd_m32_m16_AVX512 endp


; From FMA3
;dotProd_m48_m16_AVX512 proc data_:dword,weights:dword,vals:dword,n:dword,len:dword,istd:dword
; data_ = rcx
; weights = rdx
; vals = r8
; n = r9d

dotProd_m48_m16_AVX512 proc public frame

len equ dword ptr[rbp+48]
istd equ qword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push r12
	.pushreg r12
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	.endprolog
		
		mov r11,rdx
		mov rax,r8
		mov ebx,r9d
		mov r12d,len
		mov r10,rcx
		
		mov r13,192
		mov r14,768
nloop2_2:
		mov rcx,r10
		vxorps ymm0,ymm0,ymm0
		vxorps ymm1,ymm1,ymm1
		vxorps ymm2,ymm2,ymm2
		vxorps ymm3,ymm3,ymm3
		mov edx,r12d
lloop2_2:
		vmovaps zmm17,ZMMWORD ptr[rcx]
		vfmadd231ps zmm0,zmm17,ZMMWORD ptr[r11]
		vfmadd231ps zmm1,zmm17,ZMMWORD ptr[r11+64]
		vfmadd231ps zmm2,zmm17,ZMMWORD ptr[r11+128]
		vfmadd231ps zmm3,zmm17,ZMMWORD ptr[r11+192]
				
		vmovaps zmm17,ZMMWORD ptr[rcx+64]
		vfmadd231ps zmm0,zmm17,ZMMWORD ptr[r11+256]
		vfmadd231ps zmm1,zmm17,ZMMWORD ptr[r11+320]
		vfmadd231ps zmm2,zmm17,ZMMWORD ptr[r11+384]
		vfmadd231ps zmm3,zmm17,ZMMWORD ptr[r11+448]
		
		vmovaps zmm17,ZMMWORD ptr[rcx+128]
		vfmadd231ps zmm0,zmm17,ZMMWORD ptr[r11+512]
		vfmadd231ps zmm1,zmm17,ZMMWORD ptr[r11+576]
		vfmadd231ps zmm2,zmm17,ZMMWORD ptr[r11+640]
		vfmadd231ps zmm3,zmm17,ZMMWORD ptr[r11+704]
						
		add rcx,r13
		add r11,r14
		sub edx,48
		jnz short lloop2_2

		vextractf32x8 ymm4,zmm0,1
		vextractf32x8 ymm5,zmm1,1
		vextractf32x8 ymm16,zmm2,1
		vextractf32x8 ymm17,zmm3,1
		
		vaddps ymm0,ymm0,ymm4
		vaddps ymm1,ymm1,ymm5
		vaddps ymm2,ymm2,ymm16
		vaddps ymm3,ymm3,ymm17

		vextractf128 xmm4,ymm0,1
		vextractf128 xmm5,ymm1,1
		vextractf32x4 xmm16,ymm2,1
		vextractf32x4 xmm17,ymm3,1
		vaddps xmm0,xmm0,xmm4
		vaddps xmm1,xmm1,xmm5
		vaddps xmm2,xmm2,xmm16
		vaddps xmm3,xmm3,xmm17								
				
		haddps xmm0,xmm1
		haddps xmm2,xmm3
		haddps xmm0,xmm2		
		
		vmovaps XMMWORD ptr[rax],xmm0
		add rax,16
		sub ebx,4
		jnz nloop2_2
		
		mov rcx,istd
		mov rax,r8
		vmovss xmm17,dword ptr[rcx]
		mov edx,r9d
		vshufps xmm17,xmm17,xmm17,0
		xor rcx,rcx
		vinsertf32x4 ymm17,ymm17,xmm17,1
aloop2_2:
		vmulps ymm0,ymm17,YMMWORD ptr[rax+rcx*4]
		vmulps ymm2,ymm17,YMMWORD ptr[rax+rcx*4+32]
		vaddps ymm0,ymm0,YMMWORD ptr[r11+rcx*4]
		vaddps ymm2,ymm2,YMMWORD ptr[r11+rcx*4+32]
		vmovaps YMMWORD ptr[rax+rcx*4],ymm0
		vmovaps YMMWORD ptr[rax+rcx*4+32],ymm2
		add rcx,16
		sub edx,16
		jnz short aloop2_2

	vzeroupper
	
	pop r14
	pop r13
	pop r12
	pop rbx
	pop rbp

	ret

dotProd_m48_m16_AVX512 endp


;dotProd_m32_m16_i16_AVX512 proc dataf:dword,weightsf:dword,vals:dword,n:dword,len:dword,istd:dword
; dataf = rcx
; weightsf = rdx
; vals = r8
; n = r9d

dotProd_m32_m16_i16_AVX512 proc public frame

len equ dword ptr[rbp+48]
istd equ qword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push r12
	.pushreg r12
	push r14
	.pushreg r14
	.endprolog
		
		mov r11,rdx
		mov rax,r8
		mov ebx,r9d
		mov r12d,len
		mov r10,rcx
		
		mov r14,256
nloop_3:
		mov rcx,r10
		vpxord zmm0,zmm0,zmm0
		vpxord zmm1,zmm1,zmm1
		vpxord zmm2,zmm2,zmm2
		vpxord zmm3,zmm3,zmm3
		mov edx,r12d
lloop_3:
		vmovdqa64 zmm17,ZMMWORD ptr [rcx]
		vpmaddwd zmm4,zmm17,ZMMWORD ptr [r11]
		vpmaddwd zmm5,zmm17,ZMMWORD ptr [r11+64]
		vpmaddwd zmm16,zmm17,ZMMWORD ptr [r11+128]
		vpmaddwd zmm17,zmm17,ZMMWORD ptr [r11+192]
		vpaddd zmm0,zmm0,zmm4
		vpaddd zmm1,zmm1,zmm5
		vpaddd zmm2,zmm2,zmm16
		vpaddd zmm3,zmm3,zmm17

		add rcx,64
		add r11,r14
		sub edx,32
		jnz short lloop_3

		vextracti32x8 ymm4,zmm0,1
		vextracti32x8 ymm5,zmm1,1
		vextracti32x8 ymm16,zmm2,1
		vextracti32x8 ymm17,zmm3,1

		vpaddd ymm0,ymm0,ymm4
		vpaddd ymm1,ymm1,ymm5
		vpaddd ymm2,ymm2,ymm16
		vpaddd ymm3,ymm3,ymm17

		vextracti128 xmm4,ymm0,1
		vextracti128 xmm5,ymm1,1
		vextracti32x4 xmm16,ymm2,1
		vextracti32x4 xmm17,ymm3,1
		vpaddd xmm0,xmm0,xmm4
		vpaddd xmm1,xmm1,xmm5
		vpaddd xmm2,xmm2,xmm16
		vpaddd xmm3,xmm3,xmm17		
		vpunpckhqdq xmm4,xmm0,xmm1
		vpunpckhqdq xmm5,xmm2,xmm3
		vpunpcklqdq xmm0,xmm0,xmm1
		vpunpcklqdq xmm2,xmm2,xmm3
		vpaddd xmm0,xmm0,xmm4
		vpaddd xmm2,xmm2,xmm5
		vshufps xmm16,xmm0,xmm2,221
		vshufps xmm0,xmm0,xmm2,136
		vpaddd xmm16,xmm16,xmm0
		vmovdqa32 XMMWORD ptr [rax],xmm16
		add rax,16
		sub ebx,4
		jnz nloop_3

		mov rcx,istd
		mov rax,r8
		vmovss xmm17,dword ptr[rcx]
		mov edx,r9d
		vpshufd xmm17,xmm17,0
		xor rcx,rcx
aloop_3:
		vmovdqa ymm0,YMMWORD ptr[rax+rcx*4]
		vmovdqa ymm2,YMMWORD ptr[rax+rcx*4+32]
		vcvtdq2ps ymm0,ymm0
		vcvtdq2ps ymm2,ymm2
		vextractf128 xmm1,ymm0,1
		vextractf128 xmm3,ymm2,1
		vmulps xmm0,xmm0,XMMWORD ptr[r11+rcx*8]
		vmulps xmm1,xmm1,XMMWORD ptr[r11+rcx*8+32]
		vmulps xmm2,xmm2,XMMWORD ptr[r11+rcx*8+64]
		vmulps xmm3,xmm3,XMMWORD ptr[r11+rcx*8+96]
		vmulps xmm0,xmm0,xmm17
		vmulps xmm1,xmm1,xmm17
		vmulps xmm2,xmm2,xmm17
		vmulps xmm3,xmm3,xmm17
		vaddps xmm0,xmm0,XMMWORD ptr[r11+rcx*8+16]
		vaddps xmm1,xmm1,XMMWORD ptr[r11+rcx*8+48]
		vaddps xmm2,xmm2,XMMWORD ptr[r11+rcx*8+80]
		vaddps xmm3,xmm3,XMMWORD ptr[r11+rcx*8+112]
		vmovaps XMMWORD ptr[rax+rcx*4],xmm0
		vmovaps XMMWORD ptr[rax+rcx*4+16],xmm1
		vmovaps XMMWORD ptr[rax+rcx*4+32],xmm2
		vmovaps XMMWORD ptr[rax+rcx*4+48],xmm3
		add rcx,16
		sub edx,16
		jnz aloop_3

	vzeroupper
		
	pop r14
	pop r12
	pop rbx
	pop rbp		

	ret

dotProd_m32_m16_i16_AVX512 endp


;dotProd_m64_m16_i16_AVX512 proc dataf:dword,weightsf:dword,vals:dword,n:dword,len:dword,istd:dword
; dataf = rcx
; weightsf = rdx
; vals = r8
; n = r9d

dotProd_m64_m16_i16_AVX512 proc public frame

len equ dword ptr[rbp+48]
istd equ qword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push r12
	.pushreg r12
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	.endprolog
		
		mov r11,rdx
		mov rax,r8
		mov ebx,r9d
		mov r12d,len
		mov r10,rcx
		
		mov r13,128
		mov r14,512
nloop_4:
		mov rcx,r10
		vpxord zmm0,zmm0,zmm0
		vpxord zmm1,zmm1,zmm1
		vpxord zmm2,zmm2,zmm2
		vpxord zmm3,zmm3,zmm3
		mov edx,r12d
lloop_4:
		vmovdqa64 zmm17,ZMMWORD ptr [rcx]
		vpmaddwd zmm4,zmm17,ZMMWORD ptr [r11]
		vpmaddwd zmm5,zmm17,ZMMWORD ptr [r11+64]
		vpmaddwd zmm16,zmm17,ZMMWORD ptr [r11+128]
		vpmaddwd zmm17,zmm17,ZMMWORD ptr [r11+192]
		vpaddd zmm0,zmm0,zmm4
		vpaddd zmm1,zmm1,zmm5
		vpaddd zmm2,zmm2,zmm16
		vpaddd zmm3,zmm3,zmm17

		vmovdqa64 zmm17,ZMMWORD ptr [rcx+64]
		vpmaddwd zmm4,zmm17,ZMMWORD ptr [r11+256]
		vpmaddwd zmm5,zmm17,ZMMWORD ptr [r11+320]
		vpmaddwd zmm16,zmm17,ZMMWORD ptr [r11+384]
		vpmaddwd zmm17,zmm17,ZMMWORD ptr [r11+448]
		vpaddd zmm0,zmm0,zmm4
		vpaddd zmm1,zmm1,zmm5
		vpaddd zmm2,zmm2,zmm16
		vpaddd zmm3,zmm3,zmm17

		add rcx,r13
		add r11,r14
		sub edx,64
		jnz short lloop_4

		vextracti32x8 ymm4,zmm0,1
		vextracti32x8 ymm5,zmm1,1
		vextracti32x8 ymm16,zmm2,1
		vextracti32x8 ymm17,zmm3,1

		vpaddd ymm0,ymm0,ymm4
		vpaddd ymm1,ymm1,ymm5
		vpaddd ymm2,ymm2,ymm16
		vpaddd ymm3,ymm3,ymm17

		vextracti128 xmm4,ymm0,1
		vextracti128 xmm5,ymm1,1
		vextracti32x4 xmm16,ymm2,1
		vextracti32x4 xmm17,ymm3,1
		vpaddd xmm0,xmm0,xmm4
		vpaddd xmm1,xmm1,xmm5
		vpaddd xmm2,xmm2,xmm16
		vpaddd xmm3,xmm3,xmm17		
		vpunpckhqdq xmm4,xmm0,xmm1
		vpunpckhqdq xmm5,xmm2,xmm3
		vpunpcklqdq xmm0,xmm0,xmm1
		vpunpcklqdq xmm2,xmm2,xmm3
		vpaddd xmm0,xmm0,xmm4
		vpaddd xmm2,xmm2,xmm5
		vshufps xmm16,xmm0,xmm2,221
		vshufps xmm0,xmm0,xmm2,136
		vpaddd xmm16,xmm16,xmm0
		vmovdqa32 XMMWORD ptr[rax],xmm16
		
		add rax,16
		sub ebx,4
		jnz nloop_4

		mov rcx,istd
		mov rax,r8
		vmovss xmm17,dword ptr[rcx]
		mov edx,r9d
		vpshufd xmm17,xmm17,0
		xor rcx,rcx
aloop_4:
		vmovdqa ymm0,YMMWORD ptr[rax+rcx*4]
		vmovdqa ymm2,YMMWORD ptr[rax+rcx*4+32]
		vcvtdq2ps ymm0,ymm0
		vcvtdq2ps ymm2,ymm2
		vextractf128 xmm1,ymm0,1
		vextractf128 xmm3,ymm2,1
		vmulps xmm0,xmm0,XMMWORD ptr[r11+rcx*8]
		vmulps xmm1,xmm1,XMMWORD ptr[r11+rcx*8+32]
		vmulps xmm2,xmm2,XMMWORD ptr[r11+rcx*8+64]
		vmulps xmm3,xmm3,XMMWORD ptr[r11+rcx*8+96]
		vmulps xmm0,xmm0,xmm17
		vmulps xmm1,xmm1,xmm17
		vmulps xmm2,xmm2,xmm17
		vmulps xmm3,xmm3,xmm17
		vaddps xmm0,xmm0,XMMWORD ptr[r11+rcx*8+16]
		vaddps xmm1,xmm1,XMMWORD ptr[r11+rcx*8+48]
		vaddps xmm2,xmm2,XMMWORD ptr[r11+rcx*8+80]
		vaddps xmm3,xmm3,XMMWORD ptr[r11+rcx*8+112]
		vmovaps XMMWORD ptr[rax+rcx*4],xmm0
		vmovaps XMMWORD ptr[rax+rcx*4+16],xmm1
		vmovaps XMMWORD ptr[rax+rcx*4+32],xmm2
		vmovaps XMMWORD ptr[rax+rcx*4+48],xmm3
		add rcx,16
		sub edx,16
		jnz aloop_4

	vzeroupper
		
	pop r14
	pop r13
	pop r12
	pop rbx
	pop rbp		

	ret

dotProd_m64_m16_i16_AVX512 endp


; From FMA3
;e0_m16_AVX512 proc ptr_s:dword,n:dword
; ptr_s = rcx
; n = edx

e0_m16_AVX512 proc public frame

	.endprolog

		vmovaps zmm2,ZMMWORD ptr exp_hi
		vmovaps zmm3,ZMMWORD ptr exp_lo
		vmovaps zmm4,ZMMWORD ptr e0_mult
		vmovaps zmm5,ZMMWORD ptr e0_bias

		mov r10,128

eloop16_2:
		vmovaps zmm0,ZMMWORD ptr[rcx]
		vmovaps zmm1,ZMMWORD ptr[rcx+64]
		vminps zmm0,zmm0,zmm2
		vminps zmm1,zmm1,zmm2
		vmaxps zmm0,zmm0,zmm3
		vmaxps zmm1,zmm1,zmm3
		
		vfmadd213ps zmm0,zmm4,zmm5
		vfmadd213ps zmm1,zmm4,zmm5
				
		vcvtps2dq zmm0,zmm0
		vcvtps2dq zmm1,zmm1
		vmovaps ZMMWORD ptr[rcx],zmm0
		vmovaps ZMMWORD ptr[rcx+64],zmm1

		add rcx,r10
		sub edx,32
		
		jnz short eloop16_2		
		
	vzeroupper
		
	ret

e0_m16_AVX512 endp


;e1_m16_AVX512 proc ptr_s:dword,n:dword
; ptr_s = rcx
; n = edx

e1_m16_AVX512 proc public frame

	.endprolog
	
		vmovaps zmm3,ZMMWORD ptr exp_hi
		vmovaps zmm4,ZMMWORD ptr exp_lo
		vmovaps zmm5,ZMMWORD ptr e1_scale
		vmovaps zmm16,ZMMWORD ptr e1_bias
		vmovaps zmm17,ZMMWORD ptr e1_c1
		vmovaps zmm18,ZMMWORD ptr e1_c2
		vmovaps zmm19,ZMMWORD ptr e1_c0

eloop8:
		vmovaps zmm0,ZMMWORD ptr[rcx]
		vminps zmm0,zmm0,zmm3
		vmaxps zmm0,zmm0,zmm4
		vmulps zmm0,zmm0,zmm5
		vmovaps zmm1,zmm0
		vaddps zmm0,zmm0,zmm16
		vpslld zmm2,zmm0,23
		vsubps zmm0,zmm0,zmm16
		vsubps zmm1,zmm1,zmm0
		vmulps zmm0,zmm1,zmm17
		vmulps zmm1,zmm1,zmm1
		vmulps zmm1,zmm1,zmm18
		vaddps zmm0,zmm0,zmm19
		vaddps zmm0,zmm0,zmm1
		vpaddd zmm0,zmm0,zmm2
		vmovaps ZMMWORD ptr[rcx],zmm0
		add rcx,64
		sub edx,16
		jnz short eloop8

	vzeroupper

	ret

e1_m16_AVX512 endp


;e2_m16_AVX512 proc ptr_s:dword,n:dword
; ptr_s = rcx
; n = edx

e2_m16_AVX512 proc public frame

	sub rsp,24
	.allocstack 24
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	.endprolog
	
		vmovaps zmm17,ZMMWORD ptr exp_hi
		vmovaps zmm18,ZMMWORD ptr exp_lo
		vmovaps zmm19,ZMMWORD ptr exp_rln2
		vmovaps zmm20,ZMMWORD ptr am_0p5
		vmovaps zmm21,ZMMWORD ptr epi32_1
		vmovaps zmm22,ZMMWORD ptr exp_c2
		vmovaps zmm23,ZMMWORD ptr exp_c1
		vmovaps zmm24,ZMMWORD ptr exp_q0
		vmovaps zmm25,ZMMWORD ptr am_1

eloop4:
		vmovaps zmm0,ZMMWORD ptr[rcx]		
		vminps zmm0,zmm0,zmm17
		vmaxps zmm0,zmm0,zmm18
		vmulps zmm1,zmm0,zmm19
		vxorps zmm2,zmm2,zmm2
		vaddps zmm1,zmm1,zmm20

		vextractf32x8 ymm3,zmm1,1		
		vextractf32x8 ymm4,zmm2,1		
		vcmpnltps ymm2,ymm2,ymm1
		vcmpnltps ymm4,ymm4,ymm3
		vinsertf32x8 zmm2,zmm2,ymm4,1

		vpandd zmm2,zmm2,zmm21
		vcvttps2dq zmm1,zmm1
		vpsubd zmm1,zmm1,zmm2
		vmovaps zmm5,zmm23
		vcvtdq2ps zmm3,zmm1
		vmulps zmm4,zmm3,zmm22
		vmulps zmm5,zmm5,zmm3
		vsubps zmm0,zmm0,zmm4
		vsubps zmm0,zmm0,zmm5
		vpaddd zmm1,zmm1,ZMMWORD ptr epi32_0x7f
		vmovaps zmm2,zmm0
		vmulps zmm0,zmm0,zmm0
		vmulps zmm6,zmm0,zmm24
		vmulps zmm4,zmm0,ZMMWORD ptr exp_p0
		vaddps zmm6,zmm6,ZMMWORD ptr exp_q1
		vaddps zmm4,zmm4,ZMMWORD ptr exp_p1
		vmulps zmm6,zmm6,zmm0
		vmulps zmm4,zmm4,zmm0
		vaddps zmm6,zmm6,ZMMWORD ptr exp_q2
		vmulps zmm4,zmm4,zmm2
		vmulps zmm6,zmm6,zmm0
		vaddps zmm2,zmm2,zmm4
		vaddps zmm6,zmm6,ZMMWORD ptr exp_q3
		vpslld zmm1,zmm1,23
		vsubps zmm6,zmm6,zmm2

		vextractf32x8 ymm5,zmm6,1		
		vrcpps ymm6,ymm6
		vrcpps ymm5,ymm5
		vinsertf32x8 zmm6,zmm6,ymm5,1

		vmulps zmm2,zmm2,zmm6
		vaddps zmm2,zmm2,zmm2
		vaddps zmm0,zmm2,zmm25
		vmulps zmm0,zmm0,zmm1		
		vmovaps ZMMWORD ptr[rcx],zmm0

		add rcx,64
		sub edx,16
		jnz eloop4

	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,24

	vzeroupper

	ret

e2_m16_AVX512 endp

end
