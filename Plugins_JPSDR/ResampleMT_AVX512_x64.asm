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

.code


;Resize_V_AVX512_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,width32:dword,src_pitch:dword,
;	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; width32 = r9d

Resize_V_AVX512_Planar_8bits_ASM proc public frame

src_pitch equ qword ptr[rbp+48]
kernel_size_2 equ dword ptr[rbp+56]
valmin equ qword ptr[rbp+64]
valmax equ qword ptr[rbp+72]
rounder equ qword ptr[rbp+80]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push rdi
	.pushreg rdi
	push rsi
	.pushreg rsi
	push r12
	.pushreg r12
	sub rsp,32
	.allocstack 32
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog		

	mov rsi,valmin
	vbroadcastss ymm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss ymm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss zmm7,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,32

Resize_V_AVX512_Planar_8bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa64 zmm0,zmm7
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_8bits_loop_2:
	vpmovzxbw zmm2,YMMWORD ptr [rsi]
	vpmovzxbw zmm4,YMMWORD ptr [rsi+rbx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[r8+rax]

	vpmaddwd zmm3,zmm3,zmm4
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
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
	vmovdqa XMMWORD ptr [rdi],xmm0
	vmovdqa XMMWORD ptr [rdi+16],xmm3
		
	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX512_Planar_8bits_loop_1

	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,32	

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX512_Planar_8bits_ASM endp


;Resize_V_AVX512_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,width32:dword,src_pitch:dword,
;	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; width32 = r9d

Resize_V_AVX512_Planar_10to14bits_ASM proc public frame

src_pitch equ qword ptr[rbp+48]
kernel_size_2 equ dword ptr[rbp+56]
valmin equ qword ptr[rbp+64]
valmax equ qword ptr[rbp+72]
rounder equ qword ptr[rbp+80]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push rdi
	.pushreg rdi
	push rsi
	.pushreg rsi
	push r12
	.pushreg r12
	sub rsp,32
	.allocstack 32
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog		

	mov rsi,valmin
	vbroadcastss zmm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss zmm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss zmm7,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,64

Resize_V_AVX512_Planar_10to14bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa64 zmm0,zmm7
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_10to14bits_loop_2:
	vmovdqa64 zmm2,ZMMWORD ptr [rsi]
	vmovdqa64 zmm4,ZMMWORD ptr [rsi+rbx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[r8+rax]

	vpmaddwd zmm3,zmm3,zmm4
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX512_Planar_10to14bits_loop_2

	vpsrad zmm0,zmm0,13 ;FPScale16bits = 13
	vpsrad zmm1,zmm1,13

	vpackusdw zmm0,zmm0,zmm1

	vpmaxuw zmm0,zmm0,zmm5
	vpminuw zmm0,zmm0,zmm6

	vmovdqa64 ZMMWORD ptr [rdi],zmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX512_Planar_10to14bits_loop_1

	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,32	

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX512_Planar_10to14bits_ASM endp


;Resize_V_AVX512_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,width32:dword,
;	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword
;	shifttosigned:dword,shiftfromsigned:dword

; src = rcx
; dst = rdx
; coeff = r8
; width32 = r9d

Resize_V_AVX512_Planar_16bits_ASM proc public frame

src_pitch equ qword ptr[rbp+48]
kernel_size_2 equ dword ptr[rbp+56]
valmin equ qword ptr[rbp+64]
valmax equ qword ptr[rbp+72]
rounder equ qword ptr[rbp+80]
shifttosigned equ qword ptr[rbp+88]
shiftfromsigned equ qword ptr[rbp+96]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push rdi
	.pushreg rdi
	push rsi
	.pushreg rsi
	push r12
	.pushreg r12
	sub rsp,64
	.allocstack 64
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	vmovdqa XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	.endprolog		

	mov rsi,valmin
	vbroadcastss zmm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss zmm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss zmm7,dword ptr[rsi]
	mov rsi,shifttosigned
	vbroadcastss zmm8,dword ptr[rsi]
	mov rsi,shiftfromsigned
	vbroadcastss zmm9,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,64

Resize_V_AVX512_Planar_16bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa64 zmm0,zmm7
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_16bits_loop_2:
	vmovdqa64 zmm2,ZMMWORD ptr [rsi]
	vmovdqa64 zmm4,ZMMWORD ptr [rsi+rbx]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vbroadcastss zmm4,dword ptr[r8+rax]

	vpaddw zmm3,zmm3,zmm8
	vpaddw zmm2,zmm2,zmm8

	vpmaddwd zmm3,zmm3,zmm4
	vpmaddwd zmm2,zmm2,zmm4

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX512_Planar_16bits_loop_2

	vpaddw zmm0,zmm0,zmm9
	vpaddw zmm1,zmm1,zmm9

	vpsrad zmm0,zmm0,13 ;FPScale16bits = 13
	vpsrad zmm1,zmm1,13

	vpackusdw zmm0,zmm0,zmm1

	vpmaxuw zmm0,zmm0,zmm5
	vpminuw zmm0,zmm0,zmm6

	vmovdqa64 ZMMWORD ptr [rdi],zmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX512_Planar_16bits_loop_1

	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]	
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,64

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX512_Planar_16bits_ASM endp


;Resize_V_AVX512_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,src_pitch:dword,
;	kernel_size_2:dword

; src = rcx
; dst = rdx
; coeff = r8
; width16 = r9d

Resize_V_AVX512_Planar_32bits_ASM proc public frame

src_pitch equ qword ptr[rbp+48]
kernel_size_2 equ dword ptr[rbp+56]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rbx
	.pushreg rbx
	push rdi
	.pushreg rdi
	push rsi
	.pushreg rsi
	push r12
	.pushreg r12
	.endprolog		

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,8 ; 2 coeff of 32 bits
	mov rbx,src_pitch
	mov r12,64

Resize_V_AVX512_Planar_32bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1

Resize_V_AVX512_Planar_32bits_loop_2:
	vmovaps zmm2,ZMMWORD ptr [rsi]
	vmovaps zmm3,ZMMWORD ptr [rsi+rbx]

	vbroadcastss zmm4,dword ptr[r8+rax]
	vbroadcastss zmm5,dword ptr[r8+rax+4]

	vfmadd231ps zmm0,zmm2,zmm4
	vfmadd231ps zmm1,zmm3,zmm5

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX512_Planar_32bits_loop_2

	vaddps zmm0,zmm0,zmm1

	vmovaps ZMMWORD ptr [rdi],zmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz short Resize_V_AVX512_Planar_32bits_loop_1

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX512_Planar_32bits_ASM endp

end
