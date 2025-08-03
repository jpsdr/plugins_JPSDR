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

.code


;Resize_V_AVX_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,src_pitch:dword,
;	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; width8 = r9d

Resize_V_AVX_Planar_8bits_ASM proc public frame

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
	sub rsp,48
	.allocstack 48
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	.endprolog		

	mov rsi,valmin
	vbroadcastss xmm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss xmm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss xmm7,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,8

Resize_V_AVX_Planar_8bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa xmm0,xmm7 ;rounder
	vmovdqa xmm1,xmm7

Resize_V_AVX_Planar_8bits_loop_2:
	vpmovzxbw xmm2,qword ptr[rsi]
	vpmovzxbw xmm4,qword ptr[rsi+rbx]

	vbroadcastss xmm8,dword ptr[r8+rax]

	vpunpckhwd xmm3,xmm2,xmm4
	vpunpcklwd xmm2,xmm2,xmm4

	vpmaddwd xmm3,xmm3,xmm8
	vpmaddwd xmm2,xmm2,xmm8

	vpaddd xmm1,xmm1,xmm3
	vpaddd xmm0,xmm0,xmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX_Planar_8bits_loop_2

	vpsrad xmm0,xmm0,14 ;FPScale8bits = 14
	vpsrad xmm1,xmm1,14

	vpackusdw xmm0,xmm0,xmm1

	vpackuswb xmm0,xmm0,xmm0
	vpmaxub xmm0,xmm0,xmm5
	vpminub xmm0,xmm0,xmm6

	vmovq qword ptr[rdi],xmm0
		
	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX_Planar_8bits_loop_1

	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,48

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX_Planar_8bits_ASM endp


;Resize_V_AVX_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,src_pitch:dword,
;	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; width8 = r9d

Resize_V_AVX_Planar_10to14bits_ASM proc public frame

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
	sub rsp,48
	.allocstack 48
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	.endprolog		

	mov rsi,valmin
	vbroadcastss xmm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss xmm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss xmm7,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,16

Resize_V_AVX_Planar_10to14bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa xmm0,xmm7 ;rounder
	vmovdqa xmm1,xmm7

Resize_V_AVX_Planar_10to14bits_loop_2:
	vmovdqa xmm2,XMMWORD ptr[rsi]
	vmovdqa xmm4,XMMWORD ptr[rsi+rbx]

	vbroadcastss xmm8,dword ptr[r8+rax]

	vpunpckhwd xmm3,xmm2,xmm4
	vpunpcklwd xmm2,xmm2,xmm4

	vpmaddwd xmm3,xmm3,xmm8
	vpmaddwd xmm2,xmm2,xmm8

	vpaddd xmm1,xmm1,xmm3
	vpaddd xmm0,xmm0,xmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX_Planar_10to14bits_loop_2

	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	vpsrad xmm1,xmm1,13

	vpackusdw xmm0,xmm0,xmm1

	vpmaxuw xmm0,xmm0,xmm5
	vpminuw xmm0,xmm0,xmm6

	vmovdqa XMMWORD ptr[rdi],xmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX_Planar_10to14bits_loop_1

	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,48

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX_Planar_10to14bits_ASM endp


;Resize_V_AVX_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,width8:dword,
;	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword
;	shifttosigned:dword,shiftfromsigned:dword

; src = rcx
; dst = rdx
; coeff = r8
; width8 = r9d

Resize_V_AVX_Planar_16bits_ASM proc public frame

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
	sub rsp,80
	.allocstack 80
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	vmovdqa XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	vmovdqa XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	.endprolog		

	mov rsi,valmin
	vbroadcastss xmm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss xmm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss xmm7,dword ptr[rsi]
	mov rsi,shifttosigned
	vbroadcastss xmm8,dword ptr[rsi]
	mov rsi,shiftfromsigned
	vbroadcastss xmm9,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,16

Resize_V_AVX_Planar_16bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa xmm0,xmm7 ;rounder
	vmovdqa xmm1,xmm7

Resize_V_AVX_Planar_16bits_loop_2:
	vmovdqa xmm2,XMMWORD ptr[rsi]
	vmovdqa xmm4,XMMWORD ptr[rsi+rbx]

	vbroadcastss xmm10,dword ptr[r8+rax]

	vpunpckhwd xmm3,xmm2,xmm4
	vpunpcklwd xmm2,xmm2,xmm4

	vpaddw xmm3,xmm3,xmm8 ;shifttosigned
	vpaddw xmm2,xmm2,xmm8

	vpmaddwd xmm3,xmm3,xmm10
	vpmaddwd xmm2,xmm2,xmm10

	vpaddd xmm1,xmm1,xmm3
	vpaddd xmm0,xmm0,xmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX_Planar_16bits_loop_2

	vpaddw xmm0,xmm0,xmm9 ;shiftfromsigned
	vpaddw xmm1,xmm1,xmm9

	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	vpsrad xmm1,xmm1,13

	vpackusdw xmm0,xmm0,xmm1

	vpmaxuw xmm0,xmm0,xmm5
	vpminuw xmm0,xmm0,xmm6

	vmovdqa XMMWORD ptr[rdi],xmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX_Planar_16bits_loop_1

	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]	
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]	
	add rsp,80

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX_Planar_16bits_ASM endp


;Resize_V_AVX_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,src_pitch:dword,
;	kernel_size_2:dword

; src = rcx
; dst = rdx
; coeff = r8
; width16 = r9d

Resize_V_AVX_Planar_32bits_ASM proc public frame

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
	mov r12,32

Resize_V_AVX_Planar_32bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1

Resize_V_AVX_Planar_32bits_loop_2:
	vbroadcastss ymm4,dword ptr[r8+rax]
	vbroadcastss ymm5,dword ptr[r8+rax+4]
	
	vmulps ymm6,ymm4,YMMWORD ptr[rsi]
	vmulps ymm7,ymm5,YMMWORD ptr[rsi+rbx]
	vaddps ymm0,ymm0,ymm6
	vaddps ymm1,ymm1,ymm7

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX_Planar_32bits_loop_2

	vaddps ymm0,ymm0,ymm1

	vmovaps YMMWORD ptr[rdi],ymm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz short Resize_V_AVX_Planar_32bits_loop_1

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX_Planar_32bits_ASM endp


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Resize_H_AVX_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_16:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX_Planar_8bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_16 equ dword ptr[rbp+56]
sizeh equ dword ptr[rbp+64]
valmin equ qword ptr[rbp+72]
valmax equ qword ptr[rbp+80]
rounder equ qword ptr[rbp+88]

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
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	push r15
	.pushreg r15
	sub rsp,152
	.allocstack 152
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	vmovdqa XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	vmovdqa XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	vmovdqa XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	vmovdqa XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	vmovdqa XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	vmovdqa XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	.endprolog

	mov rsi,valmin
	vbroadcastss xmm9,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss xmm10,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss xmm11,dword ptr[rsi]
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_16 ;kernel_size_16 = (kernel_size + 7) >> 3
	mov r13d,sizeh
	mov r14,8
	mov r15,16
	shr r13d,2 ;r13d = sizeh/4
	jz Resize_H_AVX_Planar_8bits_1
	
Resize_H_AVX_Planar_8bits_loop_1:
	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1
	vpxor xmm2,xmm2,xmm2
	vpxor xmm3,xmm3,xmm3
	
	shr rcx,1
	jz Resize_H_AVX_Planar_8bits_loop_2b

Resize_H_AVX_Planar_8bits_loop_2:
	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff
	vmovdqa xmm12,XMMWORD ptr[rdi+r15]	;coeff+8(x2)

	vpmovzxbw xmm4,qword ptr[rsi]		;src
	vpmovzxbw xmm13,qword ptr[rsi+r14]	;src+8(x1)
	vpmovzxbw xmm5,qword ptr[r9]		;src+src_pitch
	vpmovzxbw xmm14,qword ptr[r9+r14]	;src+src_pitch+8(x1)
	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm13,xmm13,xmm12
	vpmovzxbw xmm6,qword ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddd xmm4,xmm4,xmm13
	add rsi,r14
	vpmaddwd xmm14,xmm14,xmm12
	vpmovzxbw xmm13,qword ptr[rsi+2*rbx]	;src+2*src_pitch+8(x1)
	vpmaddwd xmm5,xmm5,xmm8
	vpmovzxbw xmm7,qword ptr[r9+2*rbx]	;src+3*src_pitch
	vpmaddwd xmm13,xmm13,xmm12
	vpaddd xmm5,xmm5,xmm14
	add r9,r14
	vpmaddwd xmm6,xmm6,xmm8
	vpmovzxbw xmm14,qword ptr[r9+2*rbx]	;src+3*src_pitch+8(x1)
	vpmaddwd xmm7,xmm7,xmm8
	vpmaddwd xmm14,xmm14,xmm12

	vpaddd xmm6,xmm6,xmm13
	vpaddd xmm7,xmm7,xmm14

	add rdi,r15

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5
	vpaddd xmm2,xmm2,xmm6
	vpaddd xmm3,xmm3,xmm7

	add rsi,r14
	add r9,r14
	add rdi,r15
	dec ecx
	jnz Resize_H_AVX_Planar_8bits_loop_2

Resize_H_AVX_Planar_8bits_loop_2b:
	test r12d,1
	jz short Resize_H_AVX_Planar_8bits_loop_2c

	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff

	vpmovzxbw xmm4,qword ptr[rsi] 		;src
	vpmovzxbw xmm5,qword ptr[r9]		;src+src_pitch
	vpmovzxbw xmm6,qword ptr[rsi+2*rbx]	;src+2*src_pitch
	vpmovzxbw xmm7,qword ptr[r9+2*rbx]	;src+3*src_pitch

	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm5,xmm5,xmm8
	vpmaddwd xmm6,xmm6,xmm8
	vpmaddwd xmm7,xmm7,xmm8

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5
	vpaddd xmm2,xmm2,xmm6
	vpaddd xmm3,xmm3,xmm7

Resize_H_AVX_Planar_8bits_loop_2c:
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	vphaddd xmm0,xmm0,xmm2
	
	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,14 ;FPScale8bits = 14
	
	mov rdx,rbx
	
	vpackusdw xmm0,xmm0,xmm0
	
	sal rdx,2
	
	vpackuswb xmm0,xmm0,xmm0
	
	add r10,rdx
	
	vpmaxub xmm0,xmm0,xmm9

	mov rdx,dst_pitch

	vpminub xmm0,xmm0,xmm10
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,1
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,2
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,3
	mov byte ptr[r11],al
	add r11,rdx

	dec r13d
	jnz Resize_H_AVX_Planar_8bits_loop_1

Resize_H_AVX_Planar_8bits_1:
	test sizeh,2
	jz Resize_H_AVX_Planar_8bits_2

	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1

	shr rcx,1
	jz short Resize_H_AVX_Planar_8bits_loop_4b

Resize_H_AVX_Planar_8bits_loop_4:
	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff
	vmovdqa xmm12,XMMWORD ptr[rdi+r15]	;coeff+8(x2)

	vpmovzxbw xmm4,qword ptr[rsi] 		;src
	vpmovzxbw xmm13,qword ptr[rsi+r14]	;src+8(x1)
	vpmovzxbw xmm5,qword ptr[rsi+rbx]	;src+src_pitch
	vpmaddwd xmm4,xmm4,xmm8
	add rsi,r14
	vpmaddwd xmm13,xmm13,xmm12
	vpmovzxbw xmm14,qword ptr[rsi+rbx]	;src+src_pitch+8(x1)

	vpmaddwd xmm5,xmm5,xmm8
	vpmaddwd xmm14,xmm14,xmm12

	vpaddd xmm4,xmm4,xmm13
	vpaddd xmm5,xmm5,xmm14

	add rdi,r15

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5

	add rsi,r14
	add rdi,r15
	loop Resize_H_AVX_Planar_8bits_loop_4

Resize_H_AVX_Planar_8bits_loop_4b:
	test r12d,1
	jz short Resize_H_AVX_Planar_8bits_loop_4c

	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff

	vpmovzxbw xmm4,qword ptr[rsi] 		;src
	vpmovzxbw xmm5,qword ptr[rsi+rbx]	;src+src_pitch

	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm5,xmm5,xmm8

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5

Resize_H_AVX_Planar_8bits_loop_4c:
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	vphaddd xmm0,xmm0,xmm1

	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,14 ;FPScale18bits = 14
	
	mov rdx,rbx
	
	vpackusdw xmm0,xmm0,xmm0
	
	sal rdx,1
	
	vpackuswb xmm0,xmm0,xmm0
	
	vpmaxub xmm0,xmm0,xmm9
	
	add r10,rdx
	
	vpminub xmm0,xmm0,xmm10

	mov rdx,dst_pitch
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,2
	mov byte ptr[r11],al
	add r11,rdx

Resize_H_AVX_Planar_8bits_2:
	test sizeh,1
	jz Resize_H_AVX_Planar_8bits_end

	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor xmm0,xmm0,xmm0

	shr rcx,1
	jz short Resize_H_AVX_Planar_8bits_loop_5b

Resize_H_AVX_Planar_8bits_loop_5:
	vpmovzxbw xmm4,qword ptr[rsi] 		;src
	vpmovzxbw xmm13,qword ptr[rsi+r14]	;src+8(x1)

	vpmaddwd xmm4,xmm4,XMMWORD ptr[rdi]			;coeff
	vpmaddwd xmm13,xmm13,XMMWORD ptr[rdi+r15]	;coeff+8(x2)

	add rsi,r14

	vpaddd xmm4,xmm4,xmm13

	add rdi,r15

	vpaddd xmm0,xmm0,xmm4

	add rsi,r14
	add rdi,r15
	loop Resize_H_AVX_Planar_8bits_loop_5

Resize_H_AVX_Planar_8bits_loop_5b:
	test r12d,1
	jz short Resize_H_AVX_Planar_8bits_loop_5c

	vpmovzxbw xmm4,qword ptr[rsi] 		;src

	vpmaddwd xmm4,xmm4,XMMWORD ptr[rdi]	;coeff

	vpaddd xmm0,xmm0,xmm4

Resize_H_AVX_Planar_8bits_loop_5c:
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0
	
	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,14 ;FPScale8bits = 14
	
	vpackusdw xmm0,xmm0,xmm0
	vpackuswb xmm0,xmm0,xmm0
	
	vpmaxub xmm0,xmm0,xmm9
	vpminub xmm0,xmm0,xmm10
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al

Resize_H_AVX_Planar_8bits_end:
	vmovdqa xmm14,XMMWORD ptr[rsp+128]
	vmovdqa xmm13,XMMWORD ptr[rsp+112]
	vmovdqa xmm12,XMMWORD ptr[rsp+96]
	vmovdqa xmm11,XMMWORD ptr[rsp+80]
	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,152
		
	vzeroupper

	pop r15
	pop r14
	pop r13
	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_H_AVX_Planar_8bits_ASM endp


;Resize_H_AVX_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_16:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX_Planar_10to14bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_16 equ dword ptr[rbp+56]
sizeh equ dword ptr[rbp+64]
valmin equ qword ptr[rbp+72]
valmax equ qword ptr[rbp+80]
rounder equ qword ptr[rbp+88]

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
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	push r15
	.pushreg r15
	sub rsp,152
	.allocstack 152
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	vmovdqa XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	vmovdqa XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	vmovdqa XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	vmovdqa XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	vmovdqa XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	vmovdqa XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	.endprolog

	mov rsi,valmin
	vbroadcastss xmm9,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss xmm10,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss xmm11,dword ptr[rsi]
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_16 ;kernel_size_16 = (kernel_size + 7) >> 3
	mov r13d,sizeh
	mov r14,16
	mov r15,dst_pitch
	shr r13d,2 ;r13d = sizeh/4
	jz Resize_H_AVX_Planar_10to14bits_1
	
Resize_H_AVX_Planar_10to14bits_loop_1:
	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1
	vpxor xmm2,xmm2,xmm2
	vpxor xmm3,xmm3,xmm3
	
	shr rcx,1
	jz short Resize_H_AVX_Planar_10to14bits_loop_2b

Resize_H_AVX_Planar_10to14bits_loop_2:
	vmovdqa xmm8,XMMWORD ptr[rdi]				;coef
	vmovdqa xmm12,XMMWORD ptr[rdi+r14]			;coef+8(x2)

	vpmaddwd xmm4,xmm8,XMMWORD ptr[rsi] 		;src
	vpmaddwd xmm13,xmm12,XMMWORD ptr[rsi+r14]	;src+8(x2)
	vpmaddwd xmm5,xmm8,XMMWORD ptr[r9]			;src+src_pitch
	vpmaddwd xmm14,xmm12,XMMWORD ptr[r9+r14]		;src+src_pitch+8(x2)

	vpmaddwd xmm6,xmm8,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddd xmm4,xmm4,xmm13
	add rsi,r14
	vpaddd xmm5,xmm5,xmm14
	vpmaddwd xmm13,xmm12,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+8(x2)
	vpmaddwd xmm7,xmm8,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch
	vpaddd xmm6,xmm6,xmm13
	add r9,r14
	vpaddd xmm0,xmm0,xmm4
	vpmaddwd xmm14,xmm12,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+8(x2)
	vpaddd xmm1,xmm1,xmm5
	vpaddd xmm7,xmm7,xmm14

	vpaddd xmm2,xmm2,xmm6
	vpaddd xmm3,xmm3,xmm7

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_10to14bits_loop_2

Resize_H_AVX_Planar_10to14bits_loop_2b:
	test r12d,1
	jz short Resize_H_AVX_Planar_10to14bits_loop_2c

	vmovdqa xmm8,XMMWORD ptr[rdi]				;coef

	vpmaddwd xmm4,xmm8,XMMWORD ptr[rsi] 		;src
	vpmaddwd xmm5,xmm8,XMMWORD ptr[r9]			;src+src_pitch
	vpmaddwd xmm6,xmm8,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpmaddwd xmm7,xmm8,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5
	vpaddd xmm2,xmm2,xmm6
	vpaddd xmm3,xmm3,xmm7

Resize_H_AVX_Planar_10to14bits_loop_2c:
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	vphaddd xmm0,xmm0,xmm2
	
	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw xmm0,xmm0,xmm0

	sal rdx,2

	vpmaxuw xmm0,xmm0,xmm9

	add r10,rdx

	vpminuw xmm0,xmm0,xmm10
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,3
	mov word ptr[r11],ax
	add r11,r15
	
	dec r13d
	jnz Resize_H_AVX_Planar_10to14bits_loop_1

Resize_H_AVX_Planar_10to14bits_1:
	test sizeh,2
	jz Resize_H_AVX_Planar_10to14bits_2

	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1

	shr rcx,1
	jz short Resize_H_AVX_Planar_10to14bits_loop_4b

Resize_H_AVX_Planar_10to14bits_loop_4:
	vmovdqa xmm8,XMMWORD ptr[rdi]			;coeff
	vmovdqa xmm12,XMMWORD ptr[rdi+r14]		;coeff

	vpmaddwd xmm4,xmm8,XMMWORD ptr[rsi] 		;src
	vpmaddwd xmm13,xmm12,XMMWORD ptr[rsi+r14]	;src+8(x2)
	vpmaddwd xmm5,xmm8,XMMWORD ptr[rsi+rbx]		;src+src_pitch
	add rsi,r14
	vpaddd xmm4,xmm4,xmm13
	vpmaddwd xmm14,xmm12,XMMWORD ptr[rsi+rbx]	;src+src_pitch+8(x2)

	add rdi,r14

	vpaddd xmm5,xmm5,xmm14

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_10to14bits_loop_4

Resize_H_AVX_Planar_10to14bits_loop_4b:
	test r12d,1
	jz short Resize_H_AVX_Planar_10to14bits_loop_4c

	vmovdqa xmm8,XMMWORD ptr[rdi]			;coeff

	vpmaddwd xmm4,xmm8,XMMWORD ptr[rsi] 		;src
	vpmaddwd xmm5,xmm8,XMMWORD ptr[rsi+rbx]	;src+src_pitch

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5

Resize_H_AVX_Planar_10to14bits_loop_4c:
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	vphaddd xmm0,xmm0,xmm1
	
	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw xmm0,xmm0,xmm0

	sal rdx,1

	vpmaxuw xmm0,xmm0,xmm9

	add r10,rdx

	vpminuw xmm0,xmm0,xmm10
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15

Resize_H_AVX_Planar_10to14bits_2:
	test sizeh,1
	jz Resize_H_AVX_Planar_10to14bits_end

	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor xmm0,xmm0,xmm0

	shr rcx,1
	jz short Resize_H_AVX_Planar_10to14bits_loop_5b

Resize_H_AVX_Planar_10to14bits_loop_5:
	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff
	vmovdqa xmm12,XMMWORD ptr[rdi+r14]	;coeff+8(x2)

	vpmaddwd xmm4,xmm8,XMMWORD ptr[rsi]			;src
	vpmaddwd xmm13,xmm12,XMMWORD ptr[rsi+r14]	;src+8(x2)

	add rdi,r14

	vpaddd xmm4,xmm4,xmm13

	add rsi,r14

	vpaddd xmm0,xmm0,xmm4

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_10to14bits_loop_5

Resize_H_AVX_Planar_10to14bits_loop_5b:
	test r12d,1
	jz short Resize_H_AVX_Planar_10to14bits_loop_5c

	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff

	vpmaddwd xmm4,xmm8,XMMWORD ptr[rsi] 	;src

	vpaddd xmm0,xmm0,xmm4

Resize_H_AVX_Planar_10to14bits_loop_5c:
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0
	
	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	
	vpackusdw xmm0,xmm0,xmm0

	vpmaxuw xmm0,xmm0,xmm9
	vpminuw xmm0,xmm0,xmm10
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax

Resize_H_AVX_Planar_10to14bits_end:
	vmovdqa xmm14,XMMWORD ptr[rsp+128]
	vmovdqa xmm13,XMMWORD ptr[rsp+112]
	vmovdqa xmm12,XMMWORD ptr[rsp+96]
	vmovdqa xmm11,XMMWORD ptr[rsp+80]
	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,152
		
	vzeroupper

	pop r15
	pop r14
	pop r13
	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_H_AVX_Planar_10to14bits_ASM endp


;Resize_H_AVX_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_16:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword
;	shifttosigned:dword,shiftfromsigned:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX_Planar_16bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_16 equ dword ptr[rbp+56]
sizeh equ dword ptr[rbp+64]
valmin equ qword ptr[rbp+72]
valmax equ qword ptr[rbp+80]
rounder equ qword ptr[rbp+88]
shifttosigned equ qword ptr[rbp+96]
shiftfromsigned equ qword ptr[rbp+104]

Xshiftfromsigned equ XMMWORD ptr[rsp+160]

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
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	push r15
	.pushreg r15
	sub rsp,184
	.allocstack 184
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	vmovdqa XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	vmovdqa XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	vmovdqa XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	vmovdqa XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	vmovdqa XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	vmovdqa XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	vmovdqa XMMWORD ptr[rsp+144],xmm15
	.savexmm128 xmm15,144
	.endprolog

	mov rsi,valmin
	vbroadcastss xmm9,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss xmm10,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss xmm11,dword ptr[rsi]
	mov rsi,shifttosigned
	vbroadcastss xmm12,dword ptr[rsi]
	mov rsi,shiftfromsigned
	vbroadcastss xmm13,dword ptr[rsi]
	vmovdqa Xshiftfromsigned,xmm13

	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_16 ;kernel_size_16 = (kernel_size + 7) >> 3
	mov r13d,sizeh
	mov r14,16
	mov r15,dst_pitch
	shr r13d,2 ;r13d = sizeh/4
	jz Resize_H_AVX_Planar_16bits_1
	
Resize_H_AVX_Planar_16bits_loop_1:
	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1
	vpxor xmm2,xmm2,xmm2
	vpxor xmm3,xmm3,xmm3

	shr rcx,1
	jz Resize_H_AVX_Planar_16bits_loop_2b

Resize_H_AVX_Planar_16bits_loop_2:
	vmovdqa xmm8,XMMWORD ptr[rdi]				;coeff
	vmovdqa xmm13,XMMWORD ptr[rdi+r14]			;coeff+8(x2)

	; shifttosigned + src
	vpaddw xmm4,xmm12,XMMWORD ptr[rsi] 			;src
	vpaddw xmm14,xmm12,XMMWORD ptr[rsi+r14]		;src+8(x2)
	vpaddw xmm5,xmm12,XMMWORD ptr[r9]			;src+src_pitch
	vpaddw xmm15,xmm12,XMMWORD ptr[r9+r14]		;src+src_pitch+8(x2)
	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm14,xmm14,xmm13
	vpmaddwd xmm5,xmm5,xmm8
	vpmaddwd xmm15,xmm15,xmm13

	vpaddw xmm6,xmm12,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddd xmm4,xmm4,xmm14
	add rsi,r14
	vpaddd xmm0,xmm0,xmm4
	vpaddw xmm14,xmm12,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+8(x2)
	vpaddw xmm7,xmm12,XMMWORD ptr[r9+2*rbx]		;src+3*src_pitch
	vpaddd xmm5,xmm5,xmm15
	add r9,r14
	vpaddd xmm1,xmm1,xmm5
	vpaddw xmm15,xmm12,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+8(x2)
	vpmaddwd xmm6,xmm6,xmm8
	vpmaddwd xmm14,xmm14,xmm13
	vpmaddwd xmm7,xmm7,xmm8
	vpmaddwd xmm15,xmm15,xmm13

	vpaddd xmm6,xmm6,xmm14
	vpaddd xmm7,xmm7,xmm15

	vpaddd xmm2,xmm2,xmm6
	vpaddd xmm3,xmm3,xmm7

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rdi,r14
	dec ecx
	jnz Resize_H_AVX_Planar_16bits_loop_2

Resize_H_AVX_Planar_16bits_loop_2b:
	test r12d,1
	jz short Resize_H_AVX_Planar_16bits_loop_2c

	vmovdqa xmm8,XMMWORD ptr[rdi]				;coeff

	; shifttosigned + src
	vpaddw xmm4,xmm12,XMMWORD ptr[rsi] 			;src
	vpaddw xmm5,xmm12,XMMWORD ptr[r9]			;src+src_pitch
	vpaddw xmm6,xmm12,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddw xmm7,xmm12,XMMWORD ptr[r9+2*rbx]		;src+3*src_pitch

	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm5,xmm5,xmm8
	vpmaddwd xmm6,xmm6,xmm8
	vpmaddwd xmm7,xmm7,xmm8

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5
	vpaddd xmm2,xmm2,xmm6
	vpaddd xmm3,xmm3,xmm7

Resize_H_AVX_Planar_16bits_loop_2c:
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	vphaddd xmm0,xmm0,xmm2

	vpaddd xmm0,xmm0,Xshiftfromsigned ;ShiftfromSigned

	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw xmm0,xmm0,xmm0

	sal rdx,2

	vpmaxuw xmm0,xmm0,xmm9

	add r10,rdx

	vpminuw xmm0,xmm0,xmm10

	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,3
	mov word ptr[r11],ax
	add r11,r15
	
	dec r13d
	jnz Resize_H_AVX_Planar_16bits_loop_1

Resize_H_AVX_Planar_16bits_1:
	test sizeh,2
	jz Resize_H_AVX_Planar_16bits_2

	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor xmm0,xmm0,xmm0
	vpxor xmm1,xmm1,xmm1

	shr rcx,1
	jz short Resize_H_AVX_Planar_16bits_loop_4b

Resize_H_AVX_Planar_16bits_loop_4:
	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff
	vmovdqa xmm13,XMMWORD ptr[rdi+r14]	;coeff+8(x2)

	; shifttosigned + src
	vpaddw xmm4,xmm12,XMMWORD ptr[rsi]		;src
	vpaddw xmm14,xmm12,XMMWORD ptr[rsi+r14]	;src+8(x2)
	vpaddw xmm5,xmm12,XMMWORD ptr[r9]		;src+src_pitch
	vpaddw xmm15,xmm12,XMMWORD ptr[r9+r14]	;src+src_pitch+8(x2)

	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm14,xmm14,xmm13
	vpmaddwd xmm5,xmm5,xmm8
	vpmaddwd xmm15,xmm15,xmm13

	add rsi,r14
	add r9,r14
	add rdi,r14

	vpaddd xmm4,xmm4,xmm14
	vpaddd xmm5,xmm5,xmm15

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_16bits_loop_4

Resize_H_AVX_Planar_16bits_loop_4b:
	test r12d,1
	jz short Resize_H_AVX_Planar_16bits_loop_4c

	vmovdqa xmm8,XMMWORD ptr[rdi]		;coeff

	; shifttosigned + src
	vpaddw xmm4,xmm12,XMMWORD ptr[rsi]	;src
	vpaddw xmm5,xmm12,XMMWORD ptr[r9]	;src+src_pitch

	vpmaddwd xmm4,xmm4,xmm8
	vpmaddwd xmm5,xmm5,xmm8

	vpaddd xmm0,xmm0,xmm4
	vpaddd xmm1,xmm1,xmm5

Resize_H_AVX_Planar_16bits_loop_4c:
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	vphaddd xmm0,xmm0,xmm1

	vpaddd xmm0,xmm0,Xshiftfromsigned ;ShiftfromSigned

	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw xmm0,xmm0,xmm0

	sal rdx,1
	
	vpmaxuw xmm0,xmm0,xmm9

	add r10,rdx

	vpminuw xmm0,xmm0,xmm10
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15

Resize_H_AVX_Planar_16bits_2:
	test sizeh,1
	jz Resize_H_AVX_Planar_16bits_end

	mov rcx,r12 ;kernel_size_16
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor xmm0,xmm0,xmm0

	shr rcx,1
	jz short Resize_H_AVX_Planar_16bits_loop_5b

Resize_H_AVX_Planar_16bits_loop_5:
	; shifttosigned + src
	vpaddw xmm4,xmm12,XMMWORD ptr[rsi]		;src
	vpaddw xmm14,xmm12,XMMWORD ptr[rsi+r14]	;src+8(x2)

	vpmaddwd xmm4,xmm4,XMMWORD ptr[rdi]			;coeff
	vpmaddwd xmm14,xmm14,XMMWORD ptr[rdi+r14]	;coeff+8(x2)

	add rsi,r14

	vpaddd xmm4,xmm4,xmm14

	add rdi,r14

	vpaddd xmm0,xmm0,xmm4

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_16bits_loop_5

Resize_H_AVX_Planar_16bits_loop_5b:
	test r12d,1
	jz short Resize_H_AVX_Planar_16bits_loop_5c

	; shifttosigned + src
	vpaddw xmm4,xmm12,XMMWORD ptr[rsi]	;src

	vpmaddwd xmm4,xmm4,XMMWORD ptr[rdi]	;coeff

	vpaddd xmm0,xmm0,xmm4

Resize_H_AVX_Planar_16bits_loop_5c:
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0

	vpaddd xmm0,xmm0,Xshiftfromsigned ;ShiftfromSigned

	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	
	mov rdx,dst_pitch
	
	vpackusdw xmm0,xmm0,xmm0
	
	vpmaxuw xmm0,xmm0,xmm9
	vpminuw xmm0,xmm0,xmm10
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax

Resize_H_AVX_Planar_16bits_end:
	vmovdqa xmm15,XMMWORD ptr[rsp+144]
	vmovdqa xmm14,XMMWORD ptr[rsp+128]
	vmovdqa xmm13,XMMWORD ptr[rsp+112]
	vmovdqa xmm12,XMMWORD ptr[rsp+96]
	vmovdqa xmm11,XMMWORD ptr[rsp+80]
	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,184
		
	vzeroupper

	pop r15
	pop r14
	pop r13
	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_H_AVX_Planar_16bits_ASM endp


;Resize_H_AVX_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_32:dword,sizeh:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX_Planar_32bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_32 equ dword ptr[rbp+56]
sizeh equ dword ptr[rbp+64]

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
	push r13
	.pushreg r13
	push r14
	.pushreg r14
	push r15
	.pushreg r15
	sub rsp,168
	.allocstack 168
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	vmovdqa XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	vmovdqa XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	vmovdqa XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	vmovdqa XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	vmovdqa XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	vmovdqa XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	vmovdqa XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	vmovdqa XMMWORD ptr[rsp+144],xmm15
	.savexmm128 xmm15,144
	.endprolog
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_32 ;kernel_size_32 = (kernel_size + 7) >> 3
	mov r13d,sizeh
	mov r14,32
	mov r15,dst_pitch
	shr r13d,3 ;r13d = sizeh/8
	jz Resize_H_AVX_Planar_32bits_1
	
Resize_H_AVX_Planar_32bits_loop_1:
	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch
	mov rax,r9
	add rax,rbx ;rax=src+2*pitch
	mov rdx,rax
	add rdx,rbx ;rdx=src+3*pitch

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1
	vxorps ymm2,ymm2,ymm2
	vxorps ymm3,ymm3,ymm3
	vxorps ymm4,ymm4,ymm4
	vxorps ymm5,ymm5,ymm5
	vxorps ymm6,ymm6,ymm6
	vxorps ymm7,ymm7,ymm7

	shr rcx,1
	jz Resize_H_AVX_Planar_32bits_loop_2b

Resize_H_AVX_Planar_32bits_loop_2:
	vmovaps ymm12,YMMWORD ptr[rdi]				;coef
	vmovaps ymm13,YMMWORD ptr[rdi+r14]			;coef+8(x4)
	
	vmulps ymm8,ymm12,YMMWORD ptr[rsi] 			;src
	vmulps ymm14,ymm13,YMMWORD ptr[rsi+r14]		;src+8(x4)
	vmulps ymm9,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vmulps ymm15,ymm13,YMMWORD ptr[r9+r14]		;src+src_pitch+8(x4)
	vaddps ymm8,ymm8,ymm14
	vaddps ymm9,ymm9,ymm15
	vmulps ymm10,ymm12,YMMWORD ptr[rax]			;src+2*src_pitch
	vmulps ymm14,ymm13,YMMWORD ptr[rax+r14]		;src+2*src_pitch+8(x4)
	vmulps ymm11,ymm12,YMMWORD ptr[rdx]			;src+3*src_pitch
	vmulps ymm15,ymm13,YMMWORD ptr[rdx+r14]		;src+3*src_pitch+8(x4)
	vaddps ymm10,ymm10,ymm14
	vaddps ymm11,ymm11,ymm15

	vaddps ymm0,ymm0,ymm8
	vaddps ymm1,ymm1,ymm9
	vaddps ymm2,ymm2,ymm10
	vaddps ymm3,ymm3,ymm11
	
	vmulps ymm8,ymm12,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	add rsi,r14
	vmulps ymm14,ymm13,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch+8(x2)
	vmulps ymm9,ymm12,YMMWORD ptr[r9+4*rbx]		;src+5*src_pitch
	add r9,r14
	vmulps ymm15,ymm13,YMMWORD ptr[r9+4*rbx]		;src+5*src_pitch+8(x2)
	vaddps ymm8,ymm8,ymm14
	vaddps ymm9,ymm9,ymm15
	vmulps ymm10,ymm12,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	add rax,r14
	vmulps ymm14,ymm13,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch+8(x2)
	vmulps ymm11,ymm12,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch
	add rdx,r14
	vmulps ymm15,ymm13,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch
	vaddps ymm10,ymm10,ymm14
	vaddps ymm11,ymm11,ymm15

	vaddps ymm4,ymm4,ymm8
	vaddps ymm5,ymm5,ymm9
	vaddps ymm6,ymm6,ymm10
	vaddps ymm7,ymm7,ymm11

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r14
	dec ecx
	jnz Resize_H_AVX_Planar_32bits_loop_2

Resize_H_AVX_Planar_32bits_loop_2b:
	test r12d,1
	jz short Resize_H_AVX_Planar_32bits_loop_2c

	vmovaps ymm12,YMMWORD ptr[rdi]				;coef
	
	vmulps ymm8,ymm12,YMMWORD ptr[rsi] 			;src
	vmulps ymm9,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vmulps ymm10,ymm12,YMMWORD ptr[rax]			;src+2*src_pitch
	vmulps ymm11,ymm12,YMMWORD ptr[rdx]			;src+3*src_pitch
	vaddps ymm0,ymm0,ymm8
	vaddps ymm1,ymm1,ymm9
	vaddps ymm2,ymm2,ymm10
	vaddps ymm3,ymm3,ymm11
	
	vmulps ymm8,ymm12,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	vmulps ymm9,ymm12,YMMWORD ptr[r9+4*rbx]		;src+5*src_pitch
	vmulps ymm10,ymm12,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	vmulps ymm11,ymm12,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch
	vaddps ymm4,ymm4,ymm8
	vaddps ymm5,ymm5,ymm9
	vaddps ymm6,ymm6,ymm10
	vaddps ymm7,ymm7,ymm11

Resize_H_AVX_Planar_32bits_loop_2c:
	vextractf128 xmm8,ymm0,1
	vextractf128 xmm9,ymm1,1
	vextractf128 xmm10,ymm2,1
	vextractf128 xmm11,ymm3,1
	vhaddps xmm0,xmm0,xmm8
	vhaddps xmm1,xmm1,xmm9
	vhaddps xmm2,xmm2,xmm10
	vhaddps xmm3,xmm3,xmm11
	vextractf128 xmm8,ymm4,1
	vextractf128 xmm9,ymm5,1
	vextractf128 xmm10,ymm6,1
	vextractf128 xmm11,ymm7,1
	vhaddps xmm4,xmm4,xmm8
	vhaddps xmm5,xmm5,xmm9
	vhaddps xmm6,xmm6,xmm10
	vhaddps xmm7,xmm7,xmm11
	
	vhaddps xmm0,xmm0,xmm0
	vhaddps xmm1,xmm1,xmm1
	vhaddps xmm2,xmm2,xmm2
	vhaddps xmm3,xmm3,xmm3
	vhaddps xmm4,xmm4,xmm4
	vhaddps xmm5,xmm5,xmm5
	vhaddps xmm6,xmm6,xmm6
	vhaddps xmm7,xmm7,xmm7
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68
	vshufps xmm4,xmm4,xmm5,68
	vshufps xmm6,xmm6,xmm7,68

	mov rdx,rbx

	vhaddps xmm0,xmm0,xmm2

	sal rdx,3

	vhaddps xmm4,xmm4,xmm6

	add r10,rdx

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,1
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,2
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,3
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm4,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm4,1
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm4,2
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm4,3
	mov dword ptr[r11],eax
	add r11,r15
	
	dec r13d
	jnz Resize_H_AVX_Planar_32bits_loop_1

Resize_H_AVX_Planar_32bits_1:
	test sizeh,4
	jz Resize_H_AVX_Planar_32bits_2

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1
	vxorps ymm2,ymm2,ymm2
	vxorps ymm3,ymm3,ymm3

	shr rcx,1
	jz short Resize_H_AVX_Planar_32bits_loop_3b

Resize_H_AVX_Planar_32bits_loop_3:
	vmovaps ymm12,YMMWORD ptr[rdi]				;coeff
	vmovaps ymm13,YMMWORD ptr[rdi+r14]			;coeff

	vmulps ymm4,ymm12,YMMWORD ptr[rsi] 			;src
	vmulps ymm8,ymm13,YMMWORD ptr[rsi+r14]		;src+8(x4)
	vmulps ymm5,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vmulps ymm9,ymm13,YMMWORD ptr[r9+r14]		;src+src_pitch
	vmulps ymm6,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vaddps ymm4,ymm4,ymm8
	add rsi,r14
	vaddps ymm5,ymm5,ymm9
	vmulps ymm10,ymm13,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+8(x4)
	vmulps ymm7,ymm12,YMMWORD ptr[r9+2*rbx]		;src+3*src_pitch
	vaddps ymm6,ymm6,ymm10
	add r9,r14
	vaddps ymm0,ymm0,ymm4
	vmulps ymm11,ymm13,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+8(x4)
	vaddps ymm1,ymm1,ymm5
	vaddps ymm7,ymm7,ymm11

	vaddps ymm2,ymm2,ymm6
	vaddps ymm3,ymm3,ymm7

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_32bits_loop_3

Resize_H_AVX_Planar_32bits_loop_3b:
	test r12d,1
	jz short Resize_H_AVX_Planar_32bits_loop_3c

	vmovaps ymm12,YMMWORD ptr[rdi]				;coeff

	vmulps ymm4,ymm12,YMMWORD ptr[rsi] 			;src
	vmulps ymm5,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vmulps ymm6,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vmulps ymm7,ymm12,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch
	vaddps ymm0,ymm0,ymm4
	vaddps ymm1,ymm1,ymm5
	vaddps ymm2,ymm2,ymm6
	vaddps ymm3,ymm3,ymm7

Resize_H_AVX_Planar_32bits_loop_3c:
	vextractf128 xmm8,ymm0,1
	vextractf128 xmm9,ymm1,1
	vextractf128 xmm10,ymm2,1
	vextractf128 xmm11,ymm3,1

	vhaddps xmm0,xmm0,xmm8
	vhaddps xmm1,xmm1,xmm9
	vhaddps xmm2,xmm2,xmm10
	vhaddps xmm3,xmm3,xmm11

	vhaddps xmm0,xmm0,xmm0
	vhaddps xmm1,xmm1,xmm1
	vhaddps xmm2,xmm2,xmm2
	vhaddps xmm3,xmm3,xmm3

	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	mov rdx,rbx

	sal rdx,2

	vhaddps xmm0,xmm0,xmm2

	add r10,rdx

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,1
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,2
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,3
	mov dword ptr[r11],eax
	add r11,r15

Resize_H_AVX_Planar_32bits_2:
	test sizeh,2
	jz Resize_H_AVX_Planar_32bits_3

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1

	shr rcx,1
	jz short Resize_H_AVX_Planar_32bits_loop_4b

Resize_H_AVX_Planar_32bits_loop_4:
	vmovaps ymm12,YMMWORD ptr[rdi]				;coeff
	vmovaps ymm13,YMMWORD ptr[rdi+r14]			;coeff+8(x4)

	vmulps ymm4,ymm12,YMMWORD ptr[rsi] 			;src
	vmulps ymm8,ymm13,YMMWORD ptr[rsi+r14]		;src+8(x4)
	vmulps ymm5,ymm12,YMMWORD ptr[rsi+rbx]		;src+src_pitch
	add rsi,r14
	vaddps ymm4,ymm4,ymm8
	vmulps ymm9,ymm13,YMMWORD ptr[rsi+rbx]		;src+src_pitch+8(x4)

	vaddps ymm0,ymm0,ymm4
	vaddps ymm5,ymm5,ymm9

	add rdi,r14

	vaddps ymm1,ymm1,ymm5

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_32bits_loop_4

Resize_H_AVX_Planar_32bits_loop_4b:
	test r12d,1
	jz short Resize_H_AVX_Planar_32bits_loop_4c

	vmovaps ymm12,YMMWORD ptr[rdi]				;coeff

	vmulps ymm8,ymm12,YMMWORD ptr[rsi] 			;src
	vmulps ymm9,ymm12,YMMWORD ptr[rsi+rbx]		;src+src_pitch
	vaddps ymm0,ymm0,ymm8
	vaddps ymm1,ymm1,ymm9

Resize_H_AVX_Planar_32bits_loop_4c:
	vextractf128 xmm8,ymm0,1
	vextractf128 xmm9,ymm1,1

	vhaddps xmm0,xmm0,xmm8
	vhaddps xmm1,xmm1,xmm9

	vhaddps xmm0,xmm0,xmm0

	mov rdx,rbx

	vhaddps xmm1,xmm1,xmm1

	sal rdx,1

	vhaddps xmm0,xmm0,xmm1

	add r10,rdx

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,2
	mov dword ptr[r11],eax
	add r11,r15

Resize_H_AVX_Planar_32bits_3:
	test sizeh,1
	jz short Resize_H_AVX_Planar_32bits_end

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vxorps ymm0,ymm0,ymm0

	shr rcx,1
	jz short Resize_H_AVX_Planar_32bits_loop_5b

Resize_H_AVX_Planar_32bits_loop_5:
	vmovaps ymm12,YMMWORD ptr[rdi]			;coeff
	vmovaps ymm13,YMMWORD ptr[rdi+r14]		;coeff+8(x4)

	vmulps ymm4,ymm12,YMMWORD ptr[rsi] 		;src
	vmulps ymm8,ymm13,YMMWORD ptr[rsi+r14]	;src+8(x4)

	add rsi,r14

	vaddps ymm4,ymm4,ymm8

	add rdi,r14

	vaddps ymm0,ymm0,ymm4

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX_Planar_32bits_loop_5

Resize_H_AVX_Planar_32bits_loop_5b:
	test r12d,1
	jz short Resize_H_AVX_Planar_32bits_loop_5c

	vmovaps ymm12,YMMWORD ptr[rdi]		;coeff

	vmulps ymm4,ymm12,YMMWORD ptr[rsi]	;src
	vaddps ymm0,ymm0,ymm4

Resize_H_AVX_Planar_32bits_loop_5c:
	vextractf128 xmm8,ymm0,1

	vhaddps xmm0,xmm0,xmm8

	vhaddps xmm0,xmm0,xmm0

	vhaddps xmm0,xmm0,xmm0

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax

Resize_H_AVX_Planar_32bits_end:
	vmovdqa xmm15,XMMWORD ptr[rsp+144]
	vmovdqa xmm14,XMMWORD ptr[rsp+128]
	vmovdqa xmm13,XMMWORD ptr[rsp+112]
	vmovdqa xmm12,XMMWORD ptr[rsp+96]
	vmovdqa xmm11,XMMWORD ptr[rsp+80]
	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,168
		
	vzeroupper

	pop r15
	pop r14
	pop r13
	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_H_AVX_Planar_32bits_ASM endp


end
