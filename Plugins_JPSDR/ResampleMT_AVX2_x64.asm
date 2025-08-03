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


;Resize_V_AVX2_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,src_pitch:dword,
;	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; width16 = r9d

Resize_V_AVX2_Planar_8bits_ASM proc public frame

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
	vbroadcastss ymm7,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,16

Resize_V_AVX2_Planar_8bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa ymm0,ymm7 ;rounder
	vmovdqa ymm1,ymm7

Resize_V_AVX2_Planar_8bits_loop_2:
	vpmovzxbw ymm2,XMMWORD ptr[rsi]
	vpmovzxbw ymm4,XMMWORD ptr[rsi+rbx]

	vbroadcastss ymm8,dword ptr[r8+rax]

	vpunpckhwd ymm3,ymm2,ymm4
	vpunpcklwd ymm2,ymm2,ymm4

	vpmaddwd ymm3,ymm3,ymm8
	vpmaddwd ymm2,ymm2,ymm8

	vpaddd ymm1,ymm1,ymm3
	vpaddd ymm0,ymm0,ymm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX2_Planar_8bits_loop_2

	vpsrad ymm0,ymm0,14 ;FPScale8bits = 14
	vpsrad ymm1,ymm1,14

	vpackusdw ymm0,ymm0,ymm1

	vextracti128 xmm1,ymm0,1
	vpackuswb xmm0,xmm0,xmm1
	vpmaxub xmm0,xmm0,xmm5
	vpminub xmm0,xmm0,xmm6

	vmovdqa XMMWORD ptr[rdi],xmm0
		
	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX2_Planar_8bits_loop_1

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

Resize_V_AVX2_Planar_8bits_ASM endp


;Resize_V_AVX2_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,src_pitch:dword,
;	kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; width16 = r9d

Resize_V_AVX2_Planar_10to14bits_ASM proc public frame

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
	vbroadcastss ymm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss ymm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss ymm7,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,32

Resize_V_AVX2_Planar_10to14bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa ymm0,ymm7 ;rounder
	vmovdqa ymm1,ymm7

Resize_V_AVX2_Planar_10to14bits_loop_2:
	vmovdqa ymm2,YMMWORD ptr[rsi]
	vmovdqa ymm4,YMMWORD ptr[rsi+rbx]

	vbroadcastss ymm8,dword ptr[r8+rax]

	vpunpckhwd ymm3,ymm2,ymm4
	vpunpcklwd ymm2,ymm2,ymm4

	vpmaddwd ymm3,ymm3,ymm8
	vpmaddwd ymm2,ymm2,ymm8

	vpaddd ymm1,ymm1,ymm3
	vpaddd ymm0,ymm0,ymm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX2_Planar_10to14bits_loop_2

	vpsrad ymm0,ymm0,13 ;FPScale16bits = 13
	vpsrad ymm1,ymm1,13

	vpackusdw ymm0,ymm0,ymm1

	vpmaxuw ymm0,ymm0,ymm5
	vpminuw ymm0,ymm0,ymm6

	vmovdqa YMMWORD ptr[rdi],ymm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX2_Planar_10to14bits_loop_1

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

Resize_V_AVX2_Planar_10to14bits_ASM endp


;Resize_V_AVX2_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,
;	src_pitch:dword,kernel_size_2:dword,valmin:dword,valmax:dword,rounder:dword
;	shifttosigned:dword,shiftfromsigned:dword

; src = rcx
; dst = rdx
; coeff = r8
; width16 = r9d

Resize_V_AVX2_Planar_16bits_ASM proc public frame

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
	vbroadcastss ymm5,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss ymm6,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss ymm7,dword ptr[rsi]
	mov rsi,shifttosigned
	vbroadcastss ymm8,dword ptr[rsi]
	mov rsi,shiftfromsigned
	vbroadcastss ymm9,dword ptr[rsi]

	mov r10,rcx ; src
	mov rdi,rdx ; dst
	
	xor r11,r11
	mov r11d,kernel_size_2 ;kernel_size_2 = (kernel_size + 1) >> 1
	mov rdx,4 ; 2 coeff of 16 bits
	mov rbx,src_pitch
	mov r12,32

Resize_V_AVX2_Planar_16bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vmovdqa ymm0,ymm7 ;rounder
	vmovdqa ymm1,ymm7

Resize_V_AVX2_Planar_16bits_loop_2:
	vmovdqa ymm2,YMMWORD ptr[rsi]
	vmovdqa ymm4,YMMWORD ptr[rsi+rbx]

	vbroadcastss ymm10,dword ptr[r8+rax]

	vpunpckhwd ymm3,ymm2,ymm4
	vpunpcklwd ymm2,ymm2,ymm4

	vpaddw ymm3,ymm3,ymm8 ;shifttosigned
	vpaddw ymm2,ymm2,ymm8

	vpmaddwd ymm3,ymm3,ymm10
	vpmaddwd ymm2,ymm2,ymm10

	vpaddd ymm1,ymm1,ymm3
	vpaddd ymm0,ymm0,ymm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX2_Planar_16bits_loop_2

	vpaddw ymm0,ymm0,ymm9 ;shiftfromsigned
	vpaddw ymm1,ymm1,ymm9

	vpsrad ymm0,ymm0,13 ;FPScale16bits = 13
	vpsrad ymm1,ymm1,13

	vpackusdw ymm0,ymm0,ymm1

	vpmaxuw ymm0,ymm0,ymm5
	vpminuw ymm0,ymm0,ymm6

	vmovdqa YMMWORD ptr[rdi],ymm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX2_Planar_16bits_loop_1

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

Resize_V_AVX2_Planar_16bits_ASM endp


;Resize_V_AVX2_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,width16:dword,src_pitch:dword,
;	kernel_size_2:dword

; src = rcx
; dst = rdx
; coeff = r8
; width16 = r9d

Resize_V_AVX2_Planar_32bits_ASM proc public frame

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

Resize_V_AVX2_Planar_32bits_loop_1:
	mov rcx,r11 ; rcx = kernel_size_2
	mov rsi,r10 ; src + x
	xor rax,rax

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1

Resize_V_AVX2_Planar_32bits_loop_2:
	vbroadcastss ymm4,dword ptr[r8+rax]
	vbroadcastss ymm5,dword ptr[r8+rax+4]

	vfmadd231ps ymm0,ymm4,YMMWORD ptr[rsi]
	vfmadd231ps ymm1,ymm5,YMMWORD ptr[rsi+rbx]

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX2_Planar_32bits_loop_2

	vaddps ymm0,ymm0,ymm1

	vmovaps YMMWORD ptr[rdi],ymm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz short Resize_V_AVX2_Planar_32bits_loop_1

	vzeroupper

	pop r12
	pop rsi
	pop rdi
	pop rbx
	pop rbp

	ret

Resize_V_AVX2_Planar_32bits_ASM endp


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Resize_H_AVX2_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_32:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX2_Planar_8bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_32 equ dword ptr[rbp+56]
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
	mov r12d,kernel_size_32 ;kernel_size_32 = (kernel_size + 15) >> 4
	mov r13d,sizeh
	mov r14,16
	mov r15,32
	shr r13d,2 ;r13d = sizeh/4
	jz Resize_H_AVX2_Planar_8bits_1
	
Resize_H_AVX2_Planar_8bits_loop_1:
	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor ymm0,ymm0,ymm0
	vpxor ymm1,ymm1,ymm1
	vpxor ymm2,ymm2,ymm2
	vpxor ymm3,ymm3,ymm3
	
	shr rcx,1
	jz Resize_H_AVX2_Planar_8bits_loop_2_b

Resize_H_AVX2_Planar_8bits_loop_2:
	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff
	vmovdqa ymm12,YMMWORD ptr[rdi+r15]	;coeff+16(x2)

	vpmovzxbw ymm4,XMMWORD ptr[rsi]			;src
	vpmovzxbw ymm13,XMMWORD ptr[rsi+r14]	;src+16(x1)
	vpmovzxbw ymm5,XMMWORD ptr[r9]			;src+src_pitch
	vpmovzxbw ymm14,XMMWORD ptr[r9+r14]		;src+src_pitch+16(x1)
	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm13,ymm13,ymm12
	vpmovzxbw ymm6,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddd ymm4,ymm4,ymm13
	add rsi,r14
	vpmaddwd ymm14,ymm14,ymm12
	vpmovzxbw ymm13,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+16(x1)
	vpmaddwd ymm5,ymm5,ymm8
	vpmovzxbw ymm7,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch
	vpmaddwd ymm13,ymm13,ymm12
	vpaddd ymm5,ymm5,ymm14
	add r9,r14
	vpmaddwd ymm6,ymm6,ymm8
	vpmovzxbw ymm14,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+16(x1)
	vpmaddwd ymm7,ymm7,ymm8
	vpmaddwd ymm14,ymm14,ymm12
	
	vpaddd ymm6,ymm6,ymm13
	vpaddd ymm7,ymm7,ymm14

	add rdi,r15

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5
	vpaddd ymm2,ymm2,ymm6
	vpaddd ymm3,ymm3,ymm7

	add rsi,r14
	add r9,r14
	add rdi,r15
	dec rcx
	jnz Resize_H_AVX2_Planar_8bits_loop_2

Resize_H_AVX2_Planar_8bits_loop_2_b:
	test kernel_size_32,1
	jz short Resize_H_AVX2_Planar_8bits_loop_2_c

	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff

	vpmovzxbw ymm4,XMMWORD ptr[rsi]			;src
	vpmovzxbw ymm5,XMMWORD ptr[r9]			;src+src_pitch
	vpmovzxbw ymm6,XMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpmovzxbw ymm7,XMMWORD ptr[r9+2*rbx]	;src+3*src_pitch

	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm5,ymm5,ymm8
	vpmaddwd ymm6,ymm6,ymm8
	vpmaddwd ymm7,ymm7,ymm8

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5
	vpaddd ymm2,ymm2,ymm6
	vpaddd ymm3,ymm3,ymm7
	
Resize_H_AVX2_Planar_8bits_loop_2_c:
	vextracti128 xmm4,ymm0,1
	vextracti128 xmm5,ymm1,1
	vextracti128 xmm6,ymm2,1
	vextracti128 xmm7,ymm3,1

	vphaddd xmm0,xmm0,xmm4
	vphaddd xmm1,xmm1,xmm5
	vphaddd xmm2,xmm2,xmm6
	vphaddd xmm3,xmm3,xmm7
	
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
	jnz Resize_H_AVX2_Planar_8bits_loop_1

Resize_H_AVX2_Planar_8bits_1:
	test sizeh,2
	jz Resize_H_AVX2_Planar_8bits_2

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor ymm0,ymm0,ymm0
	vpxor ymm1,ymm1,ymm1

	shr rcx,1
	jz short Resize_H_AVX2_Planar_8bits_loop_4_b

Resize_H_AVX2_Planar_8bits_loop_4:
	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff
	vmovdqa ymm12,YMMWORD ptr[rdi+r15]	;coeff+16(x2)

	vpmovzxbw ymm4,XMMWORD ptr[rsi] 	;src
	vpmovzxbw ymm13,XMMWORD ptr[rsi+r14]	;src+16(x1)
	vpmovzxbw ymm5,XMMWORD ptr[rsi+rbx]	;src+src_pitch
	vpmaddwd ymm4,ymm4,ymm8
	add rsi,r14
	vpmaddwd ymm5,ymm5,ymm8
	vpmovzxbw ymm14,XMMWORD ptr[rsi+rbx]	;src+src_pitch+16(x1)
	vpmaddwd ymm13,ymm13,ymm12
	vpmaddwd ymm14,ymm14,ymm12

	vpaddd ymm4,ymm4,ymm13
	vpaddd ymm5,ymm5,ymm14

	add rdi,r15

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5

	add rsi,r14
	add rdi,r15
	loop Resize_H_AVX2_Planar_8bits_loop_4
	
Resize_H_AVX2_Planar_8bits_loop_4_b:
	test kernel_size_32,1
	jz short Resize_H_AVX2_Planar_8bits_loop_4_c

	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff

	vpmovzxbw ymm4,XMMWORD ptr[rsi] 	;src
	vpmovzxbw ymm5,XMMWORD ptr[rsi+rbx]	;src+src_pitch

	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm5,ymm5,ymm8

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5

Resize_H_AVX2_Planar_8bits_loop_4_c:
	vextracti128 xmm4,ymm0,1
	vextracti128 xmm5,ymm1,1

	vphaddd xmm0,xmm0,xmm4
	vphaddd xmm1,xmm1,xmm5
	
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

Resize_H_AVX2_Planar_8bits_2:
	test sizeh,1
	jz Resize_H_AVX2_Planar_8bits_end

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor ymm0,ymm0,ymm0

	shr rcx,1
	jz short Resize_H_AVX2_Planar_8bits_loop_5_b

Resize_H_AVX2_Planar_8bits_loop_5:
	vpmovzxbw ymm4,XMMWORD ptr[rsi] 		;src
	vpmovzxbw ymm13,XMMWORD ptr[rsi+r14]	;src+16(x1)

	vpmaddwd ymm4,ymm4,YMMWORD ptr[rdi]			;coeff
	vpmaddwd ymm13,ymm13,YMMWORD ptr[rdi+r15]	;coeff+16(x2)

	add rsi,r14

	vpaddd ymm4,ymm4,ymm13

	add rdi,r15

	vpaddd ymm0,ymm0,ymm4

	add rsi,r14
	add rdi,r15
	loop Resize_H_AVX2_Planar_8bits_loop_5

Resize_H_AVX2_Planar_8bits_loop_5_b:
	test kernel_size_32,1
	jz short Resize_H_AVX2_Planar_8bits_loop_5_c

	vpmovzxbw ymm4,XMMWORD ptr[rsi] 		;src

	vpmaddwd ymm4,ymm4,YMMWORD ptr[rdi]		;coeff

	vpaddd ymm0,ymm0,ymm4

Resize_H_AVX2_Planar_8bits_loop_5_c:
	vextracti128 xmm4,ymm0,1

	vphaddd xmm0,xmm0,xmm4
	
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

Resize_H_AVX2_Planar_8bits_end:
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

Resize_H_AVX2_Planar_8bits_ASM endp


;Resize_H_AVX2_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_32:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX2_Planar_10to14bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_32 equ dword ptr[rbp+56]
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
	mov r12d,kernel_size_32 ;kernel_size_32 = (kernel_size + 15) >> 4
	mov r13d,sizeh
	mov r14,32
	mov r15,dst_pitch
	shr r13d,2 ;r13d = sizeh/4
	jz Resize_H_AVX2_Planar_10to14bits_1
	
Resize_H_AVX2_Planar_10to14bits_loop_1:
	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor ymm0,ymm0,ymm0
	vpxor ymm1,ymm1,ymm1
	vpxor ymm2,ymm2,ymm2
	vpxor ymm3,ymm3,ymm3
	
	shr rcx,1
	jz short Resize_H_AVX2_Planar_10to14bits_loop_2_b

Resize_H_AVX2_Planar_10to14bits_loop_2:
	vmovdqa ymm8,YMMWORD ptr[rdi]			;coef
	vmovdqa ymm12,YMMWORD ptr[rdi+r14]		;coef+16(x2)

	vpmaddwd ymm4,ymm8,YMMWORD ptr[rsi] 		;src
	vpmaddwd ymm13,ymm12,YMMWORD ptr[rsi+r14]	;src+16(x2)
	vpmaddwd ymm5,ymm8,YMMWORD ptr[r9]			;src+src_pitch
	vpmaddwd ymm14,ymm12,YMMWORD ptr[r9+r14]	;src+src_pitch+16(x2)

	vpmaddwd ymm6,ymm8,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddd ymm4,ymm4,ymm13
	add rsi,r14
	vpaddd ymm5,ymm5,ymm14
	vpmaddwd ymm13,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+16(x2)
	vpmaddwd ymm7,ymm8,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch
	vpaddd ymm6,ymm6,ymm13
	add r9,r14
	vpaddd ymm0,ymm0,ymm4
	vpmaddwd ymm14,ymm12,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+16(x2)
	vpaddd ymm1,ymm1,ymm5
	vpaddd ymm7,ymm7,ymm14

	vpaddd ymm2,ymm2,ymm6
	vpaddd ymm3,ymm3,ymm7

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_10to14bits_loop_2

Resize_H_AVX2_Planar_10to14bits_loop_2_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_10to14bits_loop_2_c

	vmovdqa ymm8,YMMWORD ptr[rdi]				;coef

	vpmaddwd ymm4,ymm8,YMMWORD ptr[rsi] 		;src
	vpmaddwd ymm5,ymm8,YMMWORD ptr[r9]			;src+src_pitch
	vpmaddwd ymm6,ymm8,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpmaddwd ymm7,ymm8,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5
	vpaddd ymm2,ymm2,ymm6
	vpaddd ymm3,ymm3,ymm7

Resize_H_AVX2_Planar_10to14bits_loop_2_c:
	vextracti128 xmm4,ymm0,1
	vextracti128 xmm5,ymm1,1
	vextracti128 xmm6,ymm2,1
	vextracti128 xmm7,ymm3,1

	vphaddd xmm0,xmm0,xmm4
	vphaddd xmm1,xmm1,xmm5
	vphaddd xmm2,xmm2,xmm6
	vphaddd xmm3,xmm3,xmm7
	
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
	jnz Resize_H_AVX2_Planar_10to14bits_loop_1

Resize_H_AVX2_Planar_10to14bits_1:
	test sizeh,2
	jz Resize_H_AVX2_Planar_10to14bits_2

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor ymm0,ymm0,ymm0
	vpxor ymm1,ymm1,ymm1

	shr rcx,1
	jz short Resize_H_AVX2_Planar_10to14bits_loop_4_b

Resize_H_AVX2_Planar_10to14bits_loop_4:
	vmovdqa ymm8,YMMWORD ptr[rdi]			;coeff
	vmovdqa ymm12,YMMWORD ptr[rdi+r14]		;coeff+16(x2)

	vpmaddwd ymm4,ymm8,YMMWORD ptr[rsi] 		;src
	vpmaddwd ymm13,ymm12,YMMWORD ptr[rsi+r14]	;src+16(x2)
	vpmaddwd ymm5,ymm8,YMMWORD ptr[rsi+rbx]		;src+src_pitch
	add rsi,r14
	vpaddd ymm4,ymm4,ymm13
	vpmaddwd ymm14,ymm12,YMMWORD ptr[rsi+rbx]	;src+src_pitch+16(x2)

	add rdi,r14

	vpaddd ymm5,ymm5,ymm14

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_10to14bits_loop_4

Resize_H_AVX2_Planar_10to14bits_loop_4_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_10to14bits_loop_4_c

	vmovdqa ymm8,YMMWORD ptr[rdi]			;coeff

	vpmaddwd ymm4,ymm8,YMMWORD ptr[rsi] 		;src
	vpmaddwd ymm5,ymm8,YMMWORD ptr[rsi+rbx]	;src+src_pitch

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5

Resize_H_AVX2_Planar_10to14bits_loop_4_c:
	vextracti128 xmm4,ymm0,1
	vextracti128 xmm5,ymm1,1

	vphaddd xmm0,xmm0,xmm4
	vphaddd xmm1,xmm1,xmm5
	
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

Resize_H_AVX2_Planar_10to14bits_2:
	test sizeh,1
	jz Resize_H_AVX2_Planar_10to14bits_end

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor ymm0,ymm0,ymm0

	shr rcx,1
	jz short Resize_H_AVX2_Planar_10to14bits_loop_5_b

Resize_H_AVX2_Planar_10to14bits_loop_5:
	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff
	vmovdqa ymm12,YMMWORD ptr[rdi+r14]	;coeff+16(x2)

	vpmaddwd ymm4,ymm8,YMMWORD ptr[rsi]			;src
	vpmaddwd ymm13,ymm12,YMMWORD ptr[rsi+r14]	;src+16(x2)

	add rdi,r14

	vpaddd ymm4,ymm4,ymm13

	add rsi,r14

	vpaddd ymm0,ymm0,ymm4

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_10to14bits_loop_5

Resize_H_AVX2_Planar_10to14bits_loop_5_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_10to14bits_loop_5_c

	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff

	vpmaddwd ymm4,ymm8,YMMWORD ptr[rsi] ;src

	vpaddd ymm0,ymm0,ymm4

Resize_H_AVX2_Planar_10to14bits_loop_5_c:

	vextracti128 xmm4,ymm0,1

	vphaddd xmm0,xmm0,xmm4
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0
	
	vpaddd xmm0,xmm0,xmm11 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	
	vpackusdw xmm0,xmm0,xmm0

	vpmaxuw xmm0,xmm0,xmm9
	vpminuw xmm0,xmm0,xmm10
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax

Resize_H_AVX2_Planar_10to14bits_end:
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

Resize_H_AVX2_Planar_10to14bits_ASM endp


;Resize_H_AVX2_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_32:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword
;	shifttosigned:dword,shiftfromsigned:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX2_Planar_16bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_32 equ dword ptr[rbp+56]
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
	vbroadcastss ymm12,dword ptr[rsi]
	mov rsi,shiftfromsigned
	vbroadcastss xmm13,dword ptr[rsi]
	vmovdqa Xshiftfromsigned,xmm13
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_32 ;kernel_size_32 = (kernel_size + 15) >> 4
	mov r13d,sizeh
	mov r14,32
	mov r15,dst_pitch
	shr r13d,2 ;r13d = sizeh/4
	jz Resize_H_AVX2_Planar_16bits_1
	
Resize_H_AVX2_Planar_16bits_loop_1:
	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor ymm0,ymm0,ymm0
	vpxor ymm1,ymm1,ymm1
	vpxor ymm2,ymm2,ymm2
	vpxor ymm3,ymm3,ymm3
	
	shr rcx,1
	jz Resize_H_AVX2_Planar_16bits_loop_2_b

Resize_H_AVX2_Planar_16bits_loop_2:
	vmovdqa ymm8,YMMWORD ptr[rdi]				;coeff
	vmovdqa ymm13,YMMWORD ptr[rdi+r14]			;coeff+16(x2)

	; shifttosigned + src
	vpaddw ymm4,ymm12,YMMWORD ptr[rsi] 			;src
	vpaddw ymm14,ymm12,YMMWORD ptr[rsi+r14]		;src+16(x2)
	vpaddw ymm5,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vpaddw ymm15,ymm12,YMMWORD ptr[r9+r14]		;src+src_pitch+16(x2)
	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm14,ymm14,ymm13
	vpmaddwd ymm5,ymm5,ymm8
	vpmaddwd ymm15,ymm15,ymm13

	vpaddw ymm6,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddd ymm4,ymm4,ymm14
	add rsi,r14
	vpaddd ymm0,ymm0,ymm4
	vpaddw ymm14,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+16(x2)
	vpaddw ymm7,ymm12,YMMWORD ptr[r9+2*rbx]		;src+3*src_pitch
	vpaddd ymm5,ymm5,ymm15
	add r9,r14
	vpaddd ymm1,ymm1,ymm5
	vpaddw ymm15,ymm12,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+16(x2)
	vpmaddwd ymm6,ymm6,ymm8
	vpmaddwd ymm14,ymm14,ymm13
	vpmaddwd ymm7,ymm7,ymm8
	vpmaddwd ymm15,ymm15,ymm13

	vpaddd ymm6,ymm6,ymm14
	vpaddd ymm7,ymm7,ymm15

	vpaddd ymm2,ymm2,ymm6
	vpaddd ymm3,ymm3,ymm7

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rdi,r14
	dec ecx
	jnz Resize_H_AVX2_Planar_16bits_loop_2
	
Resize_H_AVX2_Planar_16bits_loop_2_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_16bits_loop_2_c

	vmovdqa ymm8,YMMWORD ptr[rdi]				;coeff

	; shifttosigned + src
	vpaddw ymm4,ymm12,YMMWORD ptr[rsi] 			;src
	vpaddw ymm5,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vpaddw ymm6,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpaddw ymm7,ymm12,YMMWORD ptr[r9+2*rbx]		;src+3*src_pitch

	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm5,ymm5,ymm8
	vpmaddwd ymm6,ymm6,ymm8
	vpmaddwd ymm7,ymm7,ymm8

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5
	vpaddd ymm2,ymm2,ymm6
	vpaddd ymm3,ymm3,ymm7
	
Resize_H_AVX2_Planar_16bits_loop_2_c:
	vextracti128 xmm4,ymm0,1
	vextracti128 xmm5,ymm1,1
	vextracti128 xmm6,ymm2,1
	vextracti128 xmm7,ymm3,1

	vphaddd xmm0,xmm0,xmm4
	vphaddd xmm1,xmm1,xmm5
	vphaddd xmm2,xmm2,xmm6
	vphaddd xmm3,xmm3,xmm7
	
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
	jnz Resize_H_AVX2_Planar_16bits_loop_1

Resize_H_AVX2_Planar_16bits_1:
	test sizeh,2
	jz Resize_H_AVX2_Planar_16bits_2

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxor ymm0,ymm0,ymm0
	vpxor ymm1,ymm1,ymm1

	shr rcx,1
	jz short Resize_H_AVX2_Planar_16bits_loop_4_b

Resize_H_AVX2_Planar_16bits_loop_4:
	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff
	vmovdqa ymm13,YMMWORD ptr[rdi+r14]	;coeff+16(x2)

	; shifttosigned + src
	vpaddw ymm4,ymm12,YMMWORD ptr[rsi]		;src
	vpaddw ymm14,ymm12,YMMWORD ptr[rsi+r14]	;src+16(x2)
	vpaddw ymm5,ymm12,YMMWORD ptr[r9]		;src+src_pitch
	vpaddw ymm15,ymm12,YMMWORD ptr[r9+r14]	;src+src_pitch+16(x2)

	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm14,ymm14,ymm13
	vpmaddwd ymm5,ymm5,ymm8
	vpmaddwd ymm15,ymm15,ymm13

	add rsi,r14
	add r9,r14
	add rdi,r14

	vpaddd ymm4,ymm4,ymm14
	vpaddd ymm5,ymm5,ymm15

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_16bits_loop_4
	
Resize_H_AVX2_Planar_16bits_loop_4_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_16bits_loop_4_c
	vmovdqa ymm8,YMMWORD ptr[rdi]		;coeff

	; shifttosigned + src
	vpaddw ymm4,ymm12,YMMWORD ptr[rsi]	;src
	vpaddw ymm5,ymm12,YMMWORD ptr[r9]	;src+src_pitch

	vpmaddwd ymm4,ymm4,ymm8
	vpmaddwd ymm5,ymm5,ymm8

	vpaddd ymm0,ymm0,ymm4
	vpaddd ymm1,ymm1,ymm5

Resize_H_AVX2_Planar_16bits_loop_4_c:	
	vextracti128 xmm4,ymm0,1
	vextracti128 xmm5,ymm1,1

	vphaddd xmm0,xmm0,xmm4
	vphaddd xmm1,xmm1,xmm5
	
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

Resize_H_AVX2_Planar_16bits_2:
	test sizeh,1
	jz Resize_H_AVX2_Planar_16bits_end

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxor ymm0,ymm0,ymm0

	shr rcx,1
	jz short Resize_H_AVX2_Planar_16bits_loop_5_b

Resize_H_AVX2_Planar_16bits_loop_5:
	; shifttosigned + src
	vpaddw ymm4,ymm12,YMMWORD ptr[rsi]		;src
	vpaddw ymm14,ymm12,YMMWORD ptr[rsi+r14]	;src+16(x2)

	vpmaddwd ymm4,ymm4,YMMWORD ptr[rdi]			;coeff
	vpmaddwd ymm14,ymm14,YMMWORD ptr[rdi+r14]	;coeff+16(x2)

	add rsi,r14

	vpaddd ymm4,ymm4,ymm14

	add rdi,r14

	vpaddd ymm0,ymm0,ymm4

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_16bits_loop_5

Resize_H_AVX2_Planar_16bits_loop_5_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_16bits_loop_5_c

	; shifttosigned + src
	vpaddw ymm4,ymm12,YMMWORD ptr[rsi]	;src

	vpmaddwd ymm4,ymm4,YMMWORD ptr[rdi]	;coeff

	vpaddd ymm0,ymm0,ymm4

Resize_H_AVX2_Planar_16bits_loop_5_c:
	vextracti128 xmm4,ymm0,1

	vphaddd xmm0,xmm0,xmm4
	
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

Resize_H_AVX2_Planar_16bits_end:
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

Resize_H_AVX2_Planar_16bits_ASM endp


;Resize_H_AVX2_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_32:dword,sizeh:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX2_Planar_32bits_ASM proc public frame

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
	sub rsp,136
	.allocstack 136
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
	jz Resize_H_AVX2_Planar_32bits_1
	
Resize_H_AVX2_Planar_32bits_loop_1:
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
	jz Resize_H_AVX2_Planar_32bits_loop_2_b

Resize_H_AVX2_Planar_32bits_loop_2:
	vmovaps ymm12,YMMWORD ptr[rdi]					;coef
	vmovaps ymm13,YMMWORD ptr[rdi+r14]				;coef+8(x4)

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi] 		;src
	vfmadd231ps ymm0,ymm13,YMMWORD ptr[rsi+r14]		;src+8(x4)
	vfmadd231ps ymm1,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vfmadd231ps ymm1,ymm13,YMMWORD ptr[r9+r14]		;src+src_pitch+8(x4)
	vfmadd231ps ymm2,ymm12,YMMWORD ptr[rax]			;src+2*src_pitch
	vfmadd231ps ymm2,ymm13,YMMWORD ptr[rax+r14]		;src+2*src_pitch+8(x4)
	vfmadd231ps ymm3,ymm12,YMMWORD ptr[rdx]			;src+3*src_pitch
	vfmadd231ps ymm3,ymm13,YMMWORD ptr[rdx+r14]		;src+3*src_pitch+8(x4)
	vfmadd231ps ymm4,ymm12,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	add rsi,r14
	vfmadd231ps ymm4,ymm13,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch+8(x4)
	vfmadd231ps ymm5,ymm12,YMMWORD ptr[r9+4*rbx]	;src+5*src_pitch
	add r9,r14
	vfmadd231ps ymm5,ymm13,YMMWORD ptr[r9+4*rbx]	;src+5*src_pitch+8(x4)
	vfmadd231ps ymm6,ymm12,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	add rax,r14
	vfmadd231ps ymm6,ymm13,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch+8(x4)
	vfmadd231ps ymm7,ymm12,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch
	add rdx,r14
	vfmadd231ps ymm7,ymm13,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch+8(x4)

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r14
	dec ecx
	jnz Resize_H_AVX2_Planar_32bits_loop_2

Resize_H_AVX2_Planar_32bits_loop_2_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_32bits_loop_2_c

	vmovaps ymm12,YMMWORD ptr[rdi]					;coef

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi] 		;src
	vfmadd231ps ymm1,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vfmadd231ps ymm2,ymm12,YMMWORD ptr[rax]			;src+2*src_pitch
	vfmadd231ps ymm3,ymm12,YMMWORD ptr[rdx]			;src+3*src_pitch
	vfmadd231ps ymm4,ymm12,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	vfmadd231ps ymm5,ymm12,YMMWORD ptr[r9+4*rbx]	;src+5*src_pitch
	vfmadd231ps ymm6,ymm12,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	vfmadd231ps ymm7,ymm12,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch

Resize_H_AVX2_Planar_32bits_loop_2_c:
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
	jnz Resize_H_AVX2_Planar_32bits_loop_1

Resize_H_AVX2_Planar_32bits_1:
	test sizeh,4
	jz Resize_H_AVX2_Planar_32bits_2

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
	jz short Resize_H_AVX2_Planar_32bits_loop_3_b

Resize_H_AVX2_Planar_32bits_loop_3:
	vmovaps ymm12,YMMWORD ptr[rdi]					;coeff
	vmovaps ymm13,YMMWORD ptr[rdi+r14]				;coeff+8(x4)

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi]			;src
	vfmadd231ps ymm0,ymm13,YMMWORD ptr[rsi+r14]		;src+8(x4)
	vfmadd231ps ymm1,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vfmadd231ps ymm1,ymm13,YMMWORD ptr[r9+r14]		;src+src_pitch+8(x4)
	vfmadd231ps ymm2,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	add rsi,r14
	vfmadd231ps ymm2,ymm13,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch+8(x4)
	vfmadd231ps ymm3,ymm12,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch
	add r9,r14
	vfmadd231ps ymm3,ymm13,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch+8(x4)

	add rdi,r14

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_32bits_loop_3
	
Resize_H_AVX2_Planar_32bits_loop_3_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_32bits_loop_3_c

	vmovaps ymm12,YMMWORD ptr[rdi]					;coeff

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi]			;src
	vfmadd231ps ymm1,ymm12,YMMWORD ptr[r9]			;src+src_pitch
	vfmadd231ps ymm2,ymm12,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vfmadd231ps ymm3,ymm12,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch

Resize_H_AVX2_Planar_32bits_loop_3_c:
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

	vhaddps xmm0,xmm0,xmm2

	sal rdx,2

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

Resize_H_AVX2_Planar_32bits_2:
	test sizeh,2
	jz Resize_H_AVX2_Planar_32bits_3

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vxorps ymm0,ymm0,ymm0
	vxorps ymm1,ymm1,ymm1

	shr rcx,1
	jz short Resize_H_AVX2_Planar_32bits_loop_4_b

Resize_H_AVX2_Planar_32bits_loop_4:
	vmovaps ymm12,YMMWORD ptr[rdi]				;coeff
	vmovaps ymm13,YMMWORD ptr[rdi+r14]			;coeff+8(x4)

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi] 	;src
	vfmadd231ps ymm0,ymm13,YMMWORD ptr[rsi+r14]	;src+8(x4)
	vfmadd231ps ymm1,ymm12,YMMWORD ptr[rsi+rbx]	;src+src_pitch
	add rsi,r14
	add rdi,r14
	vfmadd231ps ymm1,ymm13,YMMWORD ptr[rsi+rbx]	;src+src_pitch+8(x4)

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_32bits_loop_4

Resize_H_AVX2_Planar_32bits_loop_4_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_32bits_loop_4_c

	vmovaps ymm12,YMMWORD ptr[rdi]				;coeff

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi] 	;src
	vfmadd231ps ymm1,ymm12,YMMWORD ptr[rsi+rbx]	;src+src_pitch

Resize_H_AVX2_Planar_32bits_loop_4_c:
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

Resize_H_AVX2_Planar_32bits_3:
	test sizeh,1
	jz short Resize_H_AVX2_Planar_32bits_end

	mov rcx,r12 ;kernel_size_32
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vxorps ymm0,ymm0,ymm0

	shr rcx,1
	jz short Resize_H_AVX2_Planar_32bits_loop_5_b

Resize_H_AVX2_Planar_32bits_loop_5:
	vmovaps ymm12,YMMWORD ptr[rdi]			;coeff
	vmovaps ymm13,YMMWORD ptr[rdi+r14]		;coeff+8(x4)

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi]	;src
	add rsi,r14
	add rdi,r14
	vfmadd231ps ymm0,ymm13,YMMWORD ptr[rsi]	;src+8(x4)

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX2_Planar_32bits_loop_5

Resize_H_AVX2_Planar_32bits_loop_5_b:
	test r12d,1
	jz short Resize_H_AVX2_Planar_32bits_loop_5_c

	vmovaps ymm12,YMMWORD ptr[rdi]			;coeff

	vfmadd231ps ymm0,ymm12,YMMWORD ptr[rsi]	;src

Resize_H_AVX2_Planar_32bits_loop_5_c:
	vextractf128 xmm8,ymm0,1

	vhaddps xmm0,xmm0,xmm8

	vhaddps xmm0,xmm0,xmm0

	vhaddps xmm0,xmm0,xmm0

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax

Resize_H_AVX2_Planar_32bits_end:
	vmovdqa xmm13,XMMWORD ptr[rsp+112]
	vmovdqa xmm12,XMMWORD ptr[rsp+96]
	vmovdqa xmm11,XMMWORD ptr[rsp+80]
	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,136
		
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

Resize_H_AVX2_Planar_32bits_ASM endp


end
