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

	vmovdqa64 zmm0,zmm7 ;rounder
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_8bits_loop_2:
	vpmovzxbw zmm2,YMMWORD ptr[rsi]
	vpmovzxbw zmm4,YMMWORD ptr[rsi+rbx]

	vbroadcastss zmm8,dword ptr[r8+rax]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vpmaddwd zmm3,zmm3,zmm8
	vpmaddwd zmm2,zmm2,zmm8

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
	vmovdqa XMMWORD ptr[rdi],xmm0
	vmovdqa XMMWORD ptr[rdi+16],xmm3
		
	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX512_Planar_8bits_loop_1

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

	vmovdqa64 zmm0,zmm7 ;rounder
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_10to14bits_loop_2:
	vmovdqa64 zmm2,ZMMWORD ptr[rsi]
	vmovdqa64 zmm4,ZMMWORD ptr[rsi+rbx]

	vbroadcastss zmm8,dword ptr[r8+rax]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vpmaddwd zmm3,zmm3,zmm8
	vpmaddwd zmm2,zmm2,zmm8

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

	vmovdqa64 ZMMWORD ptr[rdi],zmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX512_Planar_10to14bits_loop_1

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

	vmovdqa64 zmm0,zmm7 ;rounder
	vmovdqa64 zmm1,zmm7

Resize_V_AVX512_Planar_16bits_loop_2:
	vmovdqa64 zmm2,ZMMWORD ptr[rsi]
	vmovdqa64 zmm4,ZMMWORD ptr[rsi+rbx]

	vbroadcastss zmm10,dword ptr[r8+rax]

	vpunpckhwd zmm3,zmm2,zmm4
	vpunpcklwd zmm2,zmm2,zmm4

	vpaddw zmm3,zmm3,zmm8 ;shifttosigned
	vpaddw zmm2,zmm2,zmm8

	vpmaddwd zmm3,zmm3,zmm10
	vpmaddwd zmm2,zmm2,zmm10

	vpaddd zmm1,zmm1,zmm3
	vpaddd zmm0,zmm0,zmm2

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX512_Planar_16bits_loop_2

	vpaddw zmm0,zmm0,zmm9 ;shiftfromsigned
	vpaddw zmm1,zmm1,zmm9

	vpsrad zmm0,zmm0,13 ;FPScale16bits = 13
	vpsrad zmm1,zmm1,13

	vpackusdw zmm0,zmm0,zmm1

	vpmaxuw zmm0,zmm0,zmm5
	vpminuw zmm0,zmm0,zmm6

	vmovdqa64 ZMMWORD ptr[rdi],zmm0

	add rdi,r12 ; dst + x
	add r10,r12 ; src + x
	dec r9d
	jnz Resize_V_AVX512_Planar_16bits_loop_1

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
	vbroadcastss zmm4,dword ptr[r8+rax]
	vbroadcastss zmm5,dword ptr[r8+rax+4]

	vfmadd231ps zmm0,zmm4,ZMMWORD ptr[rsi]
	vfmadd231ps zmm1,zmm5,ZMMWORD ptr[rsi+rbx]

	add rsi,rbx
	add rax,rdx
	add rsi,rbx
	loop Resize_V_AVX512_Planar_32bits_loop_2

	vaddps zmm0,zmm0,zmm1

	vmovaps ZMMWORD ptr[rdi],zmm0

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Resize_H_AVX512_Planar_8bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_64:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX512_Planar_8bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_64 equ dword ptr[rbp+56]
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

	mov rsi,valmin
	vbroadcastss ymm17,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss ymm18,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss ymm19,dword ptr[rsi]
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_64 ;kernel_size_64 = (kernel_size + 31) >> 5
	mov r13d,sizeh
	mov r14,32
	mov r15,64
	shr r13d,3 ;r13d = sizeh/8
	jz Resize_H_AVX512_Planar_8bits_1
	
Resize_H_AVX512_Planar_8bits_loop_1:
	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch
	mov rax,r9
	add rax,rbx ;rax=src+2*pitch
	mov rdx,rax
	add rdx,rbx ;rdx=src+3*pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1
	vpxord zmm2,zmm2,zmm2
	vpxord zmm3,zmm3,zmm3
	vpxord zmm4,zmm4,zmm4
	vpxord zmm5,zmm5,zmm5
	vpxord zmm6,zmm6,zmm6
	vpxord zmm7,zmm7,zmm7

Resize_H_AVX512_Planar_8bits_loop_2:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]		;coeff

	vpmovzxbw zmm8,YMMWORD ptr[rsi] 		;src
	vpmovzxbw zmm9,YMMWORD ptr[r9]			;src+src_pitch
	vpmovzxbw zmm10,YMMWORD ptr[rax]		;src+2*src_pitch
	vpmovzxbw zmm11,YMMWORD ptr[rdx]		;src+3*src_pitch
	vpmovzxbw zmm12,YMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	vpmovzxbw zmm13,YMMWORD ptr[r9+4*rbx]	;src+5*src_pitch
	vpmovzxbw zmm14,YMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	vpmovzxbw zmm15,YMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch

	vpmaddwd zmm8,zmm8,zmm16
	vpmaddwd zmm9,zmm9,zmm16
	vpmaddwd zmm10,zmm10,zmm16
	vpmaddwd zmm11,zmm11,zmm16
	vpmaddwd zmm12,zmm12,zmm16
	vpmaddwd zmm13,zmm13,zmm16
	vpmaddwd zmm14,zmm14,zmm16
	vpmaddwd zmm15,zmm15,zmm16

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9
	vpaddd zmm2,zmm2,zmm10
	vpaddd zmm3,zmm3,zmm11
	vpaddd zmm4,zmm4,zmm12
	vpaddd zmm5,zmm5,zmm13
	vpaddd zmm6,zmm6,zmm14
	vpaddd zmm7,zmm7,zmm15

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r15
	dec ecx
	jnz Resize_H_AVX512_Planar_8bits_loop_2

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	vextracti32x8 ymm10,zmm2,1
	vextracti32x8 ymm11,zmm3,1
	vextracti32x8 ymm12,zmm4,1
	vextracti32x8 ymm13,zmm5,1
	vextracti32x8 ymm14,zmm6,1
	vextracti32x8 ymm15,zmm7,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	vphaddd ymm2,ymm2,ymm10
	vphaddd ymm3,ymm3,ymm11
	vphaddd ymm4,ymm4,ymm12
	vphaddd ymm5,ymm5,ymm13
	vphaddd ymm6,ymm6,ymm14
	vphaddd ymm7,ymm7,ymm15
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1
	vextracti128 xmm10,ymm2,1
	vextracti128 xmm11,ymm3,1
	vextracti128 xmm12,ymm4,1
	vextracti128 xmm13,ymm5,1
	vextracti128 xmm14,ymm6,1
	vextracti128 xmm15,ymm7,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	vphaddd xmm2,xmm2,xmm10
	vphaddd xmm3,xmm3,xmm11
	vphaddd xmm4,xmm4,xmm12
	vphaddd xmm5,xmm5,xmm13
	vphaddd xmm6,xmm6,xmm14
	vphaddd xmm7,xmm7,xmm15
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	vphaddd xmm4,xmm4,xmm4
	vphaddd xmm5,xmm5,xmm5
	vphaddd xmm6,xmm6,xmm6
	vphaddd xmm7,xmm7,xmm7
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68
	vshufps xmm4,xmm4,xmm5,68
	vshufps xmm6,xmm6,xmm7,68

	vinserti128 ymm0,ymm0,xmm2,1
	vinserti128 ymm4,ymm4,xmm6,1

	vphaddd ymm0,ymm0,ymm4
	
	vpaddd ymm0,ymm0,ymm19 ;rounder
	vpsrad ymm0,ymm0,14 ;FPScale8bits = 14
	
	mov rdx,rbx
	
	vpackusdw ymm0,ymm0,ymm0
	
	sal rdx,3
	
	vpackuswb ymm0,ymm0,ymm0
	
	add r10,rdx
	
	vpmaxub ymm0,ymm0,ymm17

	mov rdx,dst_pitch

	vpminub ymm0,ymm0,ymm18
	
	vextracti128 xmm1,ymm0,1
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,1
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm1,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm1,1
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,2
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,3
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm1,2
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm1,3
	mov byte ptr[r11],al
	add r11,rdx
	
	dec r13d
	jnz Resize_H_AVX512_Planar_8bits_loop_1

Resize_H_AVX512_Planar_8bits_1:
	test sizeh,4
	jz Resize_H_AVX512_Planar_8bits_2

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1
	vpxord zmm2,zmm2,zmm2
	vpxord zmm3,zmm3,zmm3

Resize_H_AVX512_Planar_8bits_loop_3:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]		;coeff

	vpmovzxbw zmm8,YMMWORD ptr[rsi] 		;src
	vpmovzxbw zmm9,YMMWORD ptr[r9]			;src+src_pitch
	vpmovzxbw zmm10,YMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpmovzxbw zmm11,YMMWORD ptr[r9+2*rbx]	;src+3*src_pitch
	
	vpmaddwd zmm8,zmm8,zmm16
	vpmaddwd zmm9,zmm9,zmm16
	vpmaddwd zmm10,zmm10,zmm16
	vpmaddwd zmm11,zmm11,zmm16

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9
	vpaddd zmm2,zmm2,zmm10
	vpaddd zmm3,zmm3,zmm11

	add rsi,r14
	add r9,r14
	add rdi,r15
	loop Resize_H_AVX512_Planar_8bits_loop_3

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	vextracti32x8 ymm10,zmm2,1
	vextracti32x8 ymm11,zmm3,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	vphaddd ymm2,ymm2,ymm10
	vphaddd ymm3,ymm3,ymm11
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1
	vextracti128 xmm10,ymm2,1
	vextracti128 xmm11,ymm3,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	vphaddd xmm2,xmm2,xmm10
	vphaddd xmm3,xmm3,xmm11
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	vinserti128 ymm0,ymm0,xmm2,1

	vphaddd ymm0,ymm0,ymm0
	
	vpaddd ymm0,ymm0,ymm19 ;rounder
	vpsrad ymm0,ymm0,14 ;FPScale8bits = 14
	
	mov rdx,rbx
	
	vpackusdw ymm0,ymm0,ymm0
	
	sal rdx,2
	
	vpackuswb ymm0,ymm0,ymm0
	
	add r10,rdx
	
	vpmaxub ymm0,ymm0,ymm17

	mov rdx,dst_pitch

	vpminub ymm0,ymm0,ymm18
	
	vextracti128 xmm1,ymm0,1
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,1
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm1,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm1,1
	mov byte ptr[r11],al
	add r11,rdx

Resize_H_AVX512_Planar_8bits_2:
	test sizeh,2
	jz Resize_H_AVX512_Planar_8bits_3

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_8bits_loop_4:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]	;coeff

	vpmovzxbw zmm8,YMMWORD ptr[rsi] 	;src
	vpmovzxbw zmm9,YMMWORD ptr[rsi+rbx]	;src+src_pitch

	vpmaddwd zmm8,zmm8,zmm16
	vpmaddwd zmm9,zmm9,zmm16

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9

	add rsi,r14
	add rdi,r15
	loop Resize_H_AVX512_Planar_8bits_loop_4

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	vphaddd xmm0,xmm0,xmm1
	
	vpaddd xmm0,xmm0,xmm19 ;rounder
	vpsrad xmm0,xmm0,14 ;FPScale18bits = 14
	
	mov rdx,rbx
	
	vpackusdw xmm0,xmm0,xmm0
	
	sal rdx,1
	
	vpackuswb xmm0,xmm0,xmm0
	
	vpmaxub xmm0,xmm0,xmm17
	
	add r10,rdx
	
	vpminub xmm0,xmm0,xmm18

	mov rdx,dst_pitch
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al
	add r11,rdx
	vpextrb eax,xmm0,2
	mov byte ptr[r11],al
	add r11,rdx

Resize_H_AVX512_Planar_8bits_3:
	test sizeh,1
	jz Resize_H_AVX512_Planar_8bits_end

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxord zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_8bits_loop_5:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]	;coeff

	vpmovzxbw zmm8,YMMWORD ptr[rsi]		;src

	vpmaddwd zmm8,zmm8,zmm16

	vpaddd zmm0,zmm0,zmm8

	add rsi,r14
	add rdi,r15
	loop Resize_H_AVX512_Planar_8bits_loop_5

	vextracti32x8 ymm8,zmm0,1
	
	vphaddd ymm0,ymm0,ymm8
	
	vextracti128 xmm8,ymm0,1

	vphaddd xmm0,xmm0,xmm8
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0
	
	vpaddd xmm0,xmm0,xmm19 ;rounder
	vpsrad xmm0,xmm0,14 ;FPScale8bits = 14
	
	vpackusdw xmm0,xmm0,xmm0
	vpackuswb xmm0,xmm0,xmm0
	
	vpmaxub xmm0,xmm0,xmm17
	vpminub xmm0,xmm0,xmm18
	
	vpextrb eax,xmm0,0
	mov byte ptr[r11],al

Resize_H_AVX512_Planar_8bits_end:
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

Resize_H_AVX512_Planar_8bits_ASM endp


;Resize_H_AVX512_Planar_10to14bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_64:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX512_Planar_10to14bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_64 equ dword ptr[rbp+56]
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

	mov rsi,valmin
	vbroadcastss ymm17,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss ymm18,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss ymm19,dword ptr[rsi]
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_64 ;kernel_size_64 = (kernel_size + 31) >> 5
	mov r13d,sizeh
	mov r14,64
	mov r15,dst_pitch
	shr r13d,3 ;r13d = sizeh/8
	jz Resize_H_AVX512_Planar_10to14bits_1
	
Resize_H_AVX512_Planar_10to14bits_loop_1:
	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch
	mov rax,r9
	add rax,rbx ;rax=src+2*pitch
	mov rdx,rax
	add rdx,rbx ;rdx=src+3*pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1
	vpxord zmm2,zmm2,zmm2
	vpxord zmm3,zmm3,zmm3
	vpxord zmm4,zmm4,zmm4
	vpxord zmm5,zmm5,zmm5
	vpxord zmm6,zmm6,zmm6
	vpxord zmm7,zmm7,zmm7

Resize_H_AVX512_Planar_10to14bits_loop_2:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]			;coef

	vpmaddwd zmm8,zmm16,ZMMWORD ptr[rsi] 		;src
	vpmaddwd zmm9,zmm16,ZMMWORD ptr[r9]			;src+src_pitch
	vpmaddwd zmm10,zmm16,ZMMWORD ptr[rax]		;src+2*src_pitch
	vpmaddwd zmm11,zmm16,ZMMWORD ptr[rdx]		;src+3*src_pitch
	vpmaddwd zmm12,zmm16,ZMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	vpmaddwd zmm13,zmm16,ZMMWORD ptr[r9+4*rbx]	;src+5*src_pitch
	vpmaddwd zmm14,zmm16,ZMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	vpmaddwd zmm15,zmm16,ZMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9
	vpaddd zmm2,zmm2,zmm10
	vpaddd zmm3,zmm3,zmm11
	vpaddd zmm4,zmm4,zmm12
	vpaddd zmm5,zmm5,zmm13
	vpaddd zmm6,zmm6,zmm14
	vpaddd zmm7,zmm7,zmm15

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_10to14bits_loop_2

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	vextracti32x8 ymm10,zmm2,1
	vextracti32x8 ymm11,zmm3,1
	vextracti32x8 ymm12,zmm4,1
	vextracti32x8 ymm13,zmm5,1
	vextracti32x8 ymm14,zmm6,1
	vextracti32x8 ymm15,zmm7,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	vphaddd ymm2,ymm2,ymm10
	vphaddd ymm3,ymm3,ymm11
	vphaddd ymm4,ymm4,ymm12
	vphaddd ymm5,ymm5,ymm13
	vphaddd ymm6,ymm6,ymm14
	vphaddd ymm7,ymm7,ymm15
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1
	vextracti128 xmm10,ymm2,1
	vextracti128 xmm11,ymm3,1
	vextracti128 xmm12,ymm4,1
	vextracti128 xmm13,ymm5,1
	vextracti128 xmm14,ymm6,1
	vextracti128 xmm15,ymm7,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	vphaddd xmm2,xmm2,xmm10
	vphaddd xmm3,xmm3,xmm11
	vphaddd xmm4,xmm4,xmm12
	vphaddd xmm5,xmm5,xmm13
	vphaddd xmm6,xmm6,xmm14
	vphaddd xmm7,xmm7,xmm15
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	vphaddd xmm4,xmm4,xmm4
	vphaddd xmm5,xmm5,xmm5
	vphaddd xmm6,xmm6,xmm6
	vphaddd xmm7,xmm7,xmm7
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68
	vshufps xmm4,xmm4,xmm5,68
	vshufps xmm6,xmm6,xmm7,68

	vinserti128 ymm0,ymm0,xmm2,1
	vinserti128 ymm4,ymm4,xmm6,1

	vphaddd ymm0,ymm0,ymm4
	
	vpaddd ymm0,ymm0,ymm19 ;rounder
	vpsrad ymm0,ymm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw ymm0,ymm0,ymm0

	sal rdx,3

	vpmaxuw ymm0,ymm0,ymm17

	add r10,rdx

	vpminuw ymm0,ymm0,ymm18
	
	vextracti128 xmm1,ymm0,1
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,3
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,2
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,3
	mov word ptr[r11],ax
	add r11,r15
	
	dec r13d
	jnz Resize_H_AVX512_Planar_10to14bits_loop_1

Resize_H_AVX512_Planar_10to14bits_1:
	test sizeh,4
	jz Resize_H_AVX512_Planar_10to14bits_2

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1
	vpxord zmm2,zmm2,zmm2
	vpxord zmm3,zmm3,zmm3

Resize_H_AVX512_Planar_10to14bits_loop_3:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]			;coeff

	vpmaddwd zmm8,zmm16,ZMMWORD ptr[rsi]		;src
	vpmaddwd zmm9,zmm16,ZMMWORD ptr[r9]			;src+src_pitch
	vpmaddwd zmm10,zmm16,ZMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vpmaddwd zmm11,zmm16,ZMMWORD ptr[r9+2*rbx]	;src+3*src_pitch

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9
	vpaddd zmm2,zmm2,zmm10
	vpaddd zmm3,zmm3,zmm11

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_10to14bits_loop_3

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	vextracti32x8 ymm10,zmm2,1
	vextracti32x8 ymm11,zmm3,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	vphaddd ymm2,ymm2,ymm10
	vphaddd ymm3,ymm3,ymm11
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1
	vextracti128 xmm10,ymm2,1
	vextracti128 xmm11,ymm3,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	vphaddd xmm2,xmm2,xmm10
	vphaddd xmm3,xmm3,xmm11
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	vinserti128 ymm0,ymm0,xmm2,1

	vphaddd ymm0,ymm0,ymm0
	
	vpaddd ymm0,ymm0,ymm19 ;rounder
	vpsrad ymm0,ymm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw ymm0,ymm0,ymm0

	sal rdx,2

	vpmaxuw ymm0,ymm0,ymm17

	add r10,rdx

	vpminuw ymm0,ymm0,ymm18
	
	vextracti128 xmm1,ymm0,1
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,1
	mov word ptr[r11],ax
	add r11,r15

Resize_H_AVX512_Planar_10to14bits_2:
	test sizeh,2
	jz Resize_H_AVX512_Planar_10to14bits_3

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_10to14bits_loop_4:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]			;coeff

	vpmaddwd zmm8,zmm16,ZMMWORD ptr[rsi] 		;src
	vpmaddwd zmm9,zmm16,ZMMWORD ptr[rsi+rbx]	;src+src_pitch

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_10to14bits_loop_4

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	vphaddd xmm0,xmm0,xmm1
	
	vpaddd xmm0,xmm0,xmm19 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw xmm0,xmm0,xmm0

	sal rdx,1

	vpmaxuw xmm0,xmm0,xmm17

	add r10,rdx

	vpminuw xmm0,xmm0,xmm18
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15

Resize_H_AVX512_Planar_10to14bits_3:
	test sizeh,1
	jz Resize_H_AVX512_Planar_10to14bits_end

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxord zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_10to14bits_loop_5:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]		;coeff

	vpmaddwd zmm8,zmm16,ZMMWORD ptr[rsi] 	;src

	vpaddd zmm0,zmm0,zmm8

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_10to14bits_loop_5

	vextracti32x8 ymm8,zmm0,1
	
	vphaddd ymm0,ymm0,ymm8
	
	vextracti128 xmm8,ymm0,1

	vphaddd xmm0,xmm0,xmm8
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0
	
	vpaddd xmm0,xmm0,xmm19 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	
	vpackusdw xmm0,xmm0,xmm0

	vpmaxuw xmm0,xmm0,xmm17
	vpminuw xmm0,xmm0,xmm18
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax

Resize_H_AVX512_Planar_10to14bits_end:
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

Resize_H_AVX512_Planar_10to14bits_ASM endp


;Resize_H_AVX512_Planar_16bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_64:dword,sizeh:dword,valmin:dword,valmax:dword,rounder:dword
;	shifttosigned:dword,shiftfromsigned:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX512_Planar_16bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_64 equ dword ptr[rbp+56]
sizeh equ dword ptr[rbp+64]
valmin equ qword ptr[rbp+72]
valmax equ qword ptr[rbp+80]
rounder equ qword ptr[rbp+88]
shifttosigned equ qword ptr[rbp+96]
shiftfromsigned equ qword ptr[rbp+104]

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

	mov rsi,valmin
	vbroadcastss ymm17,dword ptr[rsi]
	mov rsi,valmax
	vbroadcastss ymm18,dword ptr[rsi]
	mov rsi,rounder
	vbroadcastss ymm19,dword ptr[rsi]
	mov rsi,shifttosigned
	vbroadcastss zmm20,dword ptr[rsi]
	mov rsi,shiftfromsigned
	vbroadcastss ymm21,dword ptr[rsi]
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_64 ;kernel_size_64 = (kernel_size + 31) >> 5
	mov r13d,sizeh
	mov r14,64
	mov r15,dst_pitch
	shr r13d,3 ;r13d = sizeh/8
	jz Resize_H_AVX512_Planar_16bits_1
	
Resize_H_AVX512_Planar_16bits_loop_1:
	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch
	mov rax,r9
	add rax,rbx ;rax=src+2*pitch
	mov rdx,rax
	add rdx,rbx ;rdx=src+3*pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1
	vpxord zmm2,zmm2,zmm2
	vpxord zmm3,zmm3,zmm3
	vpxord zmm4,zmm4,zmm4
	vpxord zmm5,zmm5,zmm5
	vpxord zmm6,zmm6,zmm6
	vpxord zmm7,zmm7,zmm7

Resize_H_AVX512_Planar_16bits_loop_2:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]			;coeff

	; shifttosigned + src
	vpaddw zmm8,zmm20,ZMMWORD ptr[rsi] 			;src
	vpaddw zmm9,zmm20,ZMMWORD ptr[r9]			;src+src_pitch
	vpaddw zmm10,zmm20,ZMMWORD ptr[rax]			;src+2*src_pitch
	vpaddw zmm11,zmm20,ZMMWORD ptr[rdx]			;src+3*src_pitch
	vpaddw zmm12,zmm20,ZMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	vpaddw zmm13,zmm20,ZMMWORD ptr[r9+4*rbx]	;src+5*src_pitch
	vpaddw zmm14,zmm20,ZMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	vpaddw zmm15,zmm20,ZMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch

	vpmaddwd zmm8,zmm8,zmm16
	vpmaddwd zmm9,zmm9,zmm16
	vpmaddwd zmm10,zmm10,zmm16
	vpmaddwd zmm11,zmm11,zmm16
	vpmaddwd zmm12,zmm12,zmm16
	vpmaddwd zmm13,zmm13,zmm16
	vpmaddwd zmm14,zmm14,zmm16
	vpmaddwd zmm15,zmm15,zmm16

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9
	vpaddd zmm2,zmm2,zmm10
	vpaddd zmm3,zmm3,zmm11
	vpaddd zmm4,zmm4,zmm12
	vpaddd zmm5,zmm5,zmm13
	vpaddd zmm6,zmm6,zmm14
	vpaddd zmm7,zmm7,zmm15

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r14
	dec ecx
	jnz Resize_H_AVX512_Planar_16bits_loop_2

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	vextracti32x8 ymm10,zmm2,1
	vextracti32x8 ymm11,zmm3,1
	vextracti32x8 ymm12,zmm4,1
	vextracti32x8 ymm13,zmm5,1
	vextracti32x8 ymm14,zmm6,1
	vextracti32x8 ymm15,zmm7,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	vphaddd ymm2,ymm2,ymm10
	vphaddd ymm3,ymm3,ymm11
	vphaddd ymm4,ymm4,ymm12
	vphaddd ymm5,ymm5,ymm13
	vphaddd ymm6,ymm6,ymm14
	vphaddd ymm7,ymm7,ymm15
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1
	vextracti128 xmm10,ymm2,1
	vextracti128 xmm11,ymm3,1
	vextracti128 xmm12,ymm4,1
	vextracti128 xmm13,ymm5,1
	vextracti128 xmm14,ymm6,1
	vextracti128 xmm15,ymm7,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	vphaddd xmm2,xmm2,xmm10
	vphaddd xmm3,xmm3,xmm11
	vphaddd xmm4,xmm4,xmm12
	vphaddd xmm5,xmm5,xmm13
	vphaddd xmm6,xmm6,xmm14
	vphaddd xmm7,xmm7,xmm15
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	vphaddd xmm4,xmm4,xmm4
	vphaddd xmm5,xmm5,xmm5
	vphaddd xmm6,xmm6,xmm6
	vphaddd xmm7,xmm7,xmm7
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68
	vshufps xmm4,xmm4,xmm5,68
	vshufps xmm6,xmm6,xmm7,68

	vinserti128 ymm0,ymm0,xmm2,1
	vinserti128 ymm4,ymm4,xmm6,1

	vphaddd ymm0,ymm0,ymm4

	vpaddd ymm0,ymm0,ymm21 ;ShiftfromSigned

	vpaddd ymm0,ymm0,ymm19 ;rounder
	vpsrad ymm0,ymm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw ymm0,ymm0,ymm0

	sal rdx,3

	vpmaxuw ymm0,ymm0,ymm17

	add r10,rdx

	vpminuw ymm0,ymm0,ymm18

	vextracti128 xmm1,ymm0,1

	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,3
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,2
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,3
	mov word ptr[r11],ax
	add r11,r15
	
	dec r13d
	jnz Resize_H_AVX512_Planar_16bits_loop_1

Resize_H_AVX512_Planar_16bits_1:
	test sizeh,4
	jz Resize_H_AVX512_Planar_16bits_2

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch
	mov rax,r9
	add rax,rbx ;rax=src+2*pitch
	mov rdx,rax
	add rdx,rbx ;rdx=src+3*pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1
	vpxord zmm2,zmm2,zmm2
	vpxord zmm3,zmm3,zmm3

Resize_H_AVX512_Planar_16bits_loop_3:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]	;coeff

	; shifttosigned + src
	vpaddw zmm8,zmm20,ZMMWORD ptr[rsi] 	;src
	vpaddw zmm9,zmm20,ZMMWORD ptr[r9]	;src+src_pitch
	vpaddw zmm10,zmm20,ZMMWORD ptr[rax]	;src+2*src_pitch
	vpaddw zmm11,zmm20,ZMMWORD ptr[rdx]	;src+3*src_pitch

	vpmaddwd zmm8,zmm8,zmm16
	vpmaddwd zmm9,zmm9,zmm16
	vpmaddwd zmm10,zmm10,zmm16
	vpmaddwd zmm11,zmm11,zmm16

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9
	vpaddd zmm2,zmm2,zmm10
	vpaddd zmm3,zmm3,zmm11

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_16bits_loop_3

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	vextracti32x8 ymm10,zmm2,1
	vextracti32x8 ymm11,zmm3,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	vphaddd ymm2,ymm2,ymm10
	vphaddd ymm3,ymm3,ymm11
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1
	vextracti128 xmm10,ymm2,1
	vextracti128 xmm11,ymm3,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	vphaddd xmm2,xmm2,xmm10
	vphaddd xmm3,xmm3,xmm11
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1
	vphaddd xmm2,xmm2,xmm2
	vphaddd xmm3,xmm3,xmm3
	
	vshufps xmm0,xmm0,xmm1,68
	vshufps xmm2,xmm2,xmm3,68

	vinserti128 ymm0,ymm0,xmm2,1

	vphaddd ymm0,ymm0,ymm0

	vpaddd ymm0,ymm0,ymm21 ;ShiftfromSigned

	vpaddd ymm0,ymm0,ymm19 ;rounder
	vpsrad ymm0,ymm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw ymm0,ymm0,ymm0

	sal rdx,2

	vpmaxuw ymm0,ymm0,ymm17

	add r10,rdx

	vpminuw ymm0,ymm0,ymm18

	vextracti128 xmm1,ymm0,1

	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,1
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm1,1
	mov word ptr[r11],ax
	add r11,r15

Resize_H_AVX512_Planar_16bits_2:
	test sizeh,2
	jz Resize_H_AVX512_Planar_16bits_3

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vpxord zmm0,zmm0,zmm0
	vpxord zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_16bits_loop_4:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]	;coeff

	; shifttosigned + src
	vpaddw zmm8,zmm20,ZMMWORD ptr[rsi]	;src
	vpaddw zmm9,zmm20,ZMMWORD ptr[r9]	;src+src_pitch

	vpmaddwd zmm8,zmm8,zmm16
	vpmaddwd zmm9,zmm9,zmm16

	vpaddd zmm0,zmm0,zmm8
	vpaddd zmm1,zmm1,zmm9

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_16bits_loop_4

	vextracti32x8 ymm8,zmm0,1
	vextracti32x8 ymm9,zmm1,1
	
	vphaddd ymm0,ymm0,ymm8
	vphaddd ymm1,ymm1,ymm9
	
	vextracti128 xmm8,ymm0,1
	vextracti128 xmm9,ymm1,1

	vphaddd xmm0,xmm0,xmm8
	vphaddd xmm1,xmm1,xmm9
	
	vphaddd xmm0,xmm0,xmm0
	vphaddd xmm1,xmm1,xmm1

	vphaddd xmm0,xmm0,xmm1

	vpaddd xmm0,xmm0,xmm21 ;ShiftfromSigned

	vpaddd xmm0,xmm0,xmm19 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13

	mov rdx,rbx

	vpackusdw xmm0,xmm0,xmm0

	sal rdx,1
	
	vpmaxuw xmm0,xmm0,xmm17

	add r10,rdx

	vpminuw xmm0,xmm0,xmm18
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax
	add r11,r15
	vpextrw eax,xmm0,2
	mov word ptr[r11],ax
	add r11,r15

Resize_H_AVX512_Planar_16bits_3:
	test sizeh,1
	jz Resize_H_AVX512_Planar_16bits_end

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vpxord zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_16bits_loop_5:
	vmovdqa64 zmm16,ZMMWORD ptr[rdi]	;coeff

	; shifttosigned + src
	vpaddw zmm8,zmm20,ZMMWORD ptr[rsi] 	;src

	vpmaddwd zmm8,zmm8,zmm16

	vpaddd zmm0,zmm0,zmm8

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_16bits_loop_5

	vextracti32x8 ymm8,zmm0,1
	
	vphaddd ymm0,ymm0,ymm8
	
	vextracti128 xmm8,ymm0,1

	vphaddd xmm0,xmm0,xmm8
	
	vphaddd xmm0,xmm0,xmm0

	vphaddd xmm0,xmm0,xmm0

	vpaddd xmm0,xmm0,xmm21 ;ShiftfromSigned

	vpaddd xmm0,xmm0,xmm19 ;rounder
	vpsrad xmm0,xmm0,13 ;FPScale16bits = 13
	
	mov rdx,dst_pitch
	
	vpackusdw xmm0,xmm0,xmm0
	
	vpmaxuw xmm0,xmm0,xmm17
	vpminuw xmm0,xmm0,xmm18
	
	vpextrw eax,xmm0,0
	mov word ptr[r11],ax

Resize_H_AVX512_Planar_16bits_end:
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

Resize_H_AVX512_Planar_16bits_ASM endp


;Resize_H_AVX512_Planar_32bits_ASM proc src:dword,dst:dword,coeff:dword,src_pitch:dword,dst_pitch:dword,
;	kernel_size_64:dword,sizeh:dword

; src = rcx
; dst = rdx
; coeff = r8
; src_pitch = r9

Resize_H_AVX512_Planar_32bits_ASM proc public frame

dst_pitch equ qword ptr[rbp+48]
kernel_size_64 equ dword ptr[rbp+56]
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
	sub rsp,120
	.allocstack 120
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
	.endprolog
	
	mov r10,rcx ;r10=src
	mov r11,rdx ;r11=dst
	mov rbx,r9 ;rbx=src_pitch
	xor r12,r12
	xor r13,r13
	mov r12d,kernel_size_64 ;kernel_size_64 = (kernel_size + 15) >> 4
	mov r13d,sizeh
	mov r14,64
	mov r15,dst_pitch
	shr r13d,3 ;r13d = sizeh/8
	jz Resize_H_AVX512_Planar_32bits_1
	
Resize_H_AVX512_Planar_32bits_loop_1:
	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch
	mov rax,r9
	add rax,rbx ;rax=src+2*pitch
	mov rdx,rax
	add rdx,rbx ;rdx=src+3*pitch

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1
	vxorps zmm2,zmm2,zmm2
	vxorps zmm3,zmm3,zmm3
	vxorps zmm4,zmm4,zmm4
	vxorps zmm5,zmm5,zmm5
	vxorps zmm6,zmm6,zmm6
	vxorps zmm7,zmm7,zmm7

Resize_H_AVX512_Planar_32bits_loop_2:
	vmovaps zmm12,ZMMWORD ptr[rdi]					;coef

	vfmadd231ps zmm0,zmm12,ZMMWORD ptr[rsi] 		;src
	vfmadd231ps zmm1,zmm12,ZMMWORD ptr[r9]			;src+src_pitch
	vfmadd231ps zmm2,zmm12,ZMMWORD ptr[rax]			;src+2*src_pitch
	vfmadd231ps zmm3,zmm12,ZMMWORD ptr[rdx]			;src+3*src_pitch
	vfmadd231ps zmm4,zmm12,ZMMWORD ptr[rsi+4*rbx]	;src+4*src_pitch
	vfmadd231ps zmm5,zmm12,ZMMWORD ptr[r9+4*rbx]	;src+5*src_pitch
	vfmadd231ps zmm6,zmm12,ZMMWORD ptr[rax+4*rbx]	;src+6*src_pitch
	vfmadd231ps zmm7,zmm12,ZMMWORD ptr[rdx+4*rbx]	;src+7*src_pitch

	add rsi,r14
	add r9,r14
	add rax,r14
	add rdx,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_32bits_loop_2

	vextractf32x8 ymm8,zmm0,1
	vextractf32x8 ymm9,zmm1,1
	vextractf32x8 ymm10,zmm2,1
	vextractf32x8 ymm11,zmm3,1
	vhaddps ymm0,ymm0,ymm8
	vhaddps ymm1,ymm1,ymm9
	vhaddps ymm2,ymm2,ymm10
	vhaddps ymm3,ymm3,ymm11
	vextractf32x8 ymm8,zmm4,1
	vextractf32x8 ymm9,zmm5,1
	vextractf32x8 ymm10,zmm6,1
	vextractf32x8 ymm11,zmm7,1	
	vhaddps ymm4,ymm4,ymm8
	vhaddps ymm5,ymm5,ymm9
	vhaddps ymm6,ymm6,ymm10
	vhaddps ymm7,ymm7,ymm11
	
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

	vinsertf128 ymm0,ymm0,xmm2,1
	vinsertf128 ymm4,ymm4,xmm6,1

	mov rdx,rbx

	vhaddps ymm0,ymm0,ymm4

	sal rdx,3

	vextractf128 xmm1,ymm0,1

	add r10,rdx

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,1
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm1,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm1,1
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,2
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,3
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm1,2
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm1,3
	mov dword ptr[r11],eax
	add r11,r15
	
	dec r13d
	jnz Resize_H_AVX512_Planar_32bits_loop_1

Resize_H_AVX512_Planar_32bits_1:
	test sizeh,4
	jz Resize_H_AVX512_Planar_32bits_2

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src
	mov r9,r10
	add r9,rbx 	;r9=src+src_pitch

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1
	vxorps zmm2,zmm2,zmm2
	vxorps zmm3,zmm3,zmm3

Resize_H_AVX512_Planar_32bits_loop_3:
	vmovaps zmm12,ZMMWORD ptr[rdi]					;coeff

	vfmadd231ps zmm0,zmm12,ZMMWORD ptr[rsi]			;src
	vfmadd231ps zmm1,zmm12,ZMMWORD ptr[r9]			;src+src_pitch
	vfmadd231ps zmm2,zmm12,ZMMWORD ptr[rsi+2*rbx]	;src+2*src_pitch
	vfmadd231ps zmm3,zmm12,ZMMWORD ptr[r9+2*rbx]	;src+3*src_pitch

	add rsi,r14
	add r9,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_32bits_loop_3

	vextractf32x8 ymm8,zmm0,1
	vextractf32x8 ymm9,zmm1,1
	vextractf32x8 ymm10,zmm2,1
	vextractf32x8 ymm11,zmm3,1
	
	vhaddps ymm0,ymm0,ymm8
	vhaddps ymm1,ymm1,ymm9
	vhaddps ymm2,ymm2,ymm10
	vhaddps ymm3,ymm3,ymm11
	
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

	vinsertf128 ymm0,ymm0,xmm2,1

	sal rdx,2

	vhaddps ymm0,ymm0,ymm0

	add r10,rdx

	vextracti128 xmm1,ymm0,1
	
	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm0,1
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm1,0
	mov dword ptr[r11],eax
	add r11,r15
	vpextrd eax,xmm1,1
	mov dword ptr[r11],eax
	add r11,r15

Resize_H_AVX512_Planar_32bits_2:
	test sizeh,2
	jz Resize_H_AVX512_Planar_32bits_3

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vxorps zmm0,zmm0,zmm0
	vxorps zmm1,zmm1,zmm1

Resize_H_AVX512_Planar_32bits_loop_4:
	vmovaps zmm12,ZMMWORD ptr[rdi]				;coeff

	vfmadd231ps zmm0,zmm12,ZMMWORD ptr[rsi] 	;src
	vfmadd231ps zmm1,zmm12,ZMMWORD ptr[rsi+rbx]	;src+src_pitch

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_32bits_loop_4

	vextractf32x8 ymm8,zmm0,1
	vextractf32x8 ymm9,zmm1,1
	
	vhaddps ymm0,ymm0,ymm8
	vhaddps ymm1,ymm1,ymm9
	
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

Resize_H_AVX512_Planar_32bits_3:
	test sizeh,1
	jz short Resize_H_AVX512_Planar_32bits_end

	mov rcx,r12 ;kernel_size_64
	mov rdi,r8 	;rdi=coeff
	mov rsi,r10	;rsi=src

	vxorps zmm0,zmm0,zmm0

Resize_H_AVX512_Planar_32bits_loop_5:
	vmovaps zmm12,ZMMWORD ptr[rdi]			;coeff

	vfmadd231ps zmm0,zmm12,ZMMWORD ptr[rsi]	;src

	add rsi,r14
	add rdi,r14
	loop Resize_H_AVX512_Planar_32bits_loop_5

	vextractf32x8 ymm8,zmm0,1

	vhaddps ymm0,ymm0,ymm8

	vextractf128 xmm8,ymm0,1

	vhaddps xmm0,xmm0,xmm8

	vhaddps xmm0,xmm0,xmm0

	vhaddps xmm0,xmm0,xmm0

	vpextrd eax,xmm0,0
	mov dword ptr[r11],eax

Resize_H_AVX512_Planar_32bits_end:
	vmovdqa xmm12,XMMWORD ptr[rsp+96]
	vmovdqa xmm11,XMMWORD ptr[rsp+80]
	vmovdqa xmm10,XMMWORD ptr[rsp+64]
	vmovdqa xmm9,XMMWORD ptr[rsp+48]
	vmovdqa xmm8,XMMWORD ptr[rsp+32]
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,120
		
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

Resize_H_AVX512_Planar_32bits_ASM endp


end
