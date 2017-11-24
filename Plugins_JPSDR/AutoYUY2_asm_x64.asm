.data

align 16

uw_2 word 8 dup(2)
uw_3 word 8 dup(3)
uw_4 word 8 dup(4)
uw_5 word 8 dup(5)
uw_7 word 8 dup(7)

.code

;JPSDR_AutoYUY2_1 proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword
; src_y = rcx
; src_u = rdx
; src_v = r8
; dst = r9
	
JPSDR_AutoYUY2_1 proc public frame

w equ dword ptr[rbp+48]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	.endprolog


	mov rsi,rcx			; rsi=src_y
	mov rdi,r9			; rdi=dst
	mov rbx,r8			; rbx=srv_v
	xor rcx,rcx
	mov ecx,w
	cld

SSE2_1_a:
	mov al,byte ptr[rsi+1]		;al=y2
	mov ah,byte ptr[rbx]		;ah=v
	inc rbx
	shl eax,16
	lodsw						;al=y1 ah=y2
	mov ah,byte ptr[rdx]		;ah=u
	inc rdx
	stosd
	loop SSE2_1_a
	
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret

JPSDR_AutoYUY2_1 endp




;JPSDR_AutoYUY2_SSE2_1 proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword
; src_y = rcx
; src_u = rdx
; src_v = r8
; dst = r9
	
JPSDR_AutoYUY2_SSE2_1 proc public frame

w equ dword ptr[rbp+48]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	.endprolog

	mov r10,rcx
	xor rcx,rcx
	mov ecx,w
	xor rax,rax
	mov r11,16
	
	or ecx,ecx
	jz short _SSE2_1_c
	
_SSE2_1_a:
	movd xmm1,dword ptr[rdx+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r8+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	movq xmm0,qword ptr[r10+8*rax]		;00000000YYYYYYYY
	inc rax
	punpcklbw xmm0,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[r9],xmm0
	add r9,r11
	loop _SSE2_1_a
	
_SSE2_1_c:	
	pop rbp

	ret

JPSDR_AutoYUY2_SSE2_1 endp


;JPSDR_AutoYUY2_AVX_1 proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword
; src_y = rcx
; src_u = rdx
; src_v = r8
; dst = r9
	
JPSDR_AutoYUY2_AVX_1 proc public frame

w equ dword ptr[rbp+48]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	.endprolog

	mov r10,rcx
	xor rcx,rcx
	mov ecx,w
	xor rax,rax
	mov r11,16
	
	or ecx,ecx
	jz short _AVX_1_c
	
_AVX_1_a:
	vmovd xmm1,dword ptr[rdx+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r8+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vmovq xmm0,qword ptr[r10+8*rax]		;00000000YYYYYYYY
	inc rax
	vpunpcklbw xmm0,xmm0,xmm1  			;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[r9],xmm0
	add r9,r11
	loop _AVX_1_a
	
_AVX_1_c:	
	pop rbp

	ret

JPSDR_AutoYUY2_AVX_1 endp


;JPSDR_AutoYUY2_SSE2_1b proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword
; src_y = rcx
; src_u = rdx
; src_v = r8
; dst = r9
	
JPSDR_AutoYUY2_SSE2_1b proc public frame

w equ dword ptr[rbp+48]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push r12
	.pushreg r12	
	push r13
	.pushreg r13	
	.endprolog

	mov r10,rcx
	xor rcx,rcx
	mov ecx,w
	xor rax,rax
	mov r11,32
	mov r12,2
	mov r13,16
	
	shr ecx,1
	jz short _SSE2_1b_b
	
_SSE2_1b_a:
	movq xmm1,qword ptr[rdx+4*rax]		;00000000UUUUUUUU
	movq xmm0,qword ptr[r8+4*rax]		;00000000VVVVVVVV
	movdqa xmm2,XMMWORD ptr[r10+8*rax]		;YYYYYYYYYYYYYYYY
	punpcklbw xmm1,xmm0					;VUVUVUVUVUVUVUVU
	movdqa xmm3,xmm2
	add rax,r12
	punpcklbw xmm2,xmm1     			;VYUYVYUYVYUYVYUY
	punpckhbw xmm3,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[r9],xmm2
	movdqa XMMWORD ptr[r9+r13],xmm3
	add r9,r11
	loop _SSE2_1b_a
	
_SSE2_1b_b:
	mov ecx,w
	and ecx,1
	jz short _SSE2_1b_c
	
	movd xmm1,dword ptr[rdx+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r8+4*rax]		;000000000000VVVV
	movq xmm2,qword ptr[r10+8*rax]		;00000000YYYYYYYY
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	punpcklbw xmm2,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[r9],xmm2
	
_SSE2_1b_c:
	pop r13
	pop r12	
	pop rbp

	ret

JPSDR_AutoYUY2_SSE2_1b endp


;JPSDR_AutoYUY2_AVX_1b proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword
; src_y = rcx
; src_u = rdx
; src_v = r8
; dst = r9
	
JPSDR_AutoYUY2_AVX_1b proc public frame

w equ dword ptr[rbp+48]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push r12
	.pushreg r12	
	push r13
	.pushreg r13	
	.endprolog

	mov r10,rcx
	xor rcx,rcx
	mov ecx,w
	xor rax,rax
	mov r11,32
	mov r12,2
	mov r13,16
	
	shr ecx,1
	jz short _AVX_1b_b
	
_AVX_1b_a:
	vmovq xmm1,qword ptr[rdx+4*rax]		;00000000UUUUUUUU
	vmovq xmm0,qword ptr[r8+4*rax]		;00000000VVVVVVVV
	vmovdqa xmm2,XMMWORD ptr[r10+8*rax]		;YYYYYYYYYYYYYYYY
	vpunpcklbw xmm1,xmm1,xmm0				;VUVUVUVUVUVUVUVU
	add rax,r12
	vpunpcklbw xmm3,xmm2,xmm1  			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm4,xmm2,xmm1  			;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[r9],xmm3
	vmovdqa XMMWORD ptr[r9+r13],xmm4
	add r9,r11
	loop _AVX_1b_a
	
_AVX_1b_b:
	mov ecx,w
	and ecx,1
	jz short _AVX_1b_c
	
	vmovd xmm1,dword ptr[rdx+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r8+4*rax]		;000000000000VVVV
	vmovq xmm2,qword ptr[r10+8*rax]		;00000000YYYYYYYY
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm1   		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[r9],xmm2
	
_AVX_1b_c:
	pop r13
	pop r12	
	pop rbp

	ret

JPSDR_AutoYUY2_AVX_1b endp


;JPSDR_AutoYUY2_SSE2_2 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_SSE2_2 proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,32
	.allocstack 32
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov rbx,8
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_3
	movdqa xmm4,XMMWORD ptr uw_5
	
	xor rax,rax
	or ecx,ecx
	jz short _SSE2_2_c

_SSE2_2_a:
	movd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	movd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc rax
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[rdi],xmm0
	add rdi,rdx
	
	loop _SSE2_2_a
		
_SSE2_2_c:		
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_2 endp


;JPSDR_AutoYUY2_AVX_2 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_AVX_2 proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,32
	.allocstack 32
	vmovdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	vpxor xmm7,xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov rbx,8
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_3
	vmovdqa xmm4,XMMWORD ptr uw_5
	
	xor rax,rax
	or ecx,ecx
	jz short _AVX_2_c

_AVX_2_a:
	vmovd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	vmovd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	

	vpmullw xmm1,xmm1,xmm5
	vpmullw xmm2,xmm2,xmm4
	vpaddsw xmm1,xmm1,xmm6
	vmovq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc rax
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[rdi],xmm0
	add rdi,rdx
	
	loop _AVX_2_a
		
_AVX_2_c:		
	vmovdqu xmm7,XMMWORD ptr[rsp+16]
	vmovdqu xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_AVX_2 endp


;JPSDR_AutoYUY2_SSE2_2b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_SSE2_2b proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push r12
	.pushreg r12
	push r13
	.pushreg r13		
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,32
	mov r12,2
	mov r13,16
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_3
	movdqa xmm4,XMMWORD ptr uw_5
	
	xor rax,rax
	shr ecx,1
	jz _SSE2_2b_b

_SSE2_2b_a:
	movq xmm0,qword ptr[r10+4*rax]		;00000000UUUUUUUU
	movq xmm1,qword ptr[r9+4*rax]		;00000000VVVVVVVV
	punpcklbw xmm0,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm1,xmm7				;0V0V0V0V0V0V0V0V
	movq xmm2,qword ptr[r8+4*rax]		;00000000UUUUUUUU
	movq xmm3,qword ptr[r11+4*rax]		;00000000VVVVVVVV
	punpcklbw xmm2,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm3,xmm7				;0V0V0V0V0V0V0V0V
	
	pmullw xmm0,xmm5
	pmullw xmm1,xmm5	
	pmullw xmm2,xmm4
	pmullw xmm3,xmm4	
	paddsw xmm2,xmm6
	paddsw xmm3,xmm6
	paddsw xmm0,xmm2
	paddsw xmm1,xmm3
	movdqa xmm2,XMMWORD ptr[rsi+8*rax]		;YYYYYYYYYYYYYYYY
	psraw xmm0,3
	psraw xmm1,3
	packuswb xmm0,xmm7				;00000000UUUUUUUU
	packuswb xmm1,xmm7				;00000000VVVVVVVV
	movdqa xmm3,xmm2
	punpcklbw xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add rax,r12
	punpcklbw xmm2,xmm0				;VYUYVYUYVYUYVYUY
	punpckhbw xmm3,xmm0				;VYUYVYUYVYUYVYUY
	movdqa XMMWORD ptr[rdi],xmm2
	movdqa XMMWORD ptr[rdi+r13],xmm3
	add rdi,rdx
	
	dec ecx
	jnz _SSE2_2b_a

_SSE2_2b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_2b_c
	
	movd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	movd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc rax
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[rdi],xmm0
	
_SSE2_2b_c:		
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop r13
	pop r12	
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_2b endp


;JPSDR_AutoYUY2_AVX_2b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_AVX_2b proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push r12
	.pushreg r12
	push r13
	.pushreg r13		
	sub rsp,32
	.allocstack 32
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	vpxor xmm7,xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,32
	mov r12,2
	mov r13,16
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_3
	vmovdqa xmm4,XMMWORD ptr uw_5
	
	xor rax,rax
	shr ecx,1
	jz _AVX_2b_b

_AVX_2b_a:
	vmovq xmm0,qword ptr[r10+4*rax]		;00000000UUUUUUUU
	vmovq xmm1,qword ptr[r9+4*rax]		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm1,xmm1,xmm7			;0V0V0V0V0V0V0V0V
	vmovq xmm2,qword ptr[r8+4*rax]		;00000000UUUUUUUU
	vmovq xmm3,qword ptr[r11+4*rax]		;00000000VVVVVVVV
	vpunpcklbw xmm2,xmm2,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm3,xmm3,xmm7			;0V0V0V0V0V0V0V0V
	
	vpmullw xmm0,xmm0,xmm5
	vpmullw xmm1,xmm1,xmm5	
	vpmullw xmm2,xmm2,xmm4
	vpmullw xmm3,xmm3,xmm4	
	vpaddsw xmm2,xmm2,xmm6
	vpaddsw xmm3,xmm3,xmm6
	vpaddsw xmm0,xmm0,xmm2
	vpaddsw xmm1,xmm1,xmm3
	vmovdqa xmm2,XMMWORD ptr[rsi+8*rax]		;YYYYYYYYYYYYYYYY
	vpsraw xmm0,xmm0,3
	vpsraw xmm1,xmm1,3
	vpackuswb xmm0,xmm0,xmm7			;00000000UUUUUUUU
	vpackuswb xmm1,xmm1,xmm7			;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add rax,r12
	vpunpcklbw xmm3,xmm2,xmm0			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm2,xmm2,xmm0			;VYUYVYUYVYUYVYUY
	vmovdqa XMMWORD ptr[rdi],xmm3
	vmovdqa XMMWORD ptr[rdi+r13],xmm2
	add rdi,rdx
	
	dec ecx
	jnz _AVX_2b_a

_AVX_2b_b:	
	mov ecx,w
	and ecx,1
	jz short _AVX_2b_c
	
	vmovd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	vmovd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	

	vpmullw xmm1,xmm1,xmm5
	vpmullw xmm2,xmm2,xmm4
	vpaddsw xmm1,xmm1,xmm6
	vmovq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc rax
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[rdi],xmm0
	
_AVX_2b_c:		
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop r13
	pop r12	
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_AVX_2b endp


;JPSDR_AutoYUY2_SSE2_3 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_SSE2_3 proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,32
	.allocstack 32
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov rbx,8
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_7
	
	xor rax,rax
	or ecx,ecx
	jz short _SSE2_3_c

_SSE2_3_a:
	movd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	movd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc rax
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[rdi],xmm0
	add rdi,rdx
	
	loop _SSE2_3_a
	
_SSE2_3_c:		
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_3 endp


;JPSDR_AutoYUY2_AVX_3 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_AVX_3 proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,32
	.allocstack 32
	vmovdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	vpxor xmm7,xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov rbx,8
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_7
	
	xor rax,rax
	or ecx,ecx
	jz short _AVX_3_c

_AVX_3_a:
	vmovd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0				;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7				;0V0U0V0U0V0U0V0U
	vmovd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0				;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc rax
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7				;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[rdi],xmm0
	add rdi,rdx
	
	loop _AVX_3_a
	
_AVX_3_c:		
	vmovdqu xmm7,XMMWORD ptr[rsp+16]
	vmovdqu xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_AVX_3 endp


;JPSDR_AutoYUY2_SSE2_3b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_SSE2_3b proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push r12
	.pushreg r12
	push r13
	.pushreg r13	
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,32
	mov r12,2
	mov r13,16
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_7
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_3b_b

_SSE2_3b_a:
	movq xmm0,qword ptr[r10+4*rax]		;00000000UUUUUUUU
	movq xmm1,qword ptr[r9+4*rax]		;00000000VVVVVVVV
	punpcklbw xmm0,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm1,xmm7				;0V0V0V0V0V0V0V0V
	movq xmm2,qword ptr[r8+4*rax]		;00000000UUUUUUUU
	movq xmm3,qword ptr[r11+4*rax]		;00000000VVVVVVVV
	punpcklbw xmm2,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm3,xmm7				;0V0V0V0V0V0V0V0V
	movdqa xmm4,XMMWORD ptr[rsi+8*rax]		;YYYYYYYYYYYYYYYY
	
	pmullw xmm0,xmm5
	pmullw xmm1,xmm5	
	paddsw xmm2,xmm6
	paddsw xmm3,xmm6
	paddsw xmm0,xmm2
	paddsw xmm1,xmm3
	psraw xmm0,3
	psraw xmm1,3
	packuswb xmm0,xmm7				;00000000UUUUUUUU
	packuswb xmm1,xmm7				;00000000VVVVVVVV
	movdqa xmm2,xmm4
	punpcklbw xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add rax,r12
	punpcklbw xmm2,xmm0				;VYUYVYUYVYUYVYUY
	punpckhbw xmm4,xmm0				;VYUYVYUYVYUYVYUY
	movdqa XMMWORD ptr[rdi],xmm2
	movdqa XMMWORD ptr[rdi+r13],xmm4
	add rdi,rdx
	
	loop _SSE2_3b_a

_SSE2_3b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_3b_c
	
	movd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	movd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[rdi],xmm0
		
_SSE2_3b_c:		
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop r13
	pop r12	
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_3b endp


;JPSDR_AutoYUY2_AVX_3b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_AVX_3b proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push r12
	.pushreg r12
	push r13
	.pushreg r13	
	sub rsp,32
	.allocstack 32
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	vpxor xmm7,xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,32
	mov r12,2
	mov r13,16
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_7
	
	xor rax,rax
	shr ecx,1
	jz short _AVX_3b_b

_AVX_3b_a:
	vmovq xmm0,qword ptr[r10+4*rax]		;00000000UUUUUUUU
	vmovq xmm1,qword ptr[r9+4*rax]		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm1,xmm1,xmm7			;0V0V0V0V0V0V0V0V
	vmovq xmm2,qword ptr[r8+4*rax]		;00000000UUUUUUUU
	vmovq xmm3,qword ptr[r11+4*rax]		;00000000VVVVVVVV
	vpunpcklbw xmm2,xmm2,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm3,xmm3,xmm7			;0V0V0V0V0V0V0V0V
	vmovdqa xmm4,XMMWORD ptr[rsi+8*rax]		;YYYYYYYYYYYYYYYY
	
	vpmullw xmm0,xmm0,xmm5
	vpmullw xmm1,xmm1,xmm5	
	vpaddsw xmm2,xmm2,xmm6
	vpaddsw xmm3,xmm3,xmm6
	vpaddsw xmm0,xmm0,xmm2
	vpaddsw xmm1,xmm1,xmm3
	vpsraw xmm0,xmm0,3
	vpsraw xmm1,xmm1,3
	vpackuswb xmm0,xmm0,xmm7			;00000000UUUUUUUU
	vpackuswb xmm1,xmm1,xmm7			;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add rax,r12
	vpunpcklbw xmm2,xmm4,xmm0			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm4,xmm4,xmm0			;VYUYVYUYVYUYVYUY
	vmovdqa XMMWORD ptr[rdi],xmm2
	vmovdqa XMMWORD ptr[rdi+r13],xmm4
	add rdi,rdx
	
	loop _AVX_3b_a

_AVX_3b_b:	
	mov ecx,w
	and ecx,1
	jz short _AVX_3b_c
	
	vmovd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	vmovd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	

	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[rdi],xmm0
		
_AVX_3b_c:		
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop r13
	pop r12	
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_AVX_3b endp


;JPSDR_AutoYUY2_SSE2_4 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_SSE2_4 proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,32
	.allocstack 32
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov rbx,8
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_2
	movdqa xmm5,XMMWORD ptr uw_3
	
	xor rax,rax
	or ecx,ecx
	jz short _SSE2_4_c

_SSE2_4_a:
	movd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	movd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc rax
	psraw xmm1,2
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[rdi],xmm0
	add rdi,rdx
	
	loop _SSE2_4_a
		
_SSE2_4_c:		
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_4 endp


;JPSDR_AutoYUY2_AVX_4 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_AVX_4 proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,32
	.allocstack 32
	vmovdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	vpxor xmm7,xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov rbx,8
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_2
	vmovdqa xmm5,XMMWORD ptr uw_3
	
	xor rax,rax
	or ecx,ecx
	jz short _AVX_4_c

_AVX_4_a:
	vmovd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	vmovd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	

	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc rax
	vpsraw xmm1,xmm1,2
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[rdi],xmm0
	add rdi,rdx
	
	loop _AVX_4_a
		
_AVX_4_c:		
	vmovdqu xmm7,XMMWORD ptr[rsp+16]
	vmovdqu xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_AVX_4 endp


;JPSDR_AutoYUY2_SSE2_4b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_SSE2_4b proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push r12
	.pushreg r12
	push r13
	.pushreg r13	
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,32
	mov r12,2
	mov r13,16
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_2
	movdqa xmm5,XMMWORD ptr uw_3
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_4b_b

_SSE2_4b_a:
	movq xmm0,qword ptr[r10+4*rax]		;00000000UUUUUUUU
	movq xmm1,qword ptr[r9+4*rax]		;00000000VVVVVVVV
	punpcklbw xmm0,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm1,xmm7				;0V0V0V0V0V0V0V0V
	movq xmm2,qword ptr[r8+4*rax]		;00000000UUUUUUUU
	movq xmm3,qword ptr[r11+4*rax]		;00000000VVVVVVVV
	punpcklbw xmm2,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm3,xmm7				;0V0V0V0V0V0V0V0V
	movdqa xmm4,XMMWORD ptr[rsi+8*rax]		;YYYYYYYYYYYYYYYY
	
	pmullw xmm0,xmm5
	pmullw xmm1,xmm5	
	paddsw xmm2,xmm6
	paddsw xmm3,xmm6
	paddsw xmm0,xmm2
	paddsw xmm1,xmm3
	psraw xmm0,2
	psraw xmm1,2
	packuswb xmm0,xmm7				;00000000UUUUUUUU
	packuswb xmm1,xmm7				;00000000VVVVVVVV
	movdqa xmm2,xmm4
	punpcklbw xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add rax,r12
	punpcklbw xmm2,xmm0				;VYUYVYUYVYUYVYUY
	punpckhbw xmm4,xmm0				;VYUYVYUYVYUYVYUY
	movdqa XMMWORD ptr[rdi],xmm2
	movdqa XMMWORD ptr[rdi+r13],xmm4
	add rdi,rdx
	
	loop _SSE2_4b_a

_SSE2_4b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4b_c
	
	movd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	movd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	

	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc rax
	psraw xmm1,2
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[rdi],xmm0
		
_SSE2_4b_c:		
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop r13
	pop r12	
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_4b endp


;JPSDR_AutoYUY2_AVX_4b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword
; src_y = rcx
; src1_u = rdx
; src2_u = r8
; src1_v = r9

JPSDR_AutoYUY2_AVX_4b proc public frame

src2_v equ qword ptr[rbp+48]
dst equ qword ptr[rbp+56]
w equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push r12
	.pushreg r12
	push r13
	.pushreg r13	
	sub rsp,32
	.allocstack 32
	vmovdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	vmovdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	vpxor xmm7,xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,32
	mov r12,2
	mov r13,16
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_2
	vmovdqa xmm5,XMMWORD ptr uw_3
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_4b_b

_SSE2_4b_a:
	vmovq xmm0,qword ptr[r10+4*rax]		;00000000UUUUUUUU
	vmovq xmm1,qword ptr[r9+4*rax]		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm1,xmm1,xmm7			;0V0V0V0V0V0V0V0V
	vmovq xmm2,qword ptr[r8+4*rax]		;00000000UUUUUUUU
	vmovq xmm3,qword ptr[r11+4*rax]		;00000000VVVVVVVV
	vpunpcklbw xmm2,xmm2,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm3,xmm3,xmm7			;0V0V0V0V0V0V0V0V
	vmovdqa xmm4,XMMWORD ptr[rsi+8*rax]		;YYYYYYYYYYYYYYYY
	
	vpmullw xmm0,xmm0,xmm5
	vpmullw xmm1,xmm1,xmm5	
	vpaddsw xmm2,xmm2,xmm6
	vpaddsw xmm3,xmm3,xmm6
	vpaddsw xmm0,xmm0,xmm2
	vpaddsw xmm1,xmm1,xmm3
	vpsraw xmm0,xmm0,2
	vpsraw xmm1,xmm1,2
	vpackuswb xmm0,xmm0,xmm7			;00000000UUUUUUUU
	vpackuswb xmm1,xmm1,xmm7			;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add rax,r12
	vpunpcklbw xmm2,xmm4,xmm0			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm4,xmm4,xmm0			;VYUYVYUYVYUYVYUY
	vmovdqa XMMWORD ptr[rdi],xmm2
	vmovdqa XMMWORD ptr[rdi+r13],xmm4
	add rdi,rdx
	
	loop _SSE2_4b_a

_SSE2_4b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4b_c
	
	vmovd xmm1,dword ptr[r10+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r9+4*rax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	vmovd xmm2,dword ptr[r8+4*rax]		;000000000000UUUU
	vmovd xmm0,dword ptr[r11+4*rax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	

	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[rsi+8*rax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc rax
	vpsraw xmm1,xmm1,2
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[rdi],xmm0
		
_SSE2_4b_c:		
	vmovdqa xmm7,XMMWORD ptr[rsp+16]
	vmovdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32
	pop r13
	pop r12	
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_AVX_4b endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2 proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2 proc public frame

	.endprolog
		
	pcmpeqb xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	mov ecx,r9d
	shr ecx,1
	jz short SSE2_2_b
		
SSE2_2_a:
	movq xmm0,qword ptr[rdx+8*rax]
	movq xmm1,qword ptr[r10+8*rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movq qword ptr[r8+8*rax],xmm0
	inc rax
	loop SSE2_2_a
	
SSE2_2_b:
	and r9d,1
	jz short SSE2_2_c
	
	movd xmm0,dword ptr[rdx+8*rax]
	movd xmm1,dword ptr[r10+8*rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movd dword ptr[r8+8*rax],xmm0
	
SSE2_2_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2 endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2 proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2 proc public frame

	.endprolog
		
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	mov ecx,r9d
	shr ecx,1
	jz short AVX_2_b
		
AVX_2_a:
	vmovq xmm0,qword ptr[rdx+8*rax]
	vmovq xmm1,qword ptr[r10+8*rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovq qword ptr[r8+8*rax],xmm2
	inc rax
	loop AVX_2_a
	
AVX_2_b:
	and r9d,1
	jz short AVX_2_c
	
	vmovd xmm0,dword ptr[rdx+8*rax]
	vmovd xmm1,dword ptr[r10+8*rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovd dword ptr[r8+8*rax],xmm2
	
AVX_2_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2 endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b proc public frame

	.endprolog
		
	pcmpeqb xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short SSE2_2b_b
	
	mov r11,16
SSE2_2b_a:
	movdqa xmm0,XMMWORD ptr[rdx+rax]
	movdqa xmm1,XMMWORD ptr[r10+rax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movdqa XMMWORD ptr[r8+rax],xmm0
	add rax,r11
	loop SSE2_2b_a

SSE2_2b_b:
	mov ecx,r9d
	and ecx,1
	jz short SSE2_2b_c

	movq xmm0,qword ptr[rdx+rax]
	movq xmm1,qword ptr[r10+rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movq qword ptr[r8+rax],xmm0	
		
SSE2_2b_c:		
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b proc public frame

	.endprolog
		
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short AVX_2b_b
	
	mov r11,16
AVX_2b_a:
	vmovdqa xmm0,XMMWORD ptr[rdx+rax]
	vmovdqa xmm1,XMMWORD ptr[r10+rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovdqa XMMWORD ptr[r8+rax],xmm2
	add rax,r11
	loop AVX_2b_a

AVX_2b_b:
	mov ecx,r9d
	and ecx,1
	jz short AVX_2b_c

	vmovq xmm0,qword ptr[rdx+rax]
	vmovq xmm1,qword ptr[r10+rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovq qword ptr[r8+rax],xmm2
		
AVX_2b_c:		
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3 proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
	
JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3 proc public frame	
	
	.endprolog
	
	pcmpeqb xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short SSE2_3_b
	
SSE2_3_a:
	movq xmm0,qword ptr[r10+8*rax]
	movq xmm1,qword ptr[rdx+8*rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movq qword ptr[r8+8*rax],xmm1
	inc rax
	loop SSE2_3_a
	
SSE2_3_b:
	and r9d,1
	jz short SSE2_3_c

	movd xmm0,dword ptr[r10+8*rax]
	movd xmm1,dword ptr[rdx+8*rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movd dword ptr[r8+8*rax],xmm1
	
SSE2_3_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3 endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3 proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
	
JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3 proc public frame	
	
	.endprolog
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short AVX_3_b
	
AVX_3_a:
	vmovq xmm0,qword ptr[r10+8*rax]
	vmovq xmm1,qword ptr[rdx+8*rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovq qword ptr[r8+8*rax],xmm1
	inc rax
	loop AVX_3_a
	
AVX_3_b:
	and r9d,1
	jz short AVX_3_c

	vmovd xmm0,dword ptr[r10+8*rax]
	vmovd xmm1,dword ptr[rdx+8*rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovd dword ptr[r8+8*rax],xmm1
	
AVX_3_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3 endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
	
JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b proc public frame	
	
	.endprolog
	
	pcmpeqb xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short SSE2_3b_b
		
	mov r11,16
SSE2_3b_a:
	movdqa xmm0,XMMWORD ptr[r10+rax]
	movdqa xmm1,XMMWORD ptr[rdx+rax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movdqa XMMWORD ptr[r8+rax],xmm1
	add rax,r11
	loop SSE2_3b_a
	
SSE2_3b_b:	
	mov ecx,r9d
	and ecx,1
	jz short SSE2_3b_c
	
	movq xmm0,qword ptr[r10+rax]
	movq xmm1,qword ptr[rdx+rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movq qword ptr[r8+rax],xmm1	
	
SSE2_3b_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d
	
JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b proc public frame	
	
	.endprolog
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short AVX_3b_b
		
	mov r11,16
AVX_3b_a:
	vmovdqa xmm0,XMMWORD ptr[r10+rax]
	vmovdqa xmm1,XMMWORD ptr[rdx+rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovdqa XMMWORD ptr[r8+rax],xmm1
	add rax,r11
	loop AVX_3b_a
	
AVX_3b_b:	
	mov ecx,r9d
	and ecx,1
	jz short AVX_3b_c
	
	vmovq xmm0,qword ptr[r10+rax]
	vmovq xmm1,qword ptr[rdx+rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovq qword ptr[r8+rax],xmm1	
	
AVX_3b_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4 proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4 proc public frame

	.endprolog
		
	pcmpeqb xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short SSE2_4_b

SSE2_4_a:
	movq xmm0,qword ptr[r10+8*rax]
	movq xmm1,qword ptr[rdx+8*rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movq qword ptr[r8+8*rax],xmm0
	inc rax
	loop SSE2_4_a
	
SSE2_4_b:
	and r9d,1
	jz short SSE2_4_c
	
	movd xmm0,dword ptr[r10+8*rax]
	movd xmm1,dword ptr[rdx+8*rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movd dword ptr[r8+8*rax],xmm0
	
SSE2_4_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4 endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4 proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4 proc public frame

	.endprolog
		
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short AVX_4_b

AVX_4_a:
	vmovq xmm0,qword ptr[r10+8*rax]
	vmovq xmm1,qword ptr[rdx+8*rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovq qword ptr[r8+8*rax],xmm2
	inc rax
	loop AVX_4_a
	
AVX_4_b:
	and r9d,1
	jz short AVX_4_c
	
	vmovd xmm0,dword ptr[r10+8*rax]
	vmovd xmm1,dword ptr[rdx+8*rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovd dword ptr[r8+8*rax],xmm2
	
AVX_4_c:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4 endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b proc public frame

	.endprolog
		
	pcmpeqb xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short SSE2_4b_b
	
	mov r11,16
SSE2_4b_a:
	movdqa xmm0,XMMWORD ptr[r10+rax]
	movdqa xmm1,XMMWORD ptr[rdx+rax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movdqa XMMWORD ptr[r8+rax],xmm0
	add rax,r11
	loop SSE2_4b_a
	
SSE2_4b_b:
	mov ecx,r9d
	and ecx,1
	jz short SSE2_4b_c
	
	movq xmm0,qword ptr[r10+rax]
	movq xmm1,qword ptr[rdx+rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movq qword ptr[r8+rax],xmm0	
	
SSE2_4b_c:		
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b endp


;JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b proc src1:dword,src2:dword,dst:dword,w:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w = r9d

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b proc public frame

	.endprolog
		
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax
	
	mov ecx,r9d
	shr ecx,1
	jz short AVX_4b_b
	
	mov r11,16
AVX_4b_a:
	vmovdqa xmm0,XMMWORD ptr[r10+rax]
	vmovdqa xmm1,XMMWORD ptr[rdx+rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovdqa XMMWORD ptr[r8+rax],xmm2
	add rax,r11
	loop AVX_4b_a
	
AVX_4b_b:
	mov ecx,r9d
	and ecx,1
	jz short AVX_4b_c
	
	vmovq xmm0,qword ptr[r10+rax]
	vmovq xmm1,qword ptr[rdx+rax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovq qword ptr[r8+rax],xmm2
	
AVX_4b_c:		
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b endp

end





