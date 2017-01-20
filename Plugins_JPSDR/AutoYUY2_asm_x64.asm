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
	push rbx
	.pushreg rbx
	.endprolog

	mov r10,rcx
	xor rcx,rcx
	mov ecx,w
	xor rax,rax
	mov r11,16
	mov rbx,8
	
	shr ecx,1
	jz short _SSE2_1_b
	
_SSE2_1_a:
	movd xmm1,dword ptr[rdx+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r8+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	movq xmm0,qword ptr[r10+8*rax]		;00000000YYYYYYYY
	inc rax
	punpcklbw xmm0,xmm1     			;VYUYVYUYVYUYVYUY
	
	movq qword ptr[r9],xmm0
	psrldq xmm0,8
	movq qword ptr[r9+rbx],xmm0
	add r9,r11
	loop _SSE2_1_a
	
_SSE2_1_b:
	mov ecx,w
	and ecx,1
	jz short _SSE2_1_c
	
	movzx ecx,word ptr[rdx+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0					;000000000000VUVU
	movd xmm0,dword ptr[r10+8*rax]		;000000000000YYYY
	punpcklbw xmm0,xmm1     			;00000000VYUYVYUY
	movq qword ptr[r9],xmm0		
	
_SSE2_1_c:	
	pop rbx
	pop rbp

	ret

JPSDR_AutoYUY2_SSE2_1 endp



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
	.endprolog

	mov r10,rcx
	xor rcx,rcx
	mov ecx,w
	xor rax,rax
	mov r11,16
	
	shr ecx,1
	jz short _SSE2_1b_b
	
_SSE2_1b_a:
	movd xmm1,dword ptr[rdx+4*rax]		;000000000000UUUU
	movd xmm0,dword ptr[r8+4*rax]		;000000000000VVVV
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	movq xmm0,qword ptr[r10+8*rax]		;00000000YYYYYYYY
	inc rax
	punpcklbw xmm0,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa oword ptr[r9],xmm0
	add r9,r11
	loop _SSE2_1b_a
	
_SSE2_1b_b:
	mov ecx,w
	and ecx,1
	jz short _SSE2_1b_c
	
	movzx ecx,word ptr[rdx+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0					;000000000000VUVU
	movd xmm0,dword ptr[r10+8*rax]		;000000000000YYYY
	punpcklbw xmm0,xmm1     			;00000000VYUYVYUY
	movq qword ptr[r9],xmm0		
	
_SSE2_1b_c:	
	pop rbp

	ret

JPSDR_AutoYUY2_SSE2_1b endp



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
	movdqu [rsp],xmm6
	.savexmm128 xmm6,0
	movdqu [rsp+16],xmm7
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
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_3
	movdqa xmm4,oword ptr uw_5
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_2_b

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
	
	movq qword ptr[rdi],xmm0
	psrldq xmm0,8
	movq qword ptr[rdi+rbx],xmm0
	add rdi,rdx
	
	loop _SSE2_2_a

_SSE2_2_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_2_c
	
	movzx ecx,word ptr[r10+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r9+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[r11+4*rax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movd xmm0,dword ptr[rsi+8*rax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[rdi],xmm0
		
_SSE2_2_c:		
	movdqu xmm7,[rsp+16]
	movdqu xmm6,[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_2 endp


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
	sub rsp,32
	.allocstack 32
	movdqa [rsp],xmm6
	.savexmm128 xmm6,0
	movdqa [rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_3
	movdqa xmm4,oword ptr uw_5
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_2b_b

_SSE2_2b_a:
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
	
	movdqa oword ptr[rdi],xmm0
	add rdi,rdx
	
	loop _SSE2_2b_a

_SSE2_2b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_2b_c
	
	movzx ecx,word ptr[r10+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r9+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[r11+4*rax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movd xmm0,dword ptr[rsi+8*rax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[rdi],xmm0
		
_SSE2_2b_c:		
	movdqa xmm7,[rsp+16]
	movdqa xmm6,[rsp]
	add rsp,32
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_2b endp



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
	movdqu [rsp],xmm6
	.savexmm128 xmm6,0
	movdqu [rsp+16],xmm7
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
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_7
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_3_b

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
	
	movq qword ptr[rdi],xmm0
	psrldq xmm0,8
	movq qword ptr[rdi+rbx],xmm0
	add rdi,rdx
	
	loop _SSE2_3_a

_SSE2_3_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_3_c
	
	movzx ecx,word ptr[r10+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r9+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[r11+4*rax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[rsi+8*rax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[rdi],xmm0
		
_SSE2_3_c:		
	movdqu xmm7,[rsp+16]
	movdqu xmm6,[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_3 endp



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
	sub rsp,32
	.allocstack 32
	movdqa [rsp],xmm6
	.savexmm128 xmm6,0
	movdqa [rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_7
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_3b_b

_SSE2_3b_a:
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
	
	movdqa oword ptr[rdi],xmm0
	add rdi,rdx
	
	loop _SSE2_3b_a

_SSE2_3b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_3b_c
	
	movzx ecx,word ptr[r10+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r9+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[r11+4*rax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[rsi+8*rax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[rdi],xmm0
		
_SSE2_3b_c:		
	movdqa xmm7,[rsp+16]
	movdqa xmm6,[rsp]
	add rsp,32
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_3b endp



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
	movdqu [rsp],xmm6
	.savexmm128 xmm6,0
	movdqu [rsp+16],xmm7
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
	
	movdqa xmm6,oword ptr uw_2
	movdqa xmm5,oword ptr uw_3
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_4_b

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
	
	movq qword ptr[rdi],xmm0
	psrldq xmm0,8
	movq qword ptr[rdi+rbx],xmm0
	add rdi,rdx
	
	loop _SSE2_4_a

_SSE2_4_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4_c
	
	movzx ecx,word ptr[r10+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r9+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[r11+4*rax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[rsi+8*rax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,2
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[rdi],xmm0
		
_SSE2_4_c:		
	movdqu xmm7,[rsp+16]
	movdqu xmm6,[rsp]
	add rsp,32
	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_4 endp



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
	sub rsp,32
	.allocstack 32
	movdqa [rsp],xmm6
	.savexmm128 xmm6,0
	movdqa [rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog
	
	pxor xmm7,xmm7
	
	mov rsi,rcx		; rsi = src_y
	mov r10,rdx		; r10=src1_u
	mov rdi,dst
	mov r11,src2_v
	xor rcx,rcx
	mov rdx,16
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_2
	movdqa xmm5,oword ptr uw_3
	
	xor rax,rax
	shr ecx,1
	jz short _SSE2_4b_b

_SSE2_4b_a:
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
	
	movdqa oword ptr[rdi],xmm0
	add rdi,rdx
	
	loop _SSE2_4b_a

_SSE2_4b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4b_c
	
	movzx ecx,word ptr[r10+4*rax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[r9+4*rax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	movzx ecx,word ptr[r8+4*rax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[r11+4*rax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[rsi+8*rax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,2
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[rdi],xmm0
		
_SSE2_4b_c:		
	movdqa xmm7,[rsp+16]
	movdqa xmm6,[rsp]
	add rsp,32
	pop rdi
	pop rsi
	pop rbp

	ret	

JPSDR_AutoYUY2_SSE2_4b endp


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
	shr ecx,2
	jz short SSE2_2b_b
	
	mov r11,16
SSE2_2b_a:
	movdqa xmm0,oword ptr[rdx+rax]
	movdqa xmm1,oword ptr[r10+rax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movdqa oword ptr[r8+rax],xmm0
	add rax,r11
	loop SSE2_2b_a

SSE2_2b_b:
	mov ecx,r9d
	and ecx,3
	jz short SSE2_2b_d
	and ecx,2
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
		
	and r9d,1
	jz short SSE2_2b_d
	add rax,8
		
SSE2_2b_c:
	movd xmm0,dword ptr[rdx+rax]
	movd xmm1,dword ptr[r10+rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movd dword ptr[r8+rax],xmm0	
		
SSE2_2b_d:		
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b endp


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
	shr ecx,2
	jz short SSE2_3b_b
		
	mov r11,16
SSE2_3b_a:
	movdqa xmm0,oword ptr[r10+rax]
	movdqa xmm1,oword ptr[rdx+rax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movdqa oword ptr[r8+rax],xmm1
	add rax,r11
	loop SSE2_3b_a
	
SSE2_3b_b:	
	mov ecx,r9d
	and ecx,3
	jz short SSE2_3b_d
	and ecx,2
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
	
	and r9d,1
	jz short SSE2_3b_d
	add rax,8
	
SSE2_3b_c:	
	movd xmm0,dword ptr[r10+rax]
	movd xmm1,dword ptr[rdx+rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movd dword ptr[r8+rax],xmm1	
	
SSE2_3b_d:	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b endp



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
	shr ecx,2
	jz short SSE2_4b_b
	
	mov r11,16
SSE2_4b_a:
	movdqa xmm0,oword ptr[r10+rax]
	movdqa xmm1,oword ptr[rdx+rax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movdqa oword ptr[r8+rax],xmm0
	add rax,r11
	loop SSE2_4b_a
	
SSE2_4b_b:
	mov ecx,r9d
	and ecx,3
	jz short SSE2_4b_d
	and ecx,2
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
	
	and r9d,1
	jz short SSE2_4b_d
	add rax,8
	
SSE2_4b_c:	
	movd xmm0,dword ptr[r10+rax]
	movd xmm1,dword ptr[rdx+rax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movd dword ptr[r8+rax],xmm0
		
SSE2_4b_d:		
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b endp


end





