.586
.model flat,c

.data

align 16

uw_2 word 8 dup(2)
uw_3 word 8 dup(3)
uw_4 word 8 dup(4)
uw_5 word 8 dup(5)
uw_7 word 8 dup(7)

.code

JPSDR_AutoYUY2_1 proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_1

	push esi
	push edi
	push ebx

	mov esi,src_y
	mov ebx,src_u
	mov edx,src_v
	mov edi,dst
	mov ecx,w
	cld

_SSE2_0_a:
	mov al,byte ptr[esi+1]		;al=y2
	mov ah,byte ptr[edx]		;ah=v
	inc edx
	shl eax,16
	lodsw						;al=y1 ah=y2
	mov ah,byte ptr[ebx]		;ah=u
	inc ebx
	stosd
	loop _SSE2_0_a
	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_1 endp


.xmm


JPSDR_AutoYUY2_SSE2_1 proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_1

	push esi
	push edi
	push ebx

	mov esi,src_y
	mov ebx,src_u
	mov edx,src_v
	mov edi,dst
	mov ecx,w
	xor eax,eax
	
	shr ecx,1
	jz short _SSE2_1_b
	
_SSE2_1_a:
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	inc eax
	punpcklbw xmm0,xmm1     			;VYUYVYUYVYUYVYUY
	
	movq qword ptr[edi],xmm0
	psrldq xmm0,8
	movq qword ptr[edi+8],xmm0
	add edi,16
	loop _SSE2_1_a
	
_SSE2_1_b:
	mov ecx,w
	and ecx,1
	jz short _SSE2_1_c
	
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0					;000000000000VUVU
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	punpcklbw xmm0,xmm1     			;00000000VYUYVYUY
	movq qword ptr[edi],xmm0		
	
_SSE2_1_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_1 endp


JPSDR_AutoYUY2_SSE2_1b proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_1b

	push esi
	push edi
	push ebx

	mov esi,src_y
	mov ebx,src_u
	mov edx,src_v
	mov edi,dst
	mov ecx,w
	xor eax,eax
	
	shr ecx,1
	jz short _SSE2_1b_b
	
_SSE2_1b_a:
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	inc eax
	punpcklbw xmm0,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa oword ptr[edi],xmm0
	add edi,16
	loop _SSE2_1b_a
	
_SSE2_1b_b:
	mov ecx,w
	and ecx,1
	jz short _SSE2_1b_c
	
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0					;000000000000VUVU
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	punpcklbw xmm0,xmm1     			;00000000VYUYVYUY
	movq qword ptr[edi],xmm0		
	
_SSE2_1b_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_1b endp


JPSDR_AutoYUY2_SSE2_2 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_2

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_3
	movdqa xmm4,oword ptr uw_5
	
	xor eax,eax
	shr ecx,1
	jz short _SSE2_2_b

_SSE2_2_a:
	mov ebx,src1_u
	mov edx,src1_v
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	
	
	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc eax
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movq qword ptr[edi],xmm0
	psrldq xmm0,8
	movq qword ptr[edi+8],xmm0
	add edi,16
	
	loop _SSE2_2_a
	
_SSE2_2_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_2_c
	
	mov ebx,src1_u
	mov edx,src1_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[edi],xmm0
		
_SSE2_2_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_2 endp


JPSDR_AutoYUY2_SSE2_2b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_2b

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_3
	movdqa xmm4,oword ptr uw_5
	
	xor eax,eax
	shr ecx,1
	jz short _SSE2_2b_b

_SSE2_2b_a:
	mov ebx,src1_u
	mov edx,src1_v
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	
	
	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	inc eax
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa oword ptr[edi],xmm0
	add edi,16
	
	loop _SSE2_2b_a
	
_SSE2_2b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_2b_c
	
	mov ebx,src1_u
	mov edx,src1_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	pmullw xmm2,xmm4
	paddsw xmm1,xmm6
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[edi],xmm0
		
_SSE2_2b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_2b endp


JPSDR_AutoYUY2_SSE2_3 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_3

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_7
	
	xor eax,eax
	shr ecx,1
	jz short _SSE2_3_b

_SSE2_3_a:
	mov ebx,src1_u
	mov edx,src1_v
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc eax
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movq qword ptr[edi],xmm0
	psrldq xmm0,8
	movq qword ptr[edi+8],xmm0
	add edi,16
	
	loop _SSE2_3_a
	
_SSE2_3_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_3_c
	
	mov ebx,src1_u
	mov edx,src1_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[edi],xmm0
		
_SSE2_3_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_3 endp


JPSDR_AutoYUY2_SSE2_3b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_3b

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_4
	movdqa xmm5,oword ptr uw_7
	
	xor eax,eax
	shr ecx,1
	jz short _SSE2_3b_b

_SSE2_3b_a:
	mov ebx,src1_u
	mov edx,src1_v
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc eax
	psraw xmm1,3
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa oword ptr[edi],xmm0
	add edi,16
	
	loop _SSE2_3b_a
	
_SSE2_3b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_3b_c
	
	mov ebx,src1_u
	mov edx,src1_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,3
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[edi],xmm0
		
_SSE2_3b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_3b endp



JPSDR_AutoYUY2_SSE2_4 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_4

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_2
	movdqa xmm5,oword ptr uw_3
	
	xor eax,eax
	shr ecx,1
	jz short _SSE2_4_b

_SSE2_4_a:
	mov ebx,src1_u
	mov edx,src1_v
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc eax
	psraw xmm1,2
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movq qword ptr[edi],xmm0
	psrldq xmm0,8
	movq qword ptr[edi+8],xmm0
	add edi,16
	
	loop _SSE2_4_a
	
_SSE2_4_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4_c
	
	mov ebx,src1_u
	mov edx,src1_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,2
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[edi],xmm0
		
_SSE2_4_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_4 endp



JPSDR_AutoYUY2_SSE2_4b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_4b

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,oword ptr uw_2
	movdqa xmm5,oword ptr uw_3
	
	xor eax,eax
	shr ecx,1
	jz short _SSE2_4b_b

_SSE2_4b_a:
	mov ebx,src1_u
	mov edx,src1_v
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0				;00000000VUVUVUVU
	punpcklbw xmm1,xmm7				;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm2,xmm0				;00000000VUVUVUVU
	punpcklbw xmm2,xmm7				;0V0U0V0U0V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	paddsw xmm1,xmm2
	inc eax
	psraw xmm1,2
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa oword ptr[edi],xmm0
	add edi,16
	
	loop _SSE2_4b_a
	
_SSE2_4b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4b_c
	
	mov ebx,src1_u
	mov edx,src1_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm1,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0
	punpcklbw xmm1,xmm0				;000000000000VUVU
	punpcklbw xmm1,xmm7				;000000000V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	movzx ecx,word ptr[ebx+4*eax]
	pinsrw xmm2,ecx,0
	movzx ecx,word ptr[edx+4*eax]
	pinsrw xmm0,ecx,0	
	punpcklbw xmm2,xmm0				;000000000000VUVU
	punpcklbw xmm2,xmm7				;000000000V0U0V0U	
	
	pmullw xmm1,xmm5
	paddsw xmm2,xmm6
	movd xmm0,dword ptr[esi+8*eax]		;000000000000YYYY
	paddsw xmm1,xmm2
	psraw xmm1,2
	packuswb xmm1,xmm7				;000000000000VUVU
	punpcklbw xmm0,xmm1     		;00000000VYUYVYUY
	
	movq qword ptr[edi],xmm0
		
_SSE2_4b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_4b endp


JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2 proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2

	push esi
	push edi
	
	pcmpeqb xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short SSE2_2_b

SSE2_2_a:	
	movq xmm0,qword ptr[edx+8*eax]
	movq xmm1,qword ptr[esi+8*eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movq qword ptr[edi+8*eax],xmm0
	inc eax
	loop SSE2_2_a

SSE2_2_b:
	mov ecx,w
	and ecx,1
	jz short SSE2_2_c
	
	movd xmm0,dword ptr[edx+8*eax]
	movd xmm1,dword ptr[esi+8*eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movd dword ptr[edi+8*eax],xmm0
	
SSE2_2_c:	
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2 endp



JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b

	push esi
	push edi
	push ebx
	
	pcmpeqb xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,2
	jz short SSE2_2b_b

	mov ebx,16
SSE2_2b_a:	
	movdqa xmm0,oword ptr[edx+eax]
	movdqa xmm1,oword ptr[esi+eax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movdqa oword ptr[edi+eax],xmm0
	add eax,ebx
	loop SSE2_2b_a
	
SSE2_2b_b:
	mov ecx,w
	and ecx,3
	jz short SSE2_2b_d
	and ecx,2
	jz short SSE2_2b_c
	
	movq xmm0,qword ptr[edx+eax]
	movq xmm1,qword ptr[esi+eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movq qword ptr[edi+eax],xmm0
	
	mov ecx,w
	and ecx,1
	jz short SSE2_2b_d
	add eax,8
	
SSE2_2b_c:	
	movd xmm0,dword ptr[edx+eax]
	movd xmm1,dword ptr[esi+eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movd dword ptr[edi+eax],xmm0
	
SSE2_2b_d:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b endp



JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3 proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3
	
	push esi
	push edi
	
	pcmpeqb xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short SSE2_3_b
	
SSE2_3_a:
	movq xmm0,qword ptr[esi+8*eax]
	movq xmm1,qword ptr[edx+8*eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movq qword ptr[edi+8*eax],xmm1
	inc eax
	loop SSE2_3_a
	
SSE2_3_b:
	mov ecx,w
	and ecx,1
	jz short SSE2_3_c
	
	movd xmm0,dword ptr[esi+8*eax]
	movd xmm1,dword ptr[edx+8*eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movd dword ptr[edi+8*eax],xmm1	

SSE2_3_c:
	pop edi
	pop esi
	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3 endp



JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b
	
	push esi
	push edi
	push ebx
	
	pcmpeqb xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,2
	jz short SSE2_3b_b
	
	mov ebx,16
SSE2_3b_a:
	movdqa xmm0,oword ptr[esi+eax]
	movdqa xmm1,oword ptr[edx+eax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movdqa oword ptr[edi+eax],xmm1
	add eax,ebx
	loop SSE2_3b_a
	
SSE2_3b_b:	
	mov ecx,w
	and ecx,3
	jz short SSE2_3b_d
	and ecx,2
	jz short SSE2_3b_c
	
	movq xmm0,qword ptr[esi+eax]
	movq xmm1,qword ptr[edx+eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movq qword ptr[edi+eax],xmm1
	
	mov ecx,w
	and ecx,1
	jz short SSE2_3b_d
	add eax,8
	
SSE2_3b_c:	
	movd xmm0,dword ptr[esi+eax]
	movd xmm1,dword ptr[edx+eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movd dword ptr[edi+eax],xmm1
	
SSE2_3b_d:	
	pop ebx
	pop edi
	pop esi
	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b endp



JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4 proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4
	
	push esi
	push edi
	
	pcmpeqb xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short SSE2_4_b
	
SSE2_4_a:
	movq xmm0,qword ptr[esi+8*eax]
	movq xmm1,qword ptr[edx+8*eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movq qword ptr[edi+8*eax],xmm0
	inc eax
	loop SSE2_4_a
	
SSE2_4_b:
	mov ecx,w
	and ecx,1
	jz short SSE2_4_c

	movd xmm0,dword ptr[esi+8*eax]
	movd xmm1,dword ptr[edx+8*eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movd dword ptr[edi+8*eax],xmm0
	
SSE2_4_c:
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4 endp



JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b
	
	push esi
	push edi
	push ebx
	
	pcmpeqb xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,2
	jz short SSE2_4b_b
	
	mov ebx,16
SSE2_4b_a:
	movdqa xmm0,oword ptr[esi+eax]
	movdqa xmm1,oword ptr[edx+eax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movdqa oword ptr[edi+eax],xmm0
	add eax,ebx
	loop SSE2_4b_a
	
SSE2_4b_b:
	mov ecx,w
	and ecx,3
	jz short SSE2_4b_d
	and ecx,2
	jz short SSE2_4b_c

	movq xmm0,qword ptr[esi+eax]
	movq xmm1,qword ptr[edx+eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movq qword ptr[edi+eax],xmm0
	
	mov ecx,w
	and ecx,1
	jz short SSE2_4b_d
	add eax,8
	
SSE2_4b_c:	
	movd xmm0,dword ptr[esi+eax]
	movd xmm1,dword ptr[edx+eax]
	movq xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movd dword ptr[edi+eax],xmm0
	
SSE2_4b_d:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b endp


end





