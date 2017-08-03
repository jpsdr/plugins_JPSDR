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
	
	or ecx,ecx
	jz short _SSE2_1_c
	
_SSE2_1_a:
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	movq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	inc eax
	punpcklbw xmm0,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[edi],xmm0
	add edi,16
	loop _SSE2_1_a
	
_SSE2_1_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_1 endp


JPSDR_AutoYUY2_AVX_1 proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_1

	push esi
	push edi
	push ebx

	mov esi,src_y
	mov ebx,src_u
	mov edx,src_v
	mov edi,dst
	mov ecx,w
	xor eax,eax
	
	or ecx,ecx
	jz short _AVX_1_c
	
_AVX_1_a:
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	inc eax
	vpunpcklbw xmm0,xmm0,xmm1    		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0
	add edi,16
	loop _AVX_1_a
	
_AVX_1_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_1 endp


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
	movq xmm1,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm0,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	movdqa xmm2,XMMWORD ptr[esi+8*eax]	;YYYYYYYYYYYYYYYY	
	punpcklbw xmm1,xmm0					;VUVUVUVUVUVUVUVU
	movdqa xmm3,xmm2
	add eax,2
	punpcklbw xmm2,xmm1     			;VYUYVYUYVYUYVYUY
	punpckhbw xmm3,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[edi],xmm2
	movdqa XMMWORD ptr[edi+16],xmm3
	add edi,32
	loop _SSE2_1b_a
	
_SSE2_1b_b:
	mov ecx,w
	and ecx,1
	jz short _SSE2_1b_c
	
	movd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	movd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	movq xmm2,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	punpcklbw xmm1,xmm0					;00000000VUVUVUVU
	punpcklbw xmm2,xmm1     			;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[edi],xmm2
	
_SSE2_1b_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_1b endp


JPSDR_AutoYUY2_AVX_1b proc src_y:dword,src_u:dword,src_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_1b

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
	jz short _AVX_1b_b
	
_AVX_1b_a:
	vmovq xmm1,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm0,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	vmovdqa xmm2,XMMWORD ptr[esi+8*eax]	;YYYYYYYYYYYYYYYY	
	vpunpcklbw xmm1,xmm1,xmm0			;VUVUVUVUVUVUVUVU
	add eax,2
	vpunpcklbw xmm3,xmm2,xmm1  			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm4,xmm2,xmm1  			;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm3
	vmovdqa XMMWORD ptr[edi+16],xmm4
	add edi,32
	loop _AVX_1b_a
	
_AVX_1b_b:
	mov ecx,w
	and ecx,1
	jz short _AVX_1b_c
	
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vmovq xmm2,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm1  			;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm2
	
_AVX_1b_c:	
	pop ebx
	pop edi
	pop esi

	ret
	
JPSDR_AutoYUY2_AVX_1b endp


JPSDR_AutoYUY2_SSE2_2 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_2

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_3
	movdqa xmm4,XMMWORD ptr uw_5
	
	xor eax,eax
	or ecx,ecx
	jz short _SSE2_2_c

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
	
	movdqa XMMWORD ptr[edi],xmm0
	add edi,16
	
	loop _SSE2_2_a
		
_SSE2_2_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_2 endp


JPSDR_AutoYUY2_AVX_2 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_2

	push esi
	push edi
	push ebx
	
	vpxor xmm7,xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_3
	vmovdqa xmm4,XMMWORD ptr uw_5
	
	xor eax,eax
	or ecx,ecx
	jz short _AVX_2_c

_AVX_2_a:
	mov ebx,src1_u
	mov edx,src1_v
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	vmovd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	
	
	vpmullw xmm1,xmm1,xmm5
	vpmullw xmm2,xmm2,xmm4
	vpaddsw xmm1,xmm1,xmm6
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc eax
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7		;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1  		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0
	add edi,16
	
	loop _AVX_2_a
		
_AVX_2_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_2 endp


JPSDR_AutoYUY2_SSE2_2b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_2b

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_3
	movdqa xmm4,XMMWORD ptr uw_5
	
	xor eax,eax
	shr ecx,1
	jz _SSE2_2b_b

_SSE2_2b_a:
	mov ebx,src1_u
	mov edx,src1_v
	movq xmm0,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm1,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	punpcklbw xmm0,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm1,xmm7				;0V0V0V0V0V0V0V0V
	mov ebx,src2_u
	mov edx,src2_v
	movq xmm2,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm3,qword ptr[edx+4*eax]		;00000000VVVVVVVV
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
	movdqa xmm2,XMMWORD ptr[esi+8*eax]		;YYYYYYYYYYYYYYYY
	psraw xmm0,3
	psraw xmm1,3
	packuswb xmm0,xmm7				;00000000UUUUUUUU
	packuswb xmm1,xmm7				;00000000VVVVVVVV
	movdqa xmm3,xmm2
	punpcklbw xmm0,xmm1     		;VUVUVUVUVUVUVUVU
	add eax,2
	punpcklbw xmm2,xmm0				;VYUYVYUYVYUYVYUY
	punpckhbw xmm3,xmm0				;VYUYVYUYVYUYVYUY
	movdqa XMMWORD ptr[edi],xmm2
	movdqa XMMWORD ptr[edi+16],xmm3
	add edi,32
	
	dec ecx
	jnz _SSE2_2b_a
	
_SSE2_2b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_2b_c	
	
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
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[edi],xmm0	
		
_SSE2_2b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_2b endp


JPSDR_AutoYUY2_AVX_2b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_2b

	push esi
	push edi
	push ebx
	
	vpxor xmm7,xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_3
	vmovdqa xmm4,XMMWORD ptr uw_5
	
	xor eax,eax
	shr ecx,1
	jz _AVX_2b_b

_AVX_2b_a:
	mov ebx,src1_u
	mov edx,src1_v
	vmovq xmm0,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm1,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm1,xmm1,xmm7			;0V0V0V0V0V0V0V0V
	mov ebx,src2_u
	mov edx,src2_v
	vmovq xmm2,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm3,qword ptr[edx+4*eax]		;00000000VVVVVVVV
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
	vmovdqa xmm2,XMMWORD ptr[esi+8*eax]		;YYYYYYYYYYYYYYYY
	vpsraw xmm0,xmm0,3
	vpsraw xmm1,xmm1,3
	vpackuswb xmm0,xmm0,xmm7		;00000000UUUUUUUU
	vpackuswb xmm1,xmm1,xmm7		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm1  		;VUVUVUVUVUVUVUVU
	add eax,2
	vpunpcklbw xmm3,xmm2,xmm0		;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm2,xmm2,xmm0		;VYUYVYUYVYUYVYUY
	vmovdqa XMMWORD ptr[edi],xmm3
	vmovdqa XMMWORD ptr[edi+16],xmm2
	add edi,32
	
	dec ecx
	jnz _AVX_2b_a
	
_AVX_2b_b:	
	mov ecx,w
	and ecx,1
	jz short _AVX_2b_c	
	
	mov ebx,src1_u
	mov edx,src1_v
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	vmovd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	
	
	vpmullw xmm1,xmm1,xmm5
	vpmullw xmm2,xmm2,xmm4
	vpaddsw xmm1,xmm1,xmm6
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1    		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0	
		
_AVX_2b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_2b endp


JPSDR_AutoYUY2_SSE2_3 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_3

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_7
	
	xor eax,eax
	or ecx,ecx
	jz short _SSE2_3_c

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
	
	movdqa XMMWORD ptr[edi],xmm0
	add edi,16
	
	loop _SSE2_3_a
		
_SSE2_3_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_3 endp


JPSDR_AutoYUY2_AVX_3 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_3

	push esi
	push edi
	push ebx
	
	vpxor xmm7,xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_7
	
	xor eax,eax
	or ecx,ecx
	jz short _AVX_3_c

_AVX_3_a:
	mov ebx,src1_u
	mov edx,src1_v
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	vmovd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	
	
	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc eax
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0
	add edi,16
	
	loop _AVX_3_a
		
_AVX_3_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_3 endp


JPSDR_AutoYUY2_SSE2_3b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_3b

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_4
	movdqa xmm5,XMMWORD ptr uw_7
	
	xor eax,eax
	shr ecx,1
	jz _SSE2_3b_b

_SSE2_3b_a:
	mov ebx,src1_u
	mov edx,src1_v
	movq xmm0,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm1,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	punpcklbw xmm0,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm1,xmm7				;0V0V0V0V0V0V0V0V
	mov ebx,src2_u
	mov edx,src2_v
	movq xmm2,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm3,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	punpcklbw xmm2,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm3,xmm7				;0V0V0V0V0V0V0V0V
	movdqa xmm4,XMMWORD ptr[esi+8*eax]		;YYYYYYYYYYYYYYYY
	
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
	add eax,2
	punpcklbw xmm2,xmm0				;VYUYVYUYVYUYVYUY
	punpckhbw xmm4,xmm0				;VYUYVYUYVYUYVYUY
	movdqa XMMWORD ptr[edi],xmm2
	movdqa XMMWORD ptr[edi+16],xmm4
	add edi,32
	
	loop _SSE2_3b_a
	
_SSE2_3b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_3b_c
	
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
	
	movdqa XMMWORD ptr[edi],xmm0
	
_SSE2_3b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_3b endp


JPSDR_AutoYUY2_AVX_3b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_3b

	push esi
	push edi
	push ebx
	
	vpxor xmm7,xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_4
	vmovdqa xmm5,XMMWORD ptr uw_7
	
	xor eax,eax
	shr ecx,1
	jz _AVX_3b_b

_AVX_3b_a:
	mov ebx,src1_u
	mov edx,src1_v
	vmovq xmm0,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm1,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm1,xmm1,xmm7			;0V0V0V0V0V0V0V0V
	mov ebx,src2_u
	mov edx,src2_v
	vmovq xmm2,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm3,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	vpunpcklbw xmm2,xmm2,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm3,xmm3,xmm7			;0V0V0V0V0V0V0V0V
	vmovdqa xmm4,XMMWORD ptr[esi+8*eax]		;YYYYYYYYYYYYYYYY
	
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
	add eax,2
	vpunpcklbw xmm4,xmm2,xmm0			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm2,xmm2,xmm0			;VYUYVYUYVYUYVYUY
	vmovdqa XMMWORD ptr[edi],xmm4
	vmovdqa XMMWORD ptr[edi+16],xmm2
	add edi,32
	
	loop _AVX_3b_a
	
_AVX_3b_b:	
	mov ecx,w
	and ecx,1
	jz short _AVX_3b_c
	
	mov ebx,src1_u
	mov edx,src1_v
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	vmovd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	
	
	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc eax
	vpsraw xmm1,xmm1,3
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0
	
_AVX_3b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_3b endp


JPSDR_AutoYUY2_SSE2_4 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_4

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_2
	movdqa xmm5,XMMWORD ptr uw_3
	
	xor eax,eax
	or ecx,ecx
	jz short _SSE2_4_c

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
	
	movdqa XMMWORD ptr[edi],xmm0
	add edi,16
	
	loop _SSE2_4_a
	
_SSE2_4_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_4 endp


JPSDR_AutoYUY2_AVX_4 proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_4

	push esi
	push edi
	push ebx
	
	vpxor xmm7,xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_2
	vmovdqa xmm5,XMMWORD ptr uw_3
	
	xor eax,eax
	or ecx,ecx
	jz short _AVX_4_c

_AVX_4_a:
	mov ebx,src1_u
	mov edx,src1_v
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	vmovd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	
	
	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	inc eax
	vpsraw xmm1,xmm1,2
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0
	add edi,16
	
	loop _AVX_4_a
	
_AVX_4_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_4 endp


JPSDR_AutoYUY2_SSE2_4b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_SSE2_4b

	push esi
	push edi
	push ebx
	
	pxor xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	movdqa xmm6,XMMWORD ptr uw_2
	movdqa xmm5,XMMWORD ptr uw_3
	
	xor eax,eax
	shr ecx,1
	jz _SSE2_4b_b

_SSE2_4b_a:
	mov ebx,src1_u
	mov edx,src1_v
	movq xmm0,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm1,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	punpcklbw xmm0,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm1,xmm7				;0V0V0V0V0V0V0V0V
	mov ebx,src2_u
	mov edx,src2_v
	movq xmm2,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	movq xmm3,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	punpcklbw xmm2,xmm7				;0U0U0U0U0U0U0U0U
	punpcklbw xmm3,xmm7				;0V0V0V0V0V0V0V0V
	movdqa xmm4,XMMWORD ptr[esi+8*eax]		;YYYYYYYYYYYYYYYY
	
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
	add eax,2
	punpcklbw xmm2,xmm0				;VYUYVYUYVYUYVYUY
	punpckhbw xmm4,xmm0				;VYUYVYUYVYUYVYUY
	movdqa XMMWORD ptr[edi],xmm2
	movdqa XMMWORD ptr[edi+16],xmm4
	add edi,32
	
	loop _SSE2_4b_a
	
_SSE2_4b_b:	
	mov ecx,w
	and ecx,1
	jz short _SSE2_4b_c
	
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
	psraw xmm1,2
	packuswb xmm1,xmm7				;00000000VUVUVUVU
	punpcklbw xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	movdqa XMMWORD ptr[edi],xmm0
		
_SSE2_4b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_SSE2_4b endp


JPSDR_AutoYUY2_AVX_4b proc src_y:dword,src1_u:dword,src2_u:dword,src1_v:dword,src2_v:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_AVX_4b

	push esi
	push edi
	push ebx
	
	vpxor xmm7,xmm7,xmm7
	
	mov edi,dst
	mov esi,src_y
	mov ecx,w
	
	vmovdqa xmm6,XMMWORD ptr uw_2
	vmovdqa xmm5,XMMWORD ptr uw_3
	
	xor eax,eax
	shr ecx,1
	jz _AVX_4b_b

_AVX_4b_a:
	mov ebx,src1_u
	mov edx,src1_v
	vmovq xmm0,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm1,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	vpunpcklbw xmm0,xmm0,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm1,xmm1,xmm7			;0V0V0V0V0V0V0V0V
	mov ebx,src2_u
	mov edx,src2_v
	vmovq xmm2,qword ptr[ebx+4*eax]		;00000000UUUUUUUU
	vmovq xmm3,qword ptr[edx+4*eax]		;00000000VVVVVVVV
	vpunpcklbw xmm2,xmm2,xmm7			;0U0U0U0U0U0U0U0U
	vpunpcklbw xmm3,xmm3,xmm7			;0V0V0V0V0V0V0V0V
	vmovdqa xmm4,XMMWORD ptr[esi+8*eax]		;YYYYYYYYYYYYYYYY
	
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
	add eax,2
	vpunpcklbw xmm4,xmm2,xmm0			;VYUYVYUYVYUYVYUY
	vpunpckhbw xmm2,xmm2,xmm0			;VYUYVYUYVYUYVYUY
	vmovdqa XMMWORD ptr[edi],xmm4
	vmovdqa XMMWORD ptr[edi+16],xmm2
	add edi,32
	
	loop _AVX_4b_a
	
_AVX_4b_b:	
	mov ecx,w
	and ecx,1
	jz short _AVX_4b_c
	
	mov ebx,src1_u
	mov edx,src1_v
	vmovd xmm1,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm1,xmm1,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm1,xmm1,xmm7			;0V0U0V0U0V0U0V0U
	mov ebx,src2_u
	mov edx,src2_v
	vmovd xmm2,dword ptr[ebx+4*eax]		;000000000000UUUU
	vmovd xmm0,dword ptr[edx+4*eax]		;000000000000VVVV
	vpunpcklbw xmm2,xmm2,xmm0			;00000000VUVUVUVU
	vpunpcklbw xmm2,xmm2,xmm7			;0V0U0V0U0V0U0V0U	
	
	vpmullw xmm1,xmm1,xmm5
	vpaddsw xmm2,xmm2,xmm6
	vmovq xmm0,qword ptr[esi+8*eax]		;00000000YYYYYYYY
	vpaddsw xmm1,xmm1,xmm2
	vpsraw xmm1,xmm1,2
	vpackuswb xmm1,xmm1,xmm7			;00000000VUVUVUVU
	vpunpcklbw xmm0,xmm0,xmm1     		;VYUYVYUYVYUYVYUY
	
	vmovdqa XMMWORD ptr[edi],xmm0
		
_AVX_4b_c:		
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_AVX_4b endp


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


JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2 proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2

	push esi
	push edi
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short AVX_2_b

AVX_2_a:	
	vmovq xmm0,qword ptr[edx+8*eax]
	vmovq xmm1,qword ptr[esi+8*eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovq qword ptr[edi+8*eax],xmm2
	inc eax
	loop AVX_2_a

AVX_2_b:
	mov ecx,w
	and ecx,1
	jz short AVX_2_c
	
	vmovd xmm0,dword ptr[edx+8*eax]
	vmovd xmm1,dword ptr[esi+8*eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovd dword ptr[edi+8*eax],xmm2
	
AVX_2_c:	
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2 endp


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
	shr ecx,1
	jz short SSE2_2b_b

	mov ebx,16
SSE2_2b_a:	
	movdqa xmm0,XMMWORD ptr[edx+eax]
	movdqa xmm1,XMMWORD ptr[esi+eax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2

	movdqa XMMWORD ptr[edi+eax],xmm0
	add eax,ebx
	loop SSE2_2b_a
	
SSE2_2b_b:
	mov ecx,w
	and ecx,1
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
	
SSE2_2b_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_2b endp


JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b

	push esi
	push edi
	push ebx
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short AVX_2b_b

	mov ebx,16
AVX_2b_a:	
	vmovdqa xmm0,XMMWORD ptr[edx+eax]
	vmovdqa xmm1,XMMWORD ptr[esi+eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovdqa XMMWORD ptr[edi+eax],xmm2
	add eax,ebx
	loop AVX_2b_a
	
AVX_2b_b:
	mov ecx,w
	and ecx,1
	jz short AVX_2b_c
	
	vmovq xmm0,qword ptr[edx+eax]
	vmovq xmm1,qword ptr[esi+eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0

	vmovq qword ptr[edi+eax],xmm2
	
AVX_2b_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_2b endp


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


JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3 proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3
	
	push esi
	push edi
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short AVX_3_b
	
AVX_3_a:
	vmovq xmm0,qword ptr[esi+8*eax]
	vmovq xmm1,qword ptr[edx+8*eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovq qword ptr[edi+8*eax],xmm1
	inc eax
	loop AVX_3_a
	
AVX_3_b:
	mov ecx,w
	and ecx,1
	jz short AVX_3_c
	
	vmovd xmm0,dword ptr[esi+8*eax]
	vmovd xmm1,dword ptr[edx+8*eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovd dword ptr[edi+8*eax],xmm1	

AVX_3_c:
	pop edi
	pop esi
	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3 endp


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
	shr ecx,1
	jz short SSE2_3b_b
	
	mov ebx,16
SSE2_3b_a:
	movdqa xmm0,XMMWORD ptr[esi+eax]
	movdqa xmm1,XMMWORD ptr[edx+eax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm1,xmm0
	pavgb xmm1,xmm0
	pxor xmm1,xmm3
	pavgb xmm1,xmm2
	
	movdqa XMMWORD ptr[edi+eax],xmm1
	add eax,ebx
	loop SSE2_3b_a
	
SSE2_3b_b:	
	mov ecx,w
	and ecx,1
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
	
SSE2_3b_c:	
	pop ebx
	pop edi
	pop esi
	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_3b endp


JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short AVX_3b_b
	
	mov ebx,16
AVX_3b_a:
	vmovdqa xmm0,XMMWORD ptr[esi+eax]
	vmovdqa xmm1,XMMWORD ptr[edx+eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovdqa XMMWORD ptr[edi+eax],xmm1
	add eax,ebx
	loop AVX_3b_a
	
AVX_3b_b:	
	mov ecx,w
	and ecx,1
	jz short AVX_3b_c
	
	vmovq xmm0,qword ptr[esi+eax]
	vmovq xmm1,qword ptr[edx+eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm2
	vpavgb xmm1,xmm1,xmm2
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm1,xmm1,xmm0
	
	vmovq qword ptr[edi+eax],xmm1
	
AVX_3b_c:	
	pop ebx
	pop edi
	pop esi
	
	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_3b endp


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


JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4 proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4
	
	push esi
	push edi
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short AVX_4_b
	
AVX_4_a:
	vmovq xmm0,qword ptr[esi+8*eax]
	vmovq xmm1,qword ptr[edx+8*eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovq qword ptr[edi+8*eax],xmm2
	inc eax
	loop AVX_4_a
	
AVX_4_b:
	mov ecx,w
	and ecx,1
	jz short AVX_4_c

	vmovd xmm0,dword ptr[esi+8*eax]
	vmovd xmm1,dword ptr[edx+8*eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovd dword ptr[edi+8*eax],xmm2
	
AVX_4_c:
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4 endp


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
	shr ecx,1
	jz short SSE2_4b_b
	
	mov ebx,16
SSE2_4b_a:
	movdqa xmm0,XMMWORD ptr[esi+eax]
	movdqa xmm1,XMMWORD ptr[edx+eax]
	movdqa xmm2,xmm0
	pxor xmm0,xmm3
	pxor xmm1,xmm3
	pavgb xmm0,xmm1
	pxor xmm0,xmm3
	pavgb xmm0,xmm2
	
	movdqa XMMWORD ptr[edi+eax],xmm0
	add eax,ebx
	loop SSE2_4b_a
	
SSE2_4b_b:
	mov ecx,w
	and ecx,1
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
		
SSE2_4b_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_SSE2_4b endp


JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b proc src1:dword,src2:dword,dst:dword,w:dword

	public JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b
	
	push esi
	push edi
	push ebx
	
	vpcmpeqb xmm3,xmm3,xmm3
	
	mov edi,dst
	mov esi,src1
	mov edx,src2
	xor eax,eax
	
	mov ecx,w
	shr ecx,1
	jz short AVX_4b_b
	
	mov ebx,16
AVX_4b_a:
	vmovdqa xmm0,XMMWORD ptr[esi+eax]
	vmovdqa xmm1,XMMWORD ptr[edx+eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovdqa XMMWORD ptr[edi+eax],xmm2
	add eax,ebx
	loop AVX_4b_a
	
AVX_4b_b:
	mov ecx,w
	and ecx,1
	jz short AVX_4b_c

	vmovq xmm0,qword ptr[esi+eax]
	vmovq xmm1,qword ptr[edx+eax]
	vpxor xmm2,xmm0,xmm3
	vpxor xmm1,xmm1,xmm3
	vpavgb xmm2,xmm2,xmm1
	vpxor xmm2,xmm2,xmm3
	vpavgb xmm2,xmm2,xmm0
	
	vmovq qword ptr[edi+eax],xmm2
		
AVX_4b_c:	
	pop ebx
	pop edi
	pop esi

	ret

JPSDR_AutoYUY2_Convert420_to_Planar422_AVX_4b endp

end





