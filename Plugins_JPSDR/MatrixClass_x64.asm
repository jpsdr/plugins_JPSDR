.code


;CoeffProductF_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffProductF_SSE2 proc public frame

	.endprolog
	
	movss xmm0,dword ptr[rcx]
	pshufd xmm0,xmm0,0

	xor rcx,rcx
	mov rax,16
	mov ecx,r9d
	
CoeffProductF_SSE2_1:	
	movaps xmm1,XMMWORD ptr[rdx]
	add rdx,rax
	mulps xmm1,xmm0
	movaps XMMWORD ptr[r8],xmm1
	add r8,rax
	loop CoeffProductF_SSE2_1
	
	ret
	
CoeffProductF_SSE2 endp	


;CoeffProductF_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffProductF_AVX proc public frame

	.endprolog
	
	vbroadcastss ymm0,dword ptr[rcx]

	xor rcx,rcx
	mov rax,32
	mov ecx,r9d
	
CoeffProductF_AVX_1:	
	vmulps ymm1,ymm0,YMMWORD ptr[rdx]
	add rdx,rax
	vmovaps YMMWORD ptr[r8],ymm1
	add r8,rax
	loop CoeffProductF_AVX_1
	
	vzeroupper
	
	ret
	
CoeffProductF_AVX endp	


;CoeffProductD_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffProductD_SSE2 proc public frame
	
	.endprolog
	
	movsd xmm0,qword ptr[rcx]
	movlhps xmm0,xmm0
	
	xor rcx,rcx
	mov rax,16
	mov ecx,r9d
	
CoeffProductD_SSE2_1:	
	movapd xmm1,XMMWORD ptr[rdx]
	add rdx,rax
	mulpd xmm1,xmm0
	movapd XMMWORD ptr[r8],xmm1
	add r8,rax
	loop CoeffProductD_SSE2_1
	
	ret
	
CoeffProductD_SSE2 endp		


;CoeffProductD_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffProductD_AVX proc public frame

	.endprolog
	
	vbroadcastsd ymm0,qword ptr[rcx]
	
	xor rcx,rcx
	mov rax,32
	mov ecx,r9d
	
CoeffProductD_AVX_1:	
	vmulpd ymm1,ymm0,YMMWORD ptr[rdx]
	add rdx,rax
	vmovapd YMMWORD ptr[r8],ymm1
	add r8,rax
	loop CoeffProductD_AVX_1
	
	vzeroupper
		
	ret
	
CoeffProductD_AVX endp	


;CoeffAddProductF_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffAddProductF_SSE2 proc public frame
	
	.endprolog
	
	movss xmm0,dword ptr[rcx]
	pshufd xmm0,xmm0,0
	
	xor rcx,rcx
	mov rax,16
	mov ecx,r9d
	
CoeffAddProductF_SSE2_1:	
	movaps xmm1,XMMWORD ptr[rdx]
	movaps xmm2,XMMWORD ptr[r8]
	mulps xmm1,xmm0
	add rdx,rax
	addps xmm2,xmm1
	movaps XMMWORD ptr[r8],xmm2
	add r8,rax
	loop CoeffAddProductF_SSE2_1
	
	ret
	
CoeffAddProductF_SSE2 endp	


;CoeffAddProductF_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffAddProductF_AVX proc public frame

	.endprolog
	
	vbroadcastss ymm0,dword ptr[rcx]
	
	xor rcx,rcx
	mov rax,32
	mov ecx,r9d
	
CoeffAddProductF_AVX_1:	
	vmulps ymm1,ymm0,YMMWORD ptr[rdx]
	add rdx,rax
	vaddps ymm2,ymm1,YMMWORD ptr[r8]
	vmovaps YMMWORD ptr[r8],ymm2
	add r8,rax
	loop CoeffAddProductF_AVX_1
	
	vzeroupper
	
	ret
	
CoeffAddProductF_AVX endp


;CoeffAddProductD_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffAddProductD_SSE2 proc public frame
	
	.endprolog
	
	movsd xmm0,qword ptr[rcx]
	movlhps xmm0,xmm0
	
	xor rcx,rcx
	mov rax,16
	mov ecx,r9d

CoeffAddProductD_SSE2_1:	
	movapd xmm1,XMMWORD ptr[rdx]
	movapd xmm2,XMMWORD ptr[r8]
	mulpd xmm1,xmm0
	add rdx,rax
	addpd xmm2,xmm1
	movapd XMMWORD ptr[r8],xmm2
	add r8,rax
	loop CoeffAddProductD_SSE2_1
	
	ret
	
CoeffAddProductD_SSE2 endp	


;CoeffAddProductD_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word
; coeff_a = rcx
; coeff_b = rdx
; coeff_c = r8
; lgth = r9d

CoeffAddProductD_AVX proc public frame

	.endprolog
	
	vbroadcastsd ymm0,qword ptr[rcx]
	
	xor rcx,rcx
	mov rax,32
	mov ecx,r9d

CoeffAddProductD_AVX_1:	
	vmulpd ymm1,ymm0,YMMWORD ptr[rdx]
	add rdx,rax
	vaddpd ymm2,ymm1,YMMWORD ptr[r8]
	vmovapd YMMWORD ptr[r8],ymm2
	add r8,rax
	loop CoeffAddProductD_AVX_1
	
	vzeroupper
	
	ret
	
CoeffAddProductD_AVX endp	


;VectorProductF_SSE2 proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word
; coeff_a = rcx
; coeff_x = rdx
; result = r8
; lgth = r9d

VectorProductF_SSE2 proc public frame
		
	.endprolog
		
	mov r10,rcx
	pxor xmm0,xmm0
	xor rcx,rcx
	mov rax,16
	mov ecx,r9d
		
VectorProductF_SSE2_1:
	movaps xmm1,XMMWORD ptr[r10]
	add r10,rax
	mulps xmm1,XMMWORD ptr[rdx]
	add rdx,rax
	addps xmm0,xmm1
	loop VectorProductF_SSE2_1
		
	movhlps xmm1,xmm0
	addps xmm0,xmm1
	movaps xmm1,xmm0
	psrldq xmm1,4
	addps xmm0,xmm1
	movss dword ptr[r8],xmm0
		
	ret
		
VectorProductF_SSE2 endp		
		
		
;VectorProductF_AVX proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word
; coeff_a = rcx
; coeff_x = rdx
; result = r8
; lgth = r9d

VectorProductF_AVX proc public frame

	.endprolog
		
	mov r10,rcx
	vxorps ymm0,ymm0,ymm0
	xor rcx,rcx
	mov rax,32
	mov ecx,r9d
		
VectorProductF_AVX_1:
	vmovaps ymm1,YMMWORD ptr[r10]
	add r10,rax
	vmulps ymm1,ymm1,YMMWORD ptr[rdx]
	add rdx,rax
	vaddps ymm0,ymm0,ymm1
	loop VectorProductF_AVX_1
		
	vextractf128 xmm1,ymm0,1
	vaddps xmm0,xmm0,xmm1
		
	vmovhlps xmm1,xmm1,xmm0
	vaddps xmm0,xmm0,xmm1
	vpsrldq xmm1,xmm0,4
	vaddps xmm0,xmm0,xmm1
	vmovss dword ptr[r8],xmm0
		
	vzeroupper
	
	ret
		
VectorProductF_AVX endp

		
;VectorProductD_SSE2 proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word
; coeff_a = rcx
; coeff_x = rdx
; result = r8
; lgth = r9d

VectorProductD_SSE2 proc public frame
		
	.endprolog
		
	mov r10,rcx
	pxor xmm0,xmm0
	xor rcx,rcx
	mov rax,16
	mov ecx,r9d
		
VectorProductD_SSE2_1:
	movapd xmm1,XMMWORD ptr[r10]
	add r10,rax
	mulpd xmm1,XMMWORD ptr[rdx]
	add rdx,rax
	addpd xmm0,xmm1
	loop VectorProductD_SSE2_1
		
	movhlps xmm1,xmm0
	addpd xmm0,xmm1
	movsd qword ptr[r8],xmm0
		
	ret
		
VectorProductD_SSE2 endp		


;VectorProductD_AVX proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word
; coeff_a = rcx
; coeff_x = rdx
; result = r8
; lgth = r9d

VectorProductD_AVX proc public frame

	.endprolog
	
	mov r10,rcx
	vxorps ymm0,ymm0,ymm0
	xor rcx,rcx
	mov rax,32
	mov ecx,r9d
		
VectorProductD_AVX_1:
	vmovapd ymm1,YMMWORD ptr[r10]
	add r10,rax
	vmulpd ymm1,ymm1,YMMWORD ptr[rdx]
	add rdx,rax
	vaddpd ymm0,ymm0,ymm1
	loop VectorProductD_AVX_1
		
	vextractf128 xmm1,ymm0,1
	vaddpd xmm0,xmm0,xmm1
		
	vmovhlps xmm1,xmm1,xmm0
	vaddps xmm0,xmm0,xmm1
	vmovsd qword ptr[r8],xmm0
		
	vzeroupper
		
	ret
		
VectorProductD_AVX endp	

		
end