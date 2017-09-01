.xmm
.model flat,c

.code


CoeffProductF_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffProductF_SSE2
	
	push ebx
	
	mov edx,coeff_a
	
	movss xmm0,dword ptr[edx]
	pshufd xmm0,xmm0,0
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,16
	movzx ecx,lgth
	
CoeffProductF_SSE2_1:	
	movaps xmm1,XMMWORD ptr[ebx]
	add ebx,eax
	mulps xmm1,xmm0
	movaps XMMWORD ptr[edx],xmm1
	add edx,eax
	loop CoeffProductF_SSE2_1
	
	pop ebx
	
	ret
	
CoeffProductF_SSE2 endp	


CoeffProductF_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffProductF_AVX
	
	push ebx
	
	mov edx,coeff_a
	
	vbroadcastss ymm0,dword ptr[edx]
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,32
	movzx ecx,lgth
	
CoeffProductF_AVX_1:	
	vmulps ymm1,ymm0,YMMWORD ptr[ebx]
	add ebx,eax
	vmovaps YMMWORD ptr[edx],ymm1
	add edx,eax
	loop CoeffProductF_AVX_1
	
	vzeroupper
	
	pop ebx
	
	ret
	
CoeffProductF_AVX endp
	
	
CoeffProductD_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffProductD_SSE2
	
	push ebx
	
	mov edx,coeff_a
	
	movsd xmm0,qword ptr[edx]
	movlhps xmm0,xmm0
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,16
	movzx ecx,lgth
	
CoeffProductD_SSE2_1:	
	movapd xmm1,XMMWORD ptr[ebx]
	add ebx,eax
	mulpd xmm1,xmm0
	movapd XMMWORD ptr[edx],xmm1
	add edx,eax
	loop CoeffProductD_SSE2_1
	
	pop ebx
	
	ret
	
CoeffProductD_SSE2 endp		


CoeffProductD_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffProductD_AVX
	
	push ebx
	
	mov edx,coeff_a
	
	vbroadcastsd ymm0,qword ptr[edx]
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,32
	movzx ecx,lgth
	
CoeffProductD_AVX_1:	
	vmulpd ymm1,ymm0,YMMWORD ptr[ebx]
	add ebx,eax
	vmovapd YMMWORD ptr[edx],ymm1
	add edx,eax
	loop CoeffProductD_AVX_1
	
	vzeroupper
	
	pop ebx
	
	ret
	
CoeffProductD_AVX endp	


CoeffAddProductF_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffAddProductF_SSE2
	
	push ebx
	
	mov edx,coeff_a
	
	movss xmm0,dword ptr[edx]
	pshufd xmm0,xmm0,0
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,16
	movzx ecx,lgth
	
CoeffAddProductF_SSE2_1:	
	movaps xmm1,XMMWORD ptr[ebx]
	movaps xmm2,XMMWORD ptr[edx]
	mulps xmm1,xmm0
	add ebx,eax
	addps xmm2,xmm1
	movaps XMMWORD ptr[edx],xmm2
	add edx,eax
	loop CoeffAddProductF_SSE2_1
	
	pop ebx
	
	ret
	
CoeffAddProductF_SSE2 endp	


CoeffAddProductF_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffAddProductF_AVX
	
	push ebx
	
	mov edx,coeff_a
	
	vbroadcastss ymm0,dword ptr[edx]
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,32
	movzx ecx,lgth
	
CoeffAddProductF_AVX_1:	
	vmulps ymm1,ymm0,YMMWORD ptr[ebx]
	add ebx,eax
	vaddps ymm2,ymm1,YMMWORD ptr[edx]
	vmovaps YMMWORD ptr[edx],ymm2
	add edx,eax
	loop CoeffAddProductF_AVX_1
	
	vzeroupper
	
	pop ebx
	
	ret
	
CoeffAddProductF_AVX endp	


CoeffAddProductD_SSE2 proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffAddProductD_SSE2
	
	push ebx
	
	mov edx,coeff_a
	
	movsd xmm0,qword ptr[edx]
	movlhps xmm0,xmm0
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,16
	movzx ecx,lgth
	
CoeffAddProductD_SSE2_1:	
	movapd xmm1,XMMWORD ptr[ebx]
	movapd xmm2,XMMWORD ptr[edx]
	mulpd xmm1,xmm0
	add ebx,eax
	addpd xmm2,xmm1
	movapd XMMWORD ptr[edx],xmm2
	add edx,eax
	loop CoeffAddProductD_SSE2_1
	
	pop ebx
	
	ret
	
CoeffAddProductD_SSE2 endp	


CoeffAddProductD_AVX proc coeff_a:dword,coeff_b:dword,coeff_c:dword,lgth:word

    public CoeffAddProductD_AVX
	
	push ebx
	
	mov edx,coeff_a
	
	vbroadcastsd ymm0,qword ptr[edx]
	
	mov ebx,coeff_b
	mov edx,coeff_c
	mov eax,32
	movzx ecx,lgth
	
CoeffAddProductD_AVX_1:	
	vmulpd ymm1,ymm0,YMMWORD ptr[ebx]
	add ebx,eax
	vaddpd ymm2,ymm1,YMMWORD ptr[edx]
	vmovapd YMMWORD ptr[edx],ymm2
	add edx,eax
	loop CoeffAddProductD_AVX_1
	
	vzeroupper
	
	pop ebx
	
	ret
	
CoeffAddProductD_AVX endp	


VectorProductF_SSE2 proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word

    public VectorProductF_SSE2
	
	push ebx
		
	mov eax,16
	mov ebx,coeff_a
	mov edx,coeff_x
		
	pxor xmm0,xmm0
	movzx ecx,lgth
		
VectorProductF_SSE2_1:
	movaps xmm1,XMMWORD ptr[ebx]
	add ebx,eax
	mulps xmm1,XMMWORD ptr[edx]
	add edx,eax
	addps xmm0,xmm1
	loop VectorProductF_SSE2_1
		
	movhlps xmm1,xmm0
	addps xmm0,xmm1
	mov edx,result
	movaps xmm1,xmm0
	psrldq xmm1,4
	addps xmm0,xmm1
	movss dword ptr[edx],xmm0
		
	pop ebx
		
	ret
		
VectorProductF_SSE2 endp		
		
		
VectorProductF_AVX proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word

    public VectorProductF_AVX
	
	push ebx
		
	mov eax,32
	mov ebx,coeff_a
	mov edx,coeff_x
		
	vxorps ymm0,ymm0,ymm0
		
	movzx ecx,lgth
		
VectorProductF_AVX_1:
	vmovaps ymm1,YMMWORD ptr[ebx]
	add ebx,eax
	vmulps ymm1,ymm1,YMMWORD ptr[edx]
	add edx,eax
	vaddps ymm0,ymm0,ymm1
	loop VectorProductF_AVX_1
		
	vextractf128 xmm1,ymm0,1
	vaddps xmm0,xmm0,xmm1
		
	vmovhlps xmm1,xmm1,xmm0
	vaddps xmm0,xmm0,xmm1
	mov edx,result
	vpsrldq xmm1,xmm0,4
	vaddps xmm0,xmm0,xmm1
	vmovss dword ptr[edx],xmm0
		
	vzeroupper
		
	pop ebx
		
	ret
		
VectorProductF_AVX endp	

		
VectorProductD_SSE2 proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word

    public VectorProductD_SSE2
	
	push ebx
		
	mov eax,16
	mov ebx,coeff_a
	mov edx,coeff_x
		
	pxor xmm0,xmm0
	movzx ecx,lgth
		
VectorProductD_SSE2_1:
	movapd xmm1,XMMWORD ptr[ebx]
	add ebx,eax
	mulpd xmm1,XMMWORD ptr[edx]
	add edx,eax
	addpd xmm0,xmm1
	loop VectorProductD_SSE2_1
		
	movhlps xmm1,xmm0
	mov edx,result
	addpd xmm0,xmm1
	movsd qword ptr[edx],xmm0
		
	pop ebx
		
	ret
		
VectorProductD_SSE2 endp		


VectorProductD_AVX proc coeff_a:dword,coeff_x:dword,result:dword,lgth:word

    public VectorProductD_AVX
	
	push ebx
		
	mov eax,32
	mov ebx,coeff_a
	mov edx,coeff_x
		
	vxorps ymm0,ymm0,ymm0
		
	movzx ecx,lgth
		
VectorProductD_AVX_1:
	vmovapd ymm1,YMMWORD ptr[ebx]
	add ebx,eax
	vmulpd ymm1,ymm1,YMMWORD ptr[edx]
	add edx,eax
	vaddpd ymm0,ymm0,ymm1
	loop VectorProductD_AVX_1
		
	vextractf128 xmm1,ymm0,1
	vaddpd xmm0,xmm0,xmm1
		
	vmovhlps xmm1,xmm1,xmm0
	mov edx,result
	vaddps xmm0,xmm0,xmm1
	vmovsd qword ptr[edx],xmm0
		
	vzeroupper
		
	pop ebx
		
	ret
		
VectorProductD_AVX endp	

		
end