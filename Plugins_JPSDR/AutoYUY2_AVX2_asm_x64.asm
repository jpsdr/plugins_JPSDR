.data

.code


;JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_8_AVX2 proc src1:dword,src2:dword,dst:dword,w32:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w32 = r9d

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_8_AVX2 proc public frame

	.endprolog
		
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax	
	mov ecx,r9d	
	mov r11,32
	
Convert_Planar420_to_Planar422_x3x1_8_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[r10+rax]
	vmovdqa ymm1,YMMWORD ptr[rdx+rax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgb ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[r8+rax],ymm2
	add rax,r11
	loop Convert_Planar420_to_Planar422_x3x1_8_AVX2_1
	
	vzeroupper
	
	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_8_AVX2 endp


;JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_16_AVX2 proc src1:dword,src2:dword,dst:dword,w16:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w16 = r9d

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_16_AVX2 proc public frame

	.endprolog
		
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax	
	mov ecx,r9d	
	mov r11,32
	
Convert_Planar420_to_Planar422_x3x1_16_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[r10+rax]
	vmovdqa ymm1,YMMWORD ptr[rdx+rax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgw ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[r8+rax],ymm2
	add rax,r11
	loop Convert_Planar420_to_Planar422_x3x1_16_AVX2_1
	
	vzeroupper
	
	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x1_16_AVX2 endp


;JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_8_AVX2 proc src1:dword,src2:dword,dst:dword,w32:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w32 = r9d

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_8_AVX2 proc public frame

	.endprolog
		
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax	
	mov ecx,r9d	
	mov r11,32
	
Convert_Planar420_to_Planar422_x3x5_8_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[r10+rax]
	vmovdqa ymm1,YMMWORD ptr[rdx+rax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm2,ymm2,ymm1
	vpavgb ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgb ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[r8+rax],ymm2
	add rax,r11
	loop Convert_Planar420_to_Planar422_x3x5_8_AVX2_1
	
	vzeroupper
	
	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_8_AVX2 endp


;JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_16_AVX2 proc src1:dword,src2:dword,dst:dword,w16:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w16 = r9d

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_16_AVX2 proc public frame

	.endprolog
		
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax	
	mov ecx,r9d	
	mov r11,32
	
Convert_Planar420_to_Planar422_x3x5_16_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[r10+rax]
	vmovdqa ymm1,YMMWORD ptr[rdx+rax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm2,ymm2,ymm1
	vpavgw ymm2,ymm2,ymm1
	vpxor ymm2,ymm2,ymm3
	vpavgw ymm2,ymm2,ymm0
	
	vmovdqa YMMWORD ptr[r8+rax],ymm2
	add rax,r11
	loop Convert_Planar420_to_Planar422_x3x5_16_AVX2_1
	
	vzeroupper
	
	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x3x5_16_AVX2 endp


;JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_8_AVX2 proc src1:dword,src2:dword,dst:dword,w32:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w32 = r9d

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_8_AVX2 proc public frame

	.endprolog
		
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax	
	mov ecx,r9d	
	mov r11,32
	
Convert_Planar420_to_Planar422_x7x1_8_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[r10+rax]
	vmovdqa ymm1,YMMWORD ptr[rdx+rax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm1,ymm1,ymm2
	vpavgb ymm1,ymm1,ymm2
	vpxor ymm1,ymm1,ymm3
	vpavgb ymm1,ymm1,ymm0
	vmovdqa YMMWORD ptr[r8+rax],ymm1
	add rax,r11
	loop Convert_Planar420_to_Planar422_x7x1_8_AVX2_1
	
	vzeroupper
	
	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_8_AVX2 endp


;JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_16_AVX2 proc src1:dword,src2:dword,dst:dword,w16:dword
; src1 = rcx
; src2 = rdx
; dst = r8
; w16 = r9d

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_16_AVX2 proc public frame

	.endprolog
		
	vpcmpeqb ymm3,ymm3,ymm3
	
	mov r10,rcx				; r10=src1
	xor rcx,rcx
	xor rax,rax	
	mov ecx,r9d	
	mov r11,32
	
Convert_Planar420_to_Planar422_x7x1_16_AVX2_1:
	vmovdqa ymm0,YMMWORD ptr[r10+rax]
	vmovdqa ymm1,YMMWORD ptr[rdx+rax]
	vpxor ymm2,ymm0,ymm3
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm1,ymm1,ymm2
	vpavgw ymm1,ymm1,ymm2
	vpxor ymm1,ymm1,ymm3
	vpavgw ymm1,ymm1,ymm0
	vmovdqa YMMWORD ptr[r8+rax],ymm1
	add rax,r11
	loop Convert_Planar420_to_Planar422_x7x1_16_AVX2_1
	
	vzeroupper
	
	ret

JPSDR_AutoYUY2_Convert_Planar420_to_Planar422_x7x1_16_AVX2 endp


end
