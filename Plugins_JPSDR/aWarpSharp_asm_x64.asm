.code


;JPSDR_Warp0_8_SSE2 proc psrc:dword,pedg:dword,pdst:dword,src_pitch:dword,edg_pitchp:dword,edg_pitchn:dword,
;	y_limit_min:dword,y_limit_max:dword,x_limit_min:dword,x_limit_max:dword,i_:dword,depthH:dword,depthV:dword
; psrc = rcx
; pedg = rdx
; pdst = r8
; src_pitch = r9d
JPSDR_Warp0_8_SSE2 proc public frame

edg_pitchp equ dword ptr[rbp+48]
edg_pitchn equ dword ptr[rbp+56]
y_limit_min equ dword ptr[rbp+64]
y_limit_max equ dword ptr[rbp+72]
x_limit_min equ qword ptr[rbp+80]
x_limit_max equ qword ptr[rbp+88]
i_ equ dword ptr[rbp+96]
depthH equ dword ptr[rbp+104]
depthV equ dword ptr[rbp+112]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,144
	.allocstack 144
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	movdqu XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	movdqu XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	movdqu XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	movdqu XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	movdqu XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	movdqu XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	movdqu XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	.endprolog

	mov r10,8
	
	mov rsi,x_limit_min
	movdqu xmm12,XMMWORD ptr[rsi]
	mov rsi,x_limit_max
	movdqu xmm13,XMMWORD ptr[rsi]
	mov	rsi,rcx
	mov	rcx,rdx
	mov	r11,r8
	movsxd rdx,r9d
	movsxd rbx,edg_pitchp
	movsxd rdi,i_ 					;signed!
    sub	r11,r10
    add	rdx,rsi
    add	rbx,rcx
	
    movd xmm10,y_limit_min
    movd xmm11,y_limit_max
    movd xmm9,depthH
    movd xmm8,r9d
    pcmpeqw xmm7,xmm7
    psrlw xmm7,15
    punpcklwd xmm8,xmm7
	movd xmm14,depthV
    pshufd xmm10,xmm10,0
    pshufd xmm11,xmm11,0
    pshufd xmm9,xmm9,0
    pshufd xmm8,xmm8,0
	pshufd xmm14,xmm14,0
    packssdw xmm10,xmm10
    packssdw xmm11,xmm11
    packssdw xmm9,xmm9
	packssdw xmm14,xmm14
	movsxd r9,edg_pitchn
    pcmpeqw xmm5,xmm5
    psllw xmm5,15
	
    add	r9,rcx
	
	movdqa xmm7,XMMWORD ptr[rdi+rcx]
	
    movdqa xmm1,xmm7
    pslldq xmm7,7
    punpcklqdq xmm7,xmm1
    psrldq xmm1,1
    psrldq xmm7,7
	psrlw xmm5,8
	
JPSDR_Warp0_8_SSE2_1:
	movq xmm4,qword ptr[rdi+rbx]
	movq xmm2,qword ptr[rdi+r9]
    pxor xmm3,xmm3
    punpcklbw xmm7,xmm3
    punpcklbw xmm1,xmm3
    punpcklbw xmm4,xmm3
    punpcklbw xmm2,xmm3
    psubw xmm7,xmm1
    psubw xmm4,xmm2
    psllw xmm7,7
    psllw xmm4,7
    pmulhw xmm7,xmm9 ;depthH
    pmulhw xmm4,xmm14 ;depthV
    pmaxsw xmm4,xmm10	;y_limit_min
    pminsw xmm4,xmm11	;y_limit_max
	
    pcmpeqw xmm3,xmm3
    psrlw xmm3,9
    movdqa xmm1,xmm7
    movdqa xmm2,xmm4
	
    pand xmm7,xmm3 ;007F
    pand xmm4,xmm3 ;007F
    psraw xmm1,7
    psraw xmm2,7

    movd xmm3,edi
    pshufd xmm3,xmm3,0
	
    packssdw xmm3,xmm3
    paddsw xmm1,xmm3
	
    movdqa xmm3,xmm13	;x_limit_max
    movdqa xmm0,xmm12	;x_limit_min
	
    pcmpgtw xmm3,xmm1
    pcmpgtw xmm0,xmm1
	
    pminsw xmm1,xmm13	;x_limit_max
    pmaxsw xmm1,xmm12	;x_limit_min
		  
    pand xmm7,xmm3
    pandn xmm0,xmm7

    movdqa xmm7,xmm2
    punpcklwd xmm2,xmm1
    punpckhwd xmm7,xmm1
    pmaddwd xmm2,xmm8 ;1 src_pitch
    pmaddwd xmm7,xmm8 ;1 src_pitch

    movd eax,xmm2
	movsxd rax,eax
    psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi],0
    pinsrw xmm1,word ptr[rax+rdx],0
    movd eax,xmm2
	movsxd rax,eax
    psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+1],1
    pinsrw xmm1,word ptr[rax+rdx+1],1
    movd eax,xmm2
	movsxd rax,eax
    psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+2],2
    pinsrw xmm1,word ptr[rax+rdx+2],2
    movd eax,xmm2
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+3],3
    pinsrw xmm1,word ptr[rax+rdx+3],3
    movd eax,xmm7
	movsxd rax,eax
    psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+4],4
    pinsrw xmm1,word ptr[rax+rdx+4],4
    movd eax,xmm7
	movsxd rax,eax
    psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+5],5
    pinsrw xmm1,word ptr[rax+rdx+5],5
    movd eax,xmm7
	movsxd rax,eax
    psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+6],6
    pinsrw xmm1,word ptr[rax+rdx+6],6
    movd eax,xmm7
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+7],7
    pinsrw xmm1,word ptr[rax+rdx+7],7
	
    pcmpeqw xmm6,xmm6
    movdqa xmm2,xmm3
    psrlw xmm6,8
    movdqa xmm7,xmm1
    pand xmm3,xmm6 ;00FF
    pand xmm1,xmm6 ;00FF
    movdqa xmm6,xmm5
    psubw xmm6,xmm0
    pmullw xmm3,xmm6
    pmullw xmm1,xmm6
    movdqa xmm6,xmm5 ;0080
    psrlw xmm5,1
    psrlw xmm2,8
    psrlw xmm7,8
    pmullw xmm2,xmm0
    pmullw xmm7,xmm0
    paddw xmm3,xmm2
    paddw xmm1,xmm7
    paddw xmm3,xmm5 ;0040
    paddw xmm1,xmm5 ;0040
    psraw xmm3,7
    psraw xmm1,7

    psubw xmm6,xmm4
    movdqu xmm7,XMMWORD ptr[rdi+rcx+7]
    pmullw xmm1,xmm4
    pmullw xmm3,xmm6
    paddw xmm3,xmm1
    movdqa xmm1,xmm7
	
    paddw xmm3,xmm5  ;0040
    psrldq xmm1,2
    psraw xmm3,7
    paddw xmm5,xmm5
    packuswb xmm3,xmm3
    add rdi,r10
    jg short JPSDR_Warp0_8_SSE2_2
    movq qword ptr[rdi+r11],xmm3
    jnz JPSDR_Warp0_8_SSE2_1
    jmp short JPSDR_Warp0_8_SSE2_Fin
	
JPSDR_Warp0_8_SSE2_2:
	movd dword ptr[rdi+r11],xmm3
	
JPSDR_Warp0_8_SSE2_Fin:
	movdqu xmm14,XMMWORD ptr[rsp+128]
	movdqu xmm13,XMMWORD ptr[rsp+112]
	movdqu xmm12,XMMWORD ptr[rsp+96]
	movdqu xmm11,XMMWORD ptr[rsp+80]
	movdqu xmm10,XMMWORD ptr[rsp+64]
	movdqu xmm9,XMMWORD ptr[rsp+48]
	movdqu xmm8,XMMWORD ptr[rsp+32]
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,144

	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_Warp0_8_SSE2 endp


;JPSDR_Warp0_8_SSE3 proc psrc:dword,pedg:dword,pdst:dword,src_pitch:dword,edg_pitchp:dword,edg_pitchn:dword,
;	y_limit_min:dword,y_limit_max:dword,x_limit_min:dword,x_limit_max:dword,i_:dword,depthH:dword,depthV:dword
; psrc = rcx
; pedg = rdx
; pdst = r8
; src_pitch = r9d
JPSDR_Warp0_8_SSE3 proc public frame

edg_pitchp equ dword ptr[rbp+48]
edg_pitchn equ dword ptr[rbp+56]
y_limit_min equ dword ptr[rbp+64]
y_limit_max equ dword ptr[rbp+72]
x_limit_min equ qword ptr[rbp+80]
x_limit_max equ qword ptr[rbp+88]
i_ equ dword ptr[rbp+96]
depthH equ dword ptr[rbp+104]
depthV equ dword ptr[rbp+112]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,144
	.allocstack 144
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	movdqu XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	movdqu XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	movdqu XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	movdqu XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	movdqu XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	movdqu XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	movdqu XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	.endprolog

	mov r10,8
	
	mov rsi,x_limit_min
	movdqu xmm12,XMMWORD ptr[rsi]
	mov rsi,x_limit_max
	movdqu xmm13,XMMWORD ptr[rsi]
	mov	rsi,rcx
	mov	rcx,rdx
	mov	r11,r8
	movsxd rdx,r9d
	movsxd rbx,edg_pitchp
	movsxd rdi,i_ 					;signed!
    sub	r11,r10
    add	rdx,rsi
    add	rbx,rcx
	
    movd xmm10,y_limit_min
    movd xmm11,y_limit_max
    movd xmm9,depthH
    movd xmm8,r9d
    pcmpeqw xmm7,xmm7
    psrlw xmm7,15
    punpcklwd xmm8,xmm7
	movd xmm14,depthV
    pshufd xmm10,xmm10,0
    pshufd xmm11,xmm11,0
    pshufd xmm9,xmm9,0
    pshufd xmm8,xmm8,0
	pshufd xmm14,xmm14,0
    packssdw xmm10,xmm10
    packssdw xmm11,xmm11
    packssdw xmm9,xmm9
	packssdw xmm14,xmm14
	movsxd r9,edg_pitchn
    pcmpeqw xmm5,xmm5
    psllw xmm5,15
	
    add	r9,rcx
	
	movdqa xmm7,XMMWORD ptr[rdi+rcx]

    movdqa xmm1,xmm7
    pslldq xmm7,7
    punpcklqdq xmm7,xmm1
    psrldq xmm1,1
    psrldq xmm7,7
	
JPSDR_Warp0_8_SSE3_1:
    movq xmm4,qword ptr[rdi+rbx]
    movq xmm2,qword ptr[rdi+r9]
    pxor xmm0,xmm0
    punpcklbw xmm7,xmm0
    punpcklbw xmm1,xmm0
    punpcklbw xmm4,xmm0
    punpcklbw xmm2,xmm0
    psubw xmm7,xmm1
    psubw xmm4,xmm2
    psllw xmm7,7
    psllw xmm4,7
    pmulhw xmm7,xmm9 ;depthH
    pmulhw xmm4,xmm14 ;depthV
    movd xmm6,edi ;preload

    pmaxsw xmm4,xmm10	;y_limit_min
    pminsw xmm4,xmm11	;y_limit_max
	
	pshufd xmm6,xmm6,0
	
	pcmpeqw xmm0,xmm0
    psrlw xmm0,9
    movdqa xmm2,xmm4
    movdqa xmm1,xmm7
	
    pand xmm4,xmm0 ;007F
    pand xmm7,xmm0 ;007F

    psraw xmm1,7
    packssdw xmm6,xmm6
    paddsw xmm1,xmm6
	
    movdqa xmm0,xmm13	;x_limit_max
    movdqa xmm3,xmm12	;x_limit_min

    pcmpgtw xmm0,xmm1
    pcmpgtw xmm3,xmm1
	
    pminsw xmm1,xmm13	;x_limit_max
    pmaxsw xmm1,xmm12	;x_limit_min

    pand xmm7,xmm0
    pandn xmm3,xmm7

    psraw xmm2,7
    movdqa xmm7,xmm2
    punpcklwd xmm2,xmm1
    punpckhwd xmm7,xmm1
    pmaddwd xmm2,xmm8 ;1 src_pitch
    pmaddwd xmm7,xmm8 ;1 src_pitch

    psignw xmm3,xmm5 ;8000
    psignw xmm4,xmm5 ;8000
    packsswb xmm5,xmm5
    packsswb xmm3,xmm3
    packsswb xmm4,xmm4

    movdqa xmm0,xmm5 ;80
    movdqa xmm1,xmm5 ;80
    psubb xmm0,xmm3
    psubb xmm1,xmm4
    psrlw xmm5,9
    punpcklbw xmm0,xmm3
    punpcklbw xmm1,xmm4

    movd eax,xmm2
	movsxd rax,eax
	psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi],0
    pinsrw xmm4,word ptr[rax+rdx],0
    movd eax,xmm2
	movsxd rax,eax
	psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+1],1
    pinsrw xmm4,word ptr[rax+rdx+1],1
    movd eax,xmm2
	movsxd rax,eax
	psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+2],2
    pinsrw xmm4,word ptr[rax+rdx+2],2
    movd eax,xmm2
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+3],3
    pinsrw xmm4,word ptr[rax+rdx+3],3
    movd eax,xmm7
	movsxd rax,eax
	psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+4],4
    pinsrw xmm4,word ptr[rax+rdx+4],4
    movd eax,xmm7
	movsxd rax,eax
	psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+5],5
    pinsrw xmm4,word ptr[rax+rdx+5],5
    movd eax,xmm7
	movsxd rax,eax
	psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+6],6
    pinsrw xmm4,word ptr[rax+rdx+6],6
    movd eax,xmm7
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+7],7
    pinsrw xmm4,word ptr[rax+rdx+7],7

    pcmpeqw xmm2,xmm2
    movdqu xmm7,XMMWORD ptr[rdi+rcx+7]
    pmaddubsw xmm3,xmm0
    pmaddubsw xmm4,xmm0
	
    psignw xmm3,xmm2
    psignw xmm4,xmm2
    paddw xmm3,xmm5 ; 0040
    paddw xmm4,xmm5 ; 0040
    psraw xmm3,7
    psraw xmm4,7
    packuswb xmm3,xmm3
    packuswb xmm4,xmm4
    punpcklbw xmm3,xmm4
    pmaddubsw xmm3,xmm1
    palignr xmm1,xmm7,2
    psignw xmm3,xmm2

    paddw xmm3,xmm5 ; 0040
    psraw xmm3,7
    psllw xmm5,9
    packuswb xmm3,xmm3
		  
    add rdi,r10
    jg short JPSDR_Warp0_8_SSE3_2
    movq qword ptr[rdi+r11],xmm3
    jnz JPSDR_Warp0_8_SSE3_1
    jmp short JPSDR_Warp0_8_SSE3_Fin
	
JPSDR_Warp0_8_SSE3_2:
	movd dword ptr[rdi+r11],xmm3
	
JPSDR_Warp0_8_SSE3_Fin:
	movdqu xmm14,XMMWORD ptr[rsp+128]
	movdqu xmm13,XMMWORD ptr[rsp+112]
	movdqu xmm12,XMMWORD ptr[rsp+96]
	movdqu xmm11,XMMWORD ptr[rsp+80]
	movdqu xmm10,XMMWORD ptr[rsp+64]
	movdqu xmm9,XMMWORD ptr[rsp+48]
	movdqu xmm8,XMMWORD ptr[rsp+32]
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,144

	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_Warp0_8_SSE3 endp	


;JPSDR_Warp2_8_SSE2 proc psrc:dword,pedg:dword,pdst:dword,src_pitch:dword,edg_pitchp:dword,edg_pitchn:dword,
;	y_limit_min:dword,y_limit_max:dword,x_limit_min:dword,x_limit_max:dword,i_:dword,depthH:dword,depthV:dword
; psrc = rcx
; pedg = rdx
; pdst = r8
; src_pitch = r9d
JPSDR_Warp2_8_SSE2 proc public frame

edg_pitchp equ dword ptr[rbp+48]
edg_pitchn equ dword ptr[rbp+56]
y_limit_min equ dword ptr[rbp+64]
y_limit_max equ dword ptr[rbp+72]
x_limit_min equ qword ptr[rbp+80]
x_limit_max equ qword ptr[rbp+88]
i_ equ dword ptr[rbp+96]
depthH equ dword ptr[rbp+104]
depthV equ dword ptr[rbp+112]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,144
	.allocstack 144
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	movdqu XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	movdqu XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	movdqu XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	movdqu XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	movdqu XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	movdqu XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	movdqu XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	.endprolog

	mov r10,8
	
	mov rsi,x_limit_min
	movdqu xmm12,XMMWORD ptr[rsi]
	mov rsi,x_limit_max
	movdqu xmm13,XMMWORD ptr[rsi]
	mov	rsi,rcx
	mov	rcx,rdx
	mov	r11,r8
	movsxd rdx,r9d
	movsxd rbx,edg_pitchp
	movsxd rdi,i_ 					;signed!
    sub	r11,r10
    add	rdx,rsi
    add	rbx,rcx
	
    movd xmm10,y_limit_min
    movd xmm11,y_limit_max
    movd xmm9,depthH
    movd xmm8,r9d
    pcmpeqw xmm7,xmm7
    psrlw xmm7,15
    punpcklwd xmm8,xmm7
	movd xmm14,depthV
    pshufd xmm10,xmm10,0
    pshufd xmm11,xmm11,0
    pshufd xmm9,xmm9,0
    pshufd xmm8,xmm8,0
	pshufd xmm14,xmm14,0
    packssdw xmm10,xmm10
    packssdw xmm11,xmm11
    packssdw xmm9,xmm9
	packssdw xmm14,xmm14
	movsxd r9,edg_pitchn
    pcmpeqw xmm5,xmm5
    psllw xmm5,15
	
    add	r9,rcx
	
	movdqa xmm7,XMMWORD ptr[rdi+rcx]
	
    movdqa xmm1,xmm7
    pslldq xmm7,7
    punpcklqdq xmm7,xmm1
    psrldq xmm1,1
    psrldq xmm7,7
	psrlw xmm5,8

JPSDR_Warp2_8_SSE2_1:
	movq xmm4,qword ptr[rdi+rbx]
	movq xmm2,qword ptr[rdi+r9]
    pxor xmm3,xmm3
    punpcklbw xmm7,xmm3
    punpcklbw xmm1,xmm3
    punpcklbw xmm4,xmm3
    punpcklbw xmm2,xmm3
    psubw xmm7,xmm1
    psubw xmm4,xmm2
    psllw xmm7,7
    psllw xmm4,7
    pmulhw xmm7,xmm9 ;depthH
    pmulhw xmm4,xmm14 ;depthV
    pmaxsw xmm4,xmm10	;y_limit_min
    pminsw xmm4,xmm11	;y_limit_max
	
    pcmpeqw xmm3,xmm3
    psrlw xmm3,9
    movdqa xmm1,xmm7
    movdqa xmm2,xmm4
	
    psllw xmm4,2 ;SMAGL=2
    psllw xmm7,2 ;SMAGL=2
	
    pand xmm7,xmm3 ;007F
    pand xmm4,xmm3 ;007F
    psraw xmm1,5	;7-SMAGL
    psraw xmm2,5	;7-SMAGL

    movd xmm3,edi
    pshufd xmm3,xmm3,0
	
	pslld xmm3,2	;SMAGL=2
	
    packssdw xmm3,xmm3
    paddsw xmm1,xmm3
	
    movdqa xmm3,xmm13	;x_limit_max
    movdqa xmm0,xmm12	;x_limit_min
	
    pcmpgtw xmm3,xmm1
    pcmpgtw xmm0,xmm1
	
    pminsw xmm1,xmm13	;x_limit_max
    pmaxsw xmm1,xmm12	;x_limit_min
		  
    pand xmm7,xmm3
    pandn xmm0,xmm7

    movdqa xmm7,xmm2
    punpcklwd xmm2,xmm1
    punpckhwd xmm7,xmm1
    pmaddwd xmm2,xmm8 ;1 src_pitch
    pmaddwd xmm7,xmm8 ;1 src_pitch

    movd eax,xmm2
	movsxd rax,eax
    psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi],0
    pinsrw xmm1,word ptr[rax+rdx],0
    movd eax,xmm2
	movsxd rax,eax
    psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+1*4],1
    pinsrw xmm1,word ptr[rax+rdx+1*4],1
    movd eax,xmm2
	movsxd rax,eax
    psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+2*4],2
    pinsrw xmm1,word ptr[rax+rdx+2*4],2
    movd eax,xmm2
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+3*4],3
    pinsrw xmm1,word ptr[rax+rdx+3*4],3
    movd eax,xmm7
	movsxd rax,eax
    psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+4*4],4
    pinsrw xmm1,word ptr[rax+rdx+4*4],4
    movd eax,xmm7
	movsxd rax,eax
    psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+5*4],5
    pinsrw xmm1,word ptr[rax+rdx+5*4],5
    movd eax,xmm7
	movsxd rax,eax
    psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+6*4],6
    pinsrw xmm1,word ptr[rax+rdx+6*4],6
    movd eax,xmm7
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+7*4],7
    pinsrw xmm1,word ptr[rax+rdx+7*4],7
	
    pcmpeqw xmm6,xmm6
    movdqa xmm2,xmm3
    psrlw xmm6,8
    movdqa xmm7,xmm1
    pand xmm3,xmm6 ;00FF
    pand xmm1,xmm6 ;00FF
    movdqa xmm6,xmm5
    psubw xmm6,xmm0
    pmullw xmm3,xmm6
    pmullw xmm1,xmm6
    movdqa xmm6,xmm5 ;0080
    psrlw xmm5,1
    psrlw xmm2,8
    psrlw xmm7,8
    pmullw xmm2,xmm0
    pmullw xmm7,xmm0
    paddw xmm3,xmm2
    paddw xmm1,xmm7
    paddw xmm3,xmm5 ;0040
    paddw xmm1,xmm5 ;0040
    psraw xmm3,7
    psraw xmm1,7

    psubw xmm6,xmm4
    movdqu xmm7,XMMWORD ptr[rdi+rcx+7]
    pmullw xmm1,xmm4
    pmullw xmm3,xmm6
    paddw xmm3,xmm1
    movdqa xmm1,xmm7
	
    paddw xmm3,xmm5  ;0040
    psrldq xmm1,2
    psraw xmm3,7
    paddw xmm5,xmm5
    packuswb xmm3,xmm3
    add rdi,r10
    jg short JPSDR_Warp2_8_SSE2_2
    movq qword ptr[rdi+r11],xmm3
    jnz JPSDR_Warp2_8_SSE2_1
    jmp short JPSDR_Warp2_8_SSE2_Fin
	
JPSDR_Warp2_8_SSE2_2:
	movd dword ptr[rdi+r11],xmm3
	
JPSDR_Warp2_8_SSE2_Fin:
	movdqu xmm14,XMMWORD ptr[rsp+128]
	movdqu xmm13,XMMWORD ptr[rsp+112]
	movdqu xmm12,XMMWORD ptr[rsp+96]
	movdqu xmm11,XMMWORD ptr[rsp+80]
	movdqu xmm10,XMMWORD ptr[rsp+64]
	movdqu xmm9,XMMWORD ptr[rsp+48]
	movdqu xmm8,XMMWORD ptr[rsp+32]
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,144

	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_Warp2_8_SSE2 endp	


;JPSDR_Warp2_8_SSE3 proc psrc:dword,pedg:dword,pdst:dword,src_pitch:dword,edg_pitchp:dword,edg_pitchn:dword,
;	y_limit_min:dword,y_limit_max:dword,x_limit_min:dword,x_limit_max:dword,i_:dword,depthH:dword,depthV:dword
; psrc = rcx
; pedg = rdx
; pdst = r8
; src_pitch = r9d
JPSDR_Warp2_8_SSE3 proc public frame

edg_pitchp equ dword ptr[rbp+48]
edg_pitchn equ dword ptr[rbp+56]
y_limit_min equ dword ptr[rbp+64]
y_limit_max equ dword ptr[rbp+72]
x_limit_min equ qword ptr[rbp+80]
x_limit_max equ qword ptr[rbp+88]
i_ equ dword ptr[rbp+96]
depthH equ dword ptr[rbp+104]
depthV equ dword ptr[rbp+112]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	push rbx
	.pushreg rbx
	sub rsp,144
	.allocstack 144
	movdqu XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqu XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	movdqu XMMWORD ptr[rsp+32],xmm8
	.savexmm128 xmm8,32
	movdqu XMMWORD ptr[rsp+48],xmm9
	.savexmm128 xmm9,48
	movdqu XMMWORD ptr[rsp+64],xmm10
	.savexmm128 xmm10,64
	movdqu XMMWORD ptr[rsp+80],xmm11
	.savexmm128 xmm11,80
	movdqu XMMWORD ptr[rsp+96],xmm12
	.savexmm128 xmm12,96
	movdqu XMMWORD ptr[rsp+112],xmm13
	.savexmm128 xmm13,112
	movdqu XMMWORD ptr[rsp+128],xmm14
	.savexmm128 xmm14,128
	.endprolog

	mov r10,8
	
	mov rsi,x_limit_min
	movdqu xmm12,XMMWORD ptr[rsi]
	mov rsi,x_limit_max
	movdqu xmm13,XMMWORD ptr[rsi]
	mov	rsi,rcx
	mov	rcx,rdx
	mov	r11,r8
	movsxd rdx,r9d
	movsxd rbx,edg_pitchp
	movsxd rdi,i_ 					;signed!
    sub	r11,r10
    add	rdx,rsi
    add	rbx,rcx
	
    movd xmm10,y_limit_min
    movd xmm11,y_limit_max
    movd xmm9,depthH
    movd xmm8,r9d
    pcmpeqw xmm7,xmm7
    psrlw xmm7,15
    punpcklwd xmm8,xmm7
	movd xmm14,depthV
    pshufd xmm10,xmm10,0
    pshufd xmm11,xmm11,0
    pshufd xmm9,xmm9,0
    pshufd xmm8,xmm8,0
	pshufd xmm14,xmm14,0
    packssdw xmm10,xmm10
    packssdw xmm11,xmm11
    packssdw xmm9,xmm9
	packssdw xmm14,xmm14
	movsxd r9,edg_pitchn
    pcmpeqw xmm5,xmm5
    psllw xmm5,15
	
    add	r9,rcx
	
	movdqa xmm7,XMMWORD ptr[rdi+rcx]

    movdqa xmm1,xmm7
    pslldq xmm7,7
    punpcklqdq xmm7,xmm1
    psrldq xmm1,1
    psrldq xmm7,7
	
JPSDR_Warp2_8_SSE3_1:
    movq xmm4,qword ptr[rdi+rbx]
    movq xmm2,qword ptr[rdi+r9]
    pxor xmm0,xmm0
    punpcklbw xmm7,xmm0
    punpcklbw xmm1,xmm0
    punpcklbw xmm4,xmm0
    punpcklbw xmm2,xmm0
    psubw xmm7,xmm1
    psubw xmm4,xmm2
    psllw xmm7,7
    psllw xmm4,7
    pmulhw xmm7,xmm9 ;depthH
    pmulhw xmm4,xmm14 ;depthV
    movd xmm6,edi ;preload

    pmaxsw xmm4,xmm10	;y_limit_min
    pminsw xmm4,xmm11	;y_limit_max
	
	pshufd xmm6,xmm6,0
	
	pslld xmm6,2	;SMAGL=2
	
	pcmpeqw xmm0,xmm0
    psrlw xmm0,9
    movdqa xmm2,xmm4
    movdqa xmm1,xmm7
	
    psllw xmm4,2	;SMAGL=2
    psllw xmm7,2	;SMAGL=2
	
    pand xmm4, xmm0 ;007F
    pand xmm7, xmm0 ;007F

    psraw xmm1,5	;7-SMAGL
    packssdw xmm6,xmm6
    paddsw xmm1,xmm6
	
    movdqa xmm0,xmm13	;x_limit_max
    movdqa xmm3,xmm12	;x_limit_min

    pcmpgtw xmm0,xmm1
    pcmpgtw xmm3,xmm1
	
    pminsw xmm1,xmm13	;x_limit_max
    pmaxsw xmm1,xmm12	;x_limit_min

    pand xmm7,xmm0
    pandn xmm3,xmm7

    psraw xmm2,5	;7-SMAGL
    movdqa xmm7,xmm2
    punpcklwd xmm2,xmm1
    punpckhwd xmm7,xmm1
    pmaddwd xmm2,xmm8 ;1 src_pitch
    pmaddwd xmm7,xmm8 ;1 src_pitch

    psignw xmm3,xmm5 ;8000
    psignw xmm4,xmm5 ;8000
    packsswb xmm5,xmm5
    packsswb xmm3,xmm3
    packsswb xmm4,xmm4

    movdqa xmm0,xmm5 ;80
    movdqa xmm1,xmm5 ;80
    psubb xmm0,xmm3
    psubb xmm1,xmm4
    psrlw xmm5,9
    punpcklbw xmm0,xmm3
    punpcklbw xmm1,xmm4

    movd eax,xmm2
	movsxd rax,eax
	psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi],0
    pinsrw xmm4,word ptr[rax+rdx],0
    movd eax,xmm2
	movsxd rax,eax
	psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+1*4],1
    pinsrw xmm4,word ptr[rax+rdx+1*4],1
    movd eax,xmm2
	movsxd rax,eax
	psrldq xmm2,4
    pinsrw xmm3,word ptr[rax+rsi+2*4],2
    pinsrw xmm4,word ptr[rax+rdx+2*4],2
    movd eax,xmm2
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+3*4],3
    pinsrw xmm4,word ptr[rax+rdx+3*4],3
    movd eax,xmm7
	movsxd rax,eax
	psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+4*4],4
    pinsrw xmm4,word ptr[rax+rdx+4*4],4
    movd eax,xmm7
	movsxd rax,eax
	psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+5*4],5
    pinsrw xmm4,word ptr[rax+rdx+5*4],5
    movd eax,xmm7
	movsxd rax,eax
	psrldq xmm7,4
    pinsrw xmm3,word ptr[rax+rsi+6*4],6
    pinsrw xmm4,word ptr[rax+rdx+6*4],6
    movd eax,xmm7
	movsxd rax,eax
    pinsrw xmm3,word ptr[rax+rsi+7*4],7
    pinsrw xmm4,word ptr[rax+rdx+7*4],7

    pcmpeqw xmm2,xmm2
    movdqu xmm7,XMMWORD ptr[rdi+rcx+7]
    pmaddubsw xmm3,xmm0
    pmaddubsw xmm4,xmm0
	
    psignw xmm3,xmm2
    psignw xmm4,xmm2
    paddw xmm3,xmm5 ; 0040
    paddw xmm4,xmm5 ; 0040
    psraw xmm3,7
    psraw xmm4,7
    packuswb xmm3,xmm3
    packuswb xmm4,xmm4
    punpcklbw xmm3,xmm4
    pmaddubsw xmm3,xmm1
    palignr xmm1,xmm7,2
    psignw xmm3,xmm2

    paddw xmm3,xmm5 ; 0040
    psraw xmm3,7
    psllw xmm5,9
    packuswb xmm3,xmm3
		  
    add rdi,r10
    jg short JPSDR_Warp2_8_SSE3_2
    movq qword ptr[rdi+r11],xmm3
    jnz JPSDR_Warp2_8_SSE3_1
    jmp short JPSDR_Warp2_8_SSE3_Fin
	
JPSDR_Warp2_8_SSE3_2:
	movd dword ptr[rdi+r11],xmm3
	
JPSDR_Warp2_8_SSE3_Fin:
	movdqu xmm14,XMMWORD ptr[rsp+128]
	movdqu xmm13,XMMWORD ptr[rsp+112]
	movdqu xmm12,XMMWORD ptr[rsp+96]
	movdqu xmm11,XMMWORD ptr[rsp+80]
	movdqu xmm10,XMMWORD ptr[rsp+64]
	movdqu xmm9,XMMWORD ptr[rsp+48]
	movdqu xmm8,XMMWORD ptr[rsp+32]
	movdqu xmm7,XMMWORD ptr[rsp+16]
	movdqu xmm6,XMMWORD ptr[rsp]
	add rsp,144

	pop rbx
	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_Warp2_8_SSE3 endp	


;JPSDR_Sobel_8_SSE2 proc psrc:dword,pdst:dword,src_pitch:dword,y_:dword,height:dword,i_:dword,thresh:dword
; psrc = rcx
; pdst = rdx
; src_pitch = r8d
; y_ = r9d
JPSDR_Sobel_8_SSE2 proc public frame

height equ dword ptr[rbp+48]
i_ equ dword ptr[rbp+56]
thresh equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog

    mov rsi,rcx
    mov rdi,rdx
    movsxd rdx,r8d
    xor	rax,rax
    movsxd rcx,r9d
    test rcx,rcx
    cmovnz rax,rdx
    inc	rcx
    add	rdx,rax
    cmp	ecx,height     ;32 bit O.K.
    cmovz rdx,rax
    sub rsi,rax
	mov r8,16
	mov r9,4
	movsxd rcx,i_
    sub	rdi,r8
    sub	rdi,rsi
    movd xmm0,thresh
    pshufd xmm0,xmm0,0
    packssdw xmm0,xmm0
    packuswb xmm0,xmm0
	
JPSDR_Sobel_8_SSE2_1:
    movdqu xmm2,XMMWORD ptr[rsi-1]
    movdqu xmm3,XMMWORD ptr[rsi]
    movdqu xmm4,XMMWORD ptr[rsi+1]
    movdqu xmm5,XMMWORD ptr[rsi+rdx-1]
    movdqu xmm6,XMMWORD ptr[rsi+rdx]
    movdqu xmm7,XMMWORD ptr[rsi+rdx+1]

    movdqa xmm1,xmm2
    pavgb xmm1,xmm4
    pavgb xmm3,xmm1

    movdqa xmm1,xmm5
    pavgb xmm1,xmm7
    pavgb xmm6,xmm1

    movdqa xmm1,xmm3
    psubusb xmm3,xmm6
    psubusb xmm6,xmm1
    por xmm6,xmm3

    movdqu xmm1,XMMWORD ptr[rsi+rax-1]
    movdqu xmm3,XMMWORD ptr[rsi+rax+1]
    pavgb xmm5,xmm2
    pavgb xmm7,xmm4
    pavgb xmm1,xmm5
    pavgb xmm3,xmm7
    movdqa xmm5,xmm1
    psubusb xmm1,xmm3
    psubusb xmm3,xmm5
    por xmm1,xmm3

    movdqa xmm2,xmm6
    paddusb xmm2,xmm1
    pmaxub xmm1,xmm6
    paddusb xmm2,xmm1

    movdqa xmm3,xmm2
    paddusb xmm2,xmm2
    paddusb xmm2,xmm3
    paddusb xmm2,xmm2
    pminub xmm2,xmm0  ;thresh
    add	rsi,r8
    sub	rcx,r9
    jb short JPSDR_Sobel_8_SSE2_2
    movntdq XMMWORD ptr[rsi+rdi],xmm2
    jnz JPSDR_Sobel_8_SSE2_1
    jmp short JPSDR_Sobel_8_SSE2_Fin

JPSDR_Sobel_8_SSE2_2:
    test rcx,2
    jz short JPSDR_Sobel_8_SSE2_3
    movq qword ptr[rsi+rdi],xmm2
    test rcx,1
    jz short JPSDR_Sobel_8_SSE2_Fin
    add rsi,8
    psrldq xmm2,8

JPSDR_Sobel_8_SSE2_3:
    movd dword ptr[rsi+rdi],xmm2

JPSDR_Sobel_8_SSE2_Fin:
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_Sobel_8_SSE2 endp


;JPSDR_H_BlurR6_8_SSE3 proc psrc2:dword,ptmp2:dword,src_row_size:dword,dq0toF:dword
; psrc2 = rcx
; ptmp2 = rdx
; src_row_size = r8d
; dq0toF = r9
JPSDR_H_BlurR6_8_SSE3 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog

	mov rax,16
    mov rsi,rcx
    mov rdi,rdx
    movsxd rcx,r8d
    add rsi,rax
    sub rdi,rsi
    movdqa xmm6,XMMWORD ptr[rsi-16]
    movdqa xmm5,xmm6
    movdqa xmm7,xmm6
    pxor xmm0,xmm0
    pshufb xmm5,xmm0
    sub rcx,rax
    jna JPSDR_H_BlurR6_8_SSE3_2
    
JPSDR_H_BlurR6_8_SSE3_1:
    movdqa xmm7,XMMWORD ptr[rsi]
    movdqa xmm0,xmm6
    movdqa xmm2,xmm7
    palignr xmm0,xmm5,10
    palignr xmm2,xmm6,6
    pavgb xmm0,xmm2
    movdqa xmm3,xmm6
    movdqa xmm4,xmm7
    palignr xmm3,xmm5,11
    palignr xmm4,xmm6,5
    pavgb xmm3,xmm4
    pavgb xmm0,xmm3
    movdqa xmm1,xmm6
    movdqa xmm2,xmm7
    palignr xmm1,xmm5,12
    palignr xmm2,xmm6,4
    pavgb xmm1,xmm2
    movdqa xmm3,xmm6
    movdqa xmm4,xmm7
    palignr xmm3,xmm5,13
    palignr xmm4,xmm6,3
    pavgb xmm3,xmm4
    pavgb xmm1,xmm3
    pavgb xmm0,xmm1
    movdqa xmm1,xmm6
    movdqa xmm2,xmm7
    palignr xmm1,xmm5,14
    palignr xmm2,xmm6,2
    pavgb xmm1,xmm2
    movdqa xmm3,xmm6
    movdqa xmm4,xmm7
    palignr xmm3,xmm5,15
    palignr xmm4,xmm6,1
    pavgb xmm3,xmm4
    pavgb xmm1,xmm3
    pavgb xmm1,xmm6
    movdqa xmm5,xmm6
    movdqa xmm6,xmm7
    pavgb xmm0,xmm1
    pavgb xmm0,xmm1
    movntdq XMMWORD ptr[rsi+rdi],xmm0
    add rsi,rax
    sub rcx,rax
    ja JPSDR_H_BlurR6_8_SSE3_1
	
JPSDR_H_BlurR6_8_SSE3_2:
    add	rcx,15
    pxor xmm0,xmm0
    movd xmm1,ecx
    pshufb xmm1,xmm0
    pminub xmm1,XMMWORD ptr[r9]	;0x0F0E..00
    pshufb xmm6,xmm1
    psrldq xmm7,15
    pshufb xmm7,xmm0
    movdqa xmm0,xmm6
    movdqa xmm2,xmm7
    palignr xmm0,xmm5,10
    palignr xmm2,xmm6,6
    pavgb xmm0,xmm2
    movdqa xmm3,xmm6
    movdqa xmm4,xmm7
    palignr xmm3,xmm5,11
    palignr xmm4,xmm6,5
    pavgb xmm3,xmm4
    pavgb xmm0,xmm3
    movdqa xmm1,xmm6
    movdqa xmm2,xmm7
    palignr xmm1,xmm5,12
    palignr xmm2,xmm6,4
    pavgb xmm1,xmm2
    movdqa xmm3,xmm6
    movdqa xmm4,xmm7
    palignr xmm3,xmm5,13
    palignr xmm4,xmm6,3
    pavgb xmm3,xmm4
    pavgb xmm1,xmm3
    pavgb xmm0,xmm1
    movdqa xmm1,xmm6
    movdqa xmm2,xmm7
    palignr xmm1,xmm5,14
    palignr xmm2,xmm6,2
    pavgb xmm1,xmm2
    movdqa xmm3,xmm6
    movdqa xmm4,xmm7
    palignr xmm3,xmm5,15
    palignr xmm4,xmm6,1
    pavgb xmm3,xmm4
    pavgb xmm1,xmm3
    pavgb xmm1,xmm6
    movdqa xmm5,xmm6
    movdqa xmm6,xmm7
    pavgb xmm0,xmm1
    pavgb xmm0,xmm1
    movntdq XMMWORD ptr[rsi+rdi],xmm0
	
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_H_BlurR6_8_SSE3 endp


;JPSDR_H_BlurR6_8_SSE2 proc proc psrc2:dword,ptmp2:dword,src_row_size_16:dword
; psrc2 = rcx
; ptmp2 = rdx
; src_row_size_16 = r8d
JPSDR_H_BlurR6_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog

    mov rsi,rcx
    mov rdi,rdx
    movsxd rcx,r8d
    sub rdi,rsi
	mov rax,16

JPSDR_H_BlurR6_8_SSE2_1:
    movdqu xmm6,XMMWORD ptr[rsi-6]
    movdqu xmm0,XMMWORD ptr[rsi+6]
    pavgb xmm6,xmm0
    movdqu xmm5,XMMWORD ptr[rsi-5]
    movdqu xmm7,XMMWORD ptr[rsi+5]
    pavgb xmm5,xmm7
    movdqu xmm4,XMMWORD ptr[rsi-4]
    movdqu xmm0,XMMWORD ptr[rsi+4]
    pavgb xmm4,xmm0
    movdqu xmm3,XMMWORD ptr[rsi-3]
    movdqu xmm7,XMMWORD ptr[rsi+3]
    pavgb xmm3,xmm7
    movdqu xmm2,XMMWORD ptr[rsi-2]
    movdqu xmm0,XMMWORD ptr[rsi+2]
    pavgb xmm2,xmm0
    movdqu xmm1,XMMWORD ptr[rsi-1]
    movdqu xmm7,XMMWORD ptr[rsi+1]
    pavgb xmm1,xmm7
    movdqa xmm0,XMMWORD ptr[rsi]
    pavgb xmm6,xmm5
    pavgb xmm4,xmm3
    pavgb xmm2,xmm1
    pavgb xmm6,xmm4
    pavgb xmm2,xmm0
    pavgb xmm6,xmm2
    pavgb xmm6,xmm2
    movntdq XMMWORD ptr[rsi+rdi],xmm6
    add rsi,rax
	loop JPSDR_H_BlurR6_8_SSE2_1
	
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_H_BlurR6_8_SSE2 endp	


;JPSDR_V_BlurR6a_8_SSE2 proc psrc2:dword,ptmp2:dword,tmp_pitch:dword,src_row_size_16:dword
; psrc2 = rcx
; ptmp2 = rdx
; tmp_pitch = r8d
; src_row_size_16 = r9d
JPSDR_V_BlurR6a_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,16
	.allocstack 16
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	.endprolog
	
	movsxd rax,r8d
	mov	rsi,rdx
	mov	rdi,rcx
	movsxd rcx,r9d
	lea	r8,[rax+2*rax]
	lea	rdx,[r8+2*rax]
	add	rdx,rsi
	sub	rdi,rsi
	mov r9,16
	
JPSDR_V_BlurR6a_8_SSE2_1:
	movdqa xmm0,XMMWORD ptr[rsi]
	movdqa xmm1,XMMWORD ptr[rsi+rax]
	movdqa xmm2,XMMWORD ptr[rsi+2*rax]
	movdqa xmm3,XMMWORD ptr[rsi+r8]
	movdqa xmm4,XMMWORD ptr[rsi+4*rax]
	movdqa xmm5,XMMWORD ptr[rdx]
	movdqa xmm6,XMMWORD ptr[rdx+rax]
	pavgb xmm6,xmm5
	pavgb xmm4,xmm3
	pavgb xmm2,xmm1
	pavgb xmm6,xmm4
	pavgb xmm2,xmm0
	pavgb xmm6,xmm2
	pavgb xmm6,xmm2
	movntdq XMMWORD ptr[rsi+rdi],xmm6
	add rsi,r9
	add rdx,r9
	loop JPSDR_V_BlurR6a_8_SSE2_1
	
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,16

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_V_BlurR6a_8_SSE2 endp	


;JPSDR_V_BlurR6b_8_SSE2 proc psrc2:dword,ptmp2:dword,tmp_pitch:dword,src_row_size_16:dword
; psrc2 = rcx
; ptmp2 = rdx
; tmp_pitch = r8d
; src_row_size_16 = r9d
JPSDR_V_BlurR6b_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,16
	.allocstack 16
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	.endprolog
	
	movsxd rax,r8d
	mov rsi,rdx
	mov rdi,rcx
	movsxd rcx,r9d
	lea r8,[rax+2*rax]
	lea rdx,[r8+2*rax]
	lea r9,[rsi+2*rdx]
	add rdx,rsi
	sub rdi,rsi
	mov r10,16
	
JPSDR_V_BlurR6b_8_SSE2_1:
	movdqa xmm6,XMMWORD ptr[rsi]
	pavgb xmm6,XMMWORD ptr[r9+2*rax]
	movdqa xmm5,XMMWORD ptr[rsi+rax]
	pavgb xmm5,XMMWORD ptr[r9+rax]
	movdqa xmm4,XMMWORD ptr[rsi+2*rax]
	pavgb xmm4,XMMWORD ptr[r9]
	movdqa xmm3,XMMWORD ptr[rsi+r8]
	pavgb xmm3,XMMWORD ptr[rdx+4*rax]
	movdqa xmm2,XMMWORD ptr[rsi+4*rax]
	pavgb xmm2,XMMWORD ptr[rsi+8*rax]
	movdqa xmm1,XMMWORD ptr[rdx]
	pavgb xmm1,XMMWORD ptr[rdx+2*rax]
	movdqa xmm0,XMMWORD ptr[rdx+rax]
	pavgb xmm6,xmm5
	pavgb xmm4,xmm3
	pavgb xmm2,xmm1
	pavgb xmm6,xmm4
	pavgb xmm2,xmm0
	pavgb xmm6,xmm2
	pavgb xmm6,xmm2
	
	movntdq XMMWORD ptr[rsi+rdi],xmm6
	add rsi,r10
	add rdx,r10
	add r9,r10
	loop JPSDR_V_BlurR6b_8_SSE2_1
	
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,16

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_V_BlurR6b_8_SSE2 endp	


;JPSDR_V_BlurR6c_8_SSE2 proc psrc2:dword,ptmp2:dword,tmp_pitch:dword,src_row_size_16:dword
; psrc2 = rcx
; ptmp2 = rdx
; tmp_pitch = r8d
; src_row_size_16 = r9d
JPSDR_V_BlurR6c_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,16
	.allocstack 16
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	.endprolog
	
	movsxd rax,r8d
	mov	rsi,rdx
	mov	rdi,rcx
	movsxd rcx,r9d
	lea	r8,[rax+2*rax]
	lea	rdx,[r8+2*rax]
	add	rdx,rsi
	sub	rdi,rsi
	mov r9,16
	
JPSDR_V_BlurR6c_8_SSE2_1:
	movdqa xmm6,XMMWORD ptr[rsi]
	movdqa xmm5,XMMWORD ptr[rsi+rax]
	movdqa xmm4,XMMWORD ptr[rsi+2*rax]
	movdqa xmm3,XMMWORD ptr[rsi+r8]
	movdqa xmm2,XMMWORD ptr[rsi+4*rax]
	movdqa xmm1,XMMWORD ptr[rdx]
	movdqa xmm0,XMMWORD ptr[rdx+rax]
	pavgb xmm6,xmm5
	pavgb xmm4,xmm3
	pavgb xmm2,xmm1
	pavgb xmm6,xmm4
	pavgb xmm2,xmm0
	pavgb xmm6,xmm2
	pavgb xmm6,xmm2
	movntdq XMMWORD ptr[rsi+rdi],xmm6
	add rsi,r9
	add rdx,r9
	loop JPSDR_V_BlurR6c_8_SSE2_1
	
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,16

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_V_BlurR6c_8_SSE2 endp


;JPSDR_H_BlurR2_8_SSE3 proc psrc2:dword,ptmp2:dword,i_:dword,dq0toF:dword
; psrc2 = rcx
; ptmp2 = rdx
; i_ = r8d
; dq0toF = r9
JPSDR_H_BlurR2_8_SSE3 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,32
	.allocstack 32
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	movdqa XMMWORD ptr[rsp+16],xmm7
	.savexmm128 xmm7,16
	.endprolog

	movdqa xmm4,XMMWORD ptr[r9]
	mov rax,16
	mov rsi,rcx
	mov rdi,rdx
	movsxd rcx,r8d
	add rsi,rax
	sub rdi,rsi
	movdqa xmm6,XMMWORD ptr[rsi-16]
	movdqa xmm5,xmm6
	movdqa xmm7,xmm6
	pxor xmm0,xmm0
	pshufb xmm5,xmm0
	sub rcx,rax
	jna short JPSDR_H_BlurR2_8_SSE3_2
	
JPSDR_H_BlurR2_8_SSE3_1:
	movdqa xmm7,XMMWORD ptr[rsi]
	movdqa xmm0,xmm6
	movdqa xmm2,xmm7
	palignr xmm0,xmm5,14
	palignr xmm2,xmm6,2
	pavgb xmm0,xmm2
	movdqa xmm1,xmm6
	movdqa xmm3,xmm7
	palignr xmm1,xmm5,15
	palignr xmm3,xmm6,1
	pavgb xmm0,xmm6
	pavgb xmm1,xmm3
	pavgb xmm0,xmm6
	movdqa xmm5,xmm6
	movdqa xmm6,xmm7
	pavgb xmm0,xmm1
	movntdq XMMWORD ptr[rsi+rdi],xmm0
	add rsi,rax
	sub rcx,rax
	ja short JPSDR_H_BlurR2_8_SSE3_1
	
JPSDR_H_BlurR2_8_SSE3_2:
	add rcx,15
	pxor xmm0,xmm0
	movd xmm1,ecx
	pshufb xmm1,xmm0
	pminub xmm1,xmm4 ;0x0F0E..00
	pshufb xmm6,xmm1
	psrldq xmm7,15
	pshufb xmm7,xmm0
	movdqa xmm0,xmm6
	movdqa xmm2,xmm7
	palignr xmm0,xmm5,14
	palignr xmm2,xmm6,2
	pavgb xmm0,xmm2
	movdqa xmm1,xmm6
	movdqa xmm3,xmm7
	palignr xmm1,xmm5,15
	palignr xmm3,xmm6,1
	pavgb xmm0,xmm6
	pavgb xmm1,xmm3
	pavgb xmm0,xmm6
	movdqa xmm5,xmm6
	movdqa xmm6,xmm7
	pavgb xmm0,xmm1
	movntdq XMMWORD ptr[rsi+rdi],xmm0
	
	movdqa xmm7,XMMWORD ptr[rsp+16]
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,32

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_H_BlurR2_8_SSE3 endp


;JPSDR_H_BlurR2_8_SSE2 proc psrc2:dword,ptmp2:dword,ia:dword
; psrc2 = rcx
; ptmp2 = rdx
; ia = r8d
JPSDR_H_BlurR2_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	.endprolog

	mov rax,16	
	mov rsi,rcx
	mov rdi,rdx
	movsxd rcx,r8d
	add rsi,rcx
	add rdi,rcx
	neg rcx
	
JPSDR_H_BlurR2_8_SSE2_1:
	movdqu xmm0,XMMWORD ptr[rcx+rsi-2]
	movdqu xmm4,XMMWORD ptr[rcx+rsi+2]
	pavgb xmm0,xmm4
	movdqu xmm1,XMMWORD ptr[rcx+rsi-1]
	movdqu xmm3,XMMWORD ptr[rcx+rsi+1]
	pavgb xmm1,xmm3
	movdqa xmm2,XMMWORD ptr[rcx+rsi]
	pavgb xmm0,xmm2
	pavgb xmm0,xmm2
	pavgb xmm0,xmm1
	movntdq XMMWORD ptr[rcx+rdi],xmm0
	add	rcx,rax
	jnz	short JPSDR_H_BlurR2_8_SSE2_1
	
	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_H_BlurR2_8_SSE2 endp
	

;JPSDR_V_BlurR2_8_SSE2 proc psrc2:dword,ptmp2:dword,ia:dword,tmp_pitchp1:dword,tmp_pitchp2:dword,tmp_pitchn1:dword,tmp_pitchn2:dword,
; psrc2 = rcx
; ptmp2 = rdx
; ia = r8d
; tmp_pitchp1 = r9d
JPSDR_V_BlurR2_8_SSE2 proc public frame

tmp_pitchp2 equ dword ptr[rbp+48]
tmp_pitchn1 equ dword ptr[rbp+56]
tmp_pitchn2 equ dword ptr[rbp+64]

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	.endprolog
	
	mov rsi,rdx
	mov rdi,rcx
	movsxd rcx,r8d
	movsxd rax,tmp_pitchp2
	movsxd rdx,tmp_pitchn2
	movsxd r8,r9d
	movsxd r9,tmp_pitchn1
	sub rdi,rsi
	mov r10,16
	
JPSDR_V_BlurR2_8_SSE2_1:
	movdqa xmm0,XMMWORD ptr[rsi+rax]
	pavgb xmm0,XMMWORD ptr[rsi+rdx]
	movdqa xmm1,XMMWORD ptr[rsi+r8]
	pavgb xmm1,XMMWORD ptr[rsi+r9]
	movdqa xmm2,XMMWORD ptr[rsi]
	pavgb xmm0,xmm2
	pavgb xmm0,xmm2
	pavgb xmm0,xmm1
	movntdq XMMWORD ptr[rsi+rdi],xmm0
	add	rsi,r10
	sub	rcx,r10
	jnz	short JPSDR_V_BlurR2_8_SSE2_1
		  
	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_V_BlurR2_8_SSE2 endp


;JPSDR_GuideChroma1_8_SSE2 proc py:dword,pu:dword,pitch_y:dword,width_uv_8:dword
; py = rcx
; pu = rdx
; pitch_y = r8d
; width_uv_8 = r9d
JPSDR_GuideChroma1_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,16
	.allocstack 16
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	.endprolog

	mov rsi,rcx
	movsxd rax,r8d
	movsxd rcx,r9d
	mov r8,16
	sub rsi,rcx
	sub rsi,rcx
	add rax,rsi
	sub rdx,rcx
	sub rdx,r8
	pcmpeqw xmm6,xmm6
	psrlw xmm6,8
	
JPSDR_GuideChroma1_8_SSE2_1:
	movdqa xmm0,XMMWORD ptr[rsi+2*rcx]
	movdqa xmm2,XMMWORD ptr[rsi+2*rcx+16]
	movdqa xmm1,xmm0
	movdqa xmm3,xmm2
	pand xmm0,xmm6
	pand xmm2,xmm6
	packuswb xmm0,xmm2
	psrlw xmm1,8
	psrlw xmm3,8
	packuswb xmm1,xmm3
	pavgb xmm0,xmm1
	movdqa xmm1,XMMWORD ptr[rax+2*rcx]
	movdqa xmm3,XMMWORD ptr[rax+2*rcx+16]
	movdqa xmm2,xmm1
	movdqa xmm4,xmm3
	pand xmm1,xmm6
	pand xmm3,xmm6
	packuswb xmm1,xmm3
	psrlw xmm2,8
	psrlw xmm4,8
	packuswb xmm2,xmm4
	pavgb xmm1,xmm2
	pavgb xmm0,xmm1
	add rcx,r8
	jg short JPSDR_GuideChroma1_8_SSE2_2
	movntdq XMMWORD ptr[rcx+rdx],xmm0
	jnz short JPSDR_GuideChroma1_8_SSE2_1
	jmp short JPSDR_GuideChroma1_8_SSE2_3
	
JPSDR_GuideChroma1_8_SSE2_2:
	movq qword ptr[rcx+rdx],xmm0
			
JPSDR_GuideChroma1_8_SSE2_3:
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,16

	pop rdi
	pop rsi
	pop rbp
	
	ret
JPSDR_GuideChroma1_8_SSE2 endp


;JPSDR_GuideChroma2_8_SSE2 proc py:dword,pu:dword,width_uv_8:dword
; py = rcx
; pu = rdx
; width_uv_8 = r8d
JPSDR_GuideChroma2_8_SSE2 proc public frame

	push rbp
	.pushreg rbp
	mov rbp,rsp
	push rsi
	.pushreg rsi
	push rdi
	.pushreg rdi
	sub rsp,16
	.allocstack 16
	movdqa XMMWORD ptr[rsp],xmm6
	.savexmm128 xmm6,0
	.endprolog
	
	mov rax,16
	mov rsi,rcx
	movsxd rcx,r8d
	sub rsi,rcx
	sub rsi,rcx
	sub rdx,rcx
	sub rdx,rax
	pcmpeqw xmm6,xmm6
	psrlw xmm6,8
	
JPSDR_GuideChroma2_8_SSE2_1:
	movdqa xmm0,XMMWORD ptr[rsi+2*rcx]
	movdqa xmm2,XMMWORD ptr[rsi+2*rcx+16]
	movdqa xmm1,xmm0
	movdqa xmm3,xmm2
	pand xmm0,xmm6
	pand xmm2,xmm6
	packuswb xmm0,xmm2
	psrlw xmm1,8
	psrlw xmm3,8
	packuswb xmm1,xmm3
	pavgb xmm0,xmm1
	add rcx,rax
	jg short JPSDR_GuideChroma2_8_SSE2_2
	movntdq XMMWORD ptr[rcx+rdx],xmm0
	jnz short JPSDR_GuideChroma2_8_SSE2_1
	jmp short JPSDR_GuideChroma2_8_SSE2_3
	
JPSDR_GuideChroma2_8_SSE2_2:
	movq qword ptr[rcx+rdx],xmm0
	
JPSDR_GuideChroma2_8_SSE2_3:	
	movdqa xmm6,XMMWORD ptr[rsp]
	add rsp,16

	pop rdi
	pop rsi
	pop rbp

	ret
JPSDR_GuideChroma2_8_SSE2 endp


end





