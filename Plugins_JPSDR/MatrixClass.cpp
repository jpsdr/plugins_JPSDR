// Test_Calcul_Matrice.cpp : définit le point d'entrée pour l'application console.
//

#include "MatrixClass.h"

#include <memory.h>

#if _MSC_VER >= 1900
#define AVX2_BUILD_POSSIBLE
#endif

extern "C" void VectorProductF_SSE2(const float *coeff_a,const float *coeff_x,float *result,uint16_t lght);
extern "C" void VectorProductD_SSE2(const double *coeff_a,const double *coeff_x,double *result,uint16_t lght);
extern "C" void CoeffProductF_SSE2(const float *coeff_a,const float *coeff_b,float *coeff_c,uint16_t lght);
extern "C" void CoeffAddProductF_SSE2(const float *coeff_a,const float *coeff_b,float *coeff_c,uint16_t lght);
extern "C" void CoeffProductD_SSE2(const double *coeff_a,const double *coeff_b,double *coeff_c,uint16_t lght);
extern "C" void CoeffAddProductD_SSE2(const double *coeff_a,const double *coeff_b,double *coeff_c,uint16_t lght);
extern "C" void VectorProductF_AVX(const float *coeff_a,const float *coeff_x,float *result,uint16_t lght);
extern "C" void VectorProductD_AVX(const double *coeff_a,const double *coeff_x,double *result,uint16_t lght);
extern "C" void CoeffProductF_AVX(const float *coeff_a,const float *coeff_b,float *coeff_c,uint16_t lght);
extern "C" void CoeffAddProductF_AVX(const float *coeff_a,const float *coeff_b,float *coeff_c,uint16_t lght);
extern "C" void CoeffProductD_AVX(const double *coeff_a,const double *coeff_b,double *coeff_c,uint16_t lght);
extern "C" void CoeffAddProductD_AVX(const double *coeff_a,const double *coeff_b,double *coeff_c,uint16_t lght);

#define MATRIX_ALIGN_SIZE 64
#define MATRIX_ALIGN_SHIFT 6

Vector::Vector(void)
{
	Coeff=NULL;
	length=0;
	size=0;
	data_type=DATA_NONE;
}


Vector::Vector(const uint16_t l,const COEFF_DATA_TYPE data)
{
	Coeff=NULL;
	length=0;
	size=0;
	data_type=DATA_NONE;

	if (l==0) return;

	size_t coeff_size;

	switch(data)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return;

	const size_t p0=((((size_t)l*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc(p0,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return;

	size=p0;
	length=l;
	data_type=data;

	const size_t n0=(size_t)l*coeff_size,n=p0-n0;

	if (n>0) memset(((uint8_t *)Coeff)+n0,0,n);
}


Vector::Vector(const Vector &x)
{
	Coeff=NULL;
	length=0;
	size=0;
	data_type=DATA_NONE;

	if (&x==NULL) return;

	const uint16_t l=x.length;

	if ((x.Coeff==NULL) || (l==0)) return;

	size_t coeff_size;

	switch(x.data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return;

	const size_t p0=((((size_t)l*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc(p0,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return;

	size=p0;
	length=l;
	data_type=x.data_type;

	const size_t n0=(size_t)l*coeff_size,n=p0-n0;

	if (n>0) memset(((uint8_t *)Coeff)+n0,0,n);

	CopyStrict(x);
}


Vector::~Vector(void)
{
	Destroy();
}


bool Vector::Create(void)
{
	if ((Coeff!=NULL) || (length==0)) return(false);

	size_t coeff_size;

	switch(data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const size_t p0=((((size_t)length*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc(p0,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=p0;

	const size_t n0=(size_t)length*coeff_size,n=p0-n0;

	if (n>0) memset(((uint8_t *)Coeff)+n0,0,n);

	return(true);
}


bool Vector::Create(const uint16_t l,const COEFF_DATA_TYPE data)
{
	if ((Coeff!=NULL) || (l==0)) return(false);

	size_t coeff_size;

	switch(data)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const size_t p0=((((size_t)l*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc(p0,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=p0;
	length=l;
	data_type=data;

	const size_t n0=(size_t)l*coeff_size,n=p0-n0;

	if (n>0) memset(((uint8_t *)Coeff)+n0,0,n);

	return(true);
}


bool Vector::Create(const Vector &x)
{
	if ((Coeff!=NULL) || (&x==NULL)) return(false);

	const uint16_t l=x.length;

	if ((x.Coeff==NULL) || (l==0)) return(false);

	size_t coeff_size;

	switch(x.data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const size_t p0=((((size_t)l*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc(p0,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=p0;
	length=l;
	data_type=x.data_type;

	const size_t n0=(size_t)l*coeff_size,n=p0-n0;

	if (n>0) memset(((uint8_t *)Coeff)+n0,0,n);

	CopyStrict(x);

	return(true);
}


void Vector::Destroy(void)
{
	if (Coeff!=NULL)
	{
		_aligned_free(Coeff);
		Coeff=NULL;
	}
	length=0;
	size=0;
	data_type=DATA_NONE;
}


bool Vector::CopyStrict(const Vector &x)
{
	if ((Coeff==NULL) || (&x==NULL) || (length==0)) return(false);

	const uint16_t l=x.length;

	if ((x.Coeff==NULL) || (l==0) || (l!=length) || (x.data_type!=data_type)) return(false);

	size_t coeff_size;

	switch(data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const size_t size_line=(size_t)l*coeff_size;

	memcpy(Coeff,x.GetPtrVector(),size_line);

	return(true);
}


bool Vector::FillD(const double data)
{
	if ((Coeff==NULL) || (length==0)) return(false);

	double *a=(double *)Coeff;

	for (uint16_t i=0; i<length; i++)
		*a++=data;

	return(true);
}


bool Vector::FillF(const float data)
{
	if ((Coeff==NULL) || (length==0)) return(false);

	float *a=(float *)Coeff;

	for (uint16_t i=0; i<length; i++)
		*a++=data;

	return(true);
}


bool Vector::FillZero(void)
{
	if ((Coeff==NULL) || (length==0)) return(false);

	memset(Coeff,0,size);

	return(true);
}


bool Vector::SetInfo(const uint16_t l,const COEFF_DATA_TYPE data)
{
	if ((Coeff!=NULL) || (length!=0) || (l==0) || (data_type==DATA_NONE)) return(false);

	length=l; data_type=data;

	return(true);
}


void Vector::GetInfo(uint16_t &l,COEFF_DATA_TYPE &data) const
{
	l=length; data=data_type;
}


bool Vector::GetSafeD(const uint16_t i,double &d) const
{
	if ((Coeff==NULL) || (length==0) || (i>=length) || (data_type!=DATA_DOUBLE)) return(false);

	d=((double *)Coeff)[i];

	return(true);
}


bool Vector::SetSafeD(const uint16_t i,const double d)
{
	if ((Coeff==NULL) || (length==0) || (i>=length) || (data_type!=DATA_DOUBLE)) return(false);

	((double *)Coeff)[i]=d;

	return(true);
}


bool Vector::GetSafeF(const uint16_t i,float &d) const
{
	if ((Coeff==NULL) || (length==0) || (i>=length) || (data_type!=DATA_FLOAT)) return(false);

	d=((float *)Coeff)[i];

	return(true);
}


bool Vector::SetSafeF(const uint16_t i,const float d)
{
	if ((Coeff==NULL) || (length==0) || (i>=length) || (data_type!=DATA_FLOAT)) return(false);

	((float *)Coeff)[i]=d;

	return(true);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Vector_Compute::Vector_Compute(void)
{
	SSE2_Enable=false;
	AVX_Enable=false;
	AVX2_Enable=false;
}


Vector_Compute::~Vector_Compute(void)
{
}


Vector_Compute::Vector_Compute(const uint16_t l,const COEFF_DATA_TYPE data):Vector(l,data)
{
	SSE2_Enable=false;
	AVX_Enable=false;
	AVX2_Enable=false;
}


Vector_Compute::Vector_Compute(const Vector_Compute &x):Vector(x)
{
	if (&x!=NULL)
	{
		SSE2_Enable=x.SSE2_Enable;
		AVX_Enable=x.AVX_Enable;
		AVX2_Enable=x.AVX2_Enable;
	}
	else
	{
		SSE2_Enable=false;
		AVX_Enable=false;
		AVX2_Enable=false;
	}
}


bool Vector_Compute::Product_AX(const Matrix &ma, const Vector &x)
{
	const uint16_t l=length;

	if ((Coeff==NULL) || (l==0) || (&ma==NULL) || (&x==NULL)) return(false);
	if (!ma.AllocCheck() || !x.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns(),lb=x.GetLength();

	if ((ca!=lb) || (l!=la) || (ma.GetDataType()!=x.GetDataType()) || (x.GetDataType()!=data_type)) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : ProductF_AX(ma,x); break;
		case DATA_DOUBLE : ProductD_AX(ma,x); break;
		default : return(false);
	}

	return(true);
}


bool Vector_Compute::Product_AX(const Matrix &ma)
{
	const uint16_t l=length;

	if ((Coeff==NULL) || (l==0) || (&ma==NULL)) return(false);
	if (!ma.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns();

	if ((ca!=la) || (l!=la) || (ma.GetDataType()!=data_type)) return(false);

	Vector b(*this);

	if (!b.AllocCheck()) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : ProductF_AX(ma,b); break;
		case DATA_DOUBLE : ProductD_AX(ma,b); break;
		default : return(false);
	}

	return(true);
}


void Vector_Compute::ProductF_AX(const Matrix &ma, const Vector &x)
{
	const uint16_t l=length,lb=x.GetLength();
	const uint8_t *a0=(uint8_t *)ma.GetPtrMatrix();
	const float *x1=(float *)x.GetPtrVector();
	float *c1=(float *)Coeff;
	const ptrdiff_t pa=ma.GetPitch();

	if (AVX_Enable)
	{
		const uint16_t n=(lb+7)>>3;

		for (int32_t i=0; i<l; i++)
		{
			VectorProductF_AVX((float *)a0,x1,c1++,n);
			a0+=pa;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(lb+3)>>2;

			for (int32_t i=0; i<l; i++)
			{
				VectorProductF_SSE2((float *)a0,x1,c1++,n);
				a0+=pa;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				const float *a1=(float *)a0;
				float s=0.0f;

				for (uint16_t k=0; k<lb; k++)
					s+=a1[k]*x1[k];
				*c1++=s;
				a0+=pa;
			}
		}
	}
}


void Vector_Compute::ProductD_AX(const Matrix &ma, const Vector &x)
{
	const uint16_t l=length,lb=x.GetLength();
	const uint8_t *a0=(uint8_t *)ma.GetPtrMatrix();
	const double *x1=(double *)x.GetPtrVector();
	double *c1=(double *)Coeff;
	const ptrdiff_t pa=ma.GetPitch();

	if (AVX_Enable)
	{
		const uint16_t n=(lb+3)>>2;

		for (int32_t i=0; i<l; i++)
		{
			VectorProductD_AVX((double *)a0,x1,c1++,lb);
			a0+=pa;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(lb+1)>>1;

			for (int32_t i=0; i<l; i++)
			{
				VectorProductD_SSE2((double *)a0,x1,c1++,lb);
				a0+=pa;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				const double *a1=(double *)a0;
				double s=0.0;

				for (uint16_t k=0; k<lb; k++)
					s+=a1[k]*x1[k];
				*c1++=s;
				a0+=pa;
			}
		}
	}
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Matrix::Matrix(void)
{
	Coeff=NULL;
	columns=0; lines=0;
	size=0;
	pitch=0;
	data_type=DATA_NONE;
}


Matrix::Matrix(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data)
{
	Coeff=NULL;
	columns=0; lines=0;
	size=0;
	pitch=0;
	data_type=DATA_NONE;

	if ((c==0) || (l==0)) return;

	size_t coeff_size;

	switch(data)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return;

	const ptrdiff_t p0=((((ptrdiff_t)c*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc((size_t)p0*(size_t)l,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return;

	size=(size_t)p0*(size_t)l;
	pitch=p0;
	columns=c; lines=l;
	data_type=data;

	const size_t n0=(size_t)c*coeff_size,n=(size_t)p0-n0;

	if (n>0) 
	{
		uint8_t *a=(uint8_t *)Coeff;

		for(uint16_t i=0; i<l; i++)
		{
			memset(a+n0,0,n);
			a+=p0;
		}
	}
}


Matrix::Matrix(const Matrix &m)
{
	Coeff=NULL;
	columns=0; lines=0;
	size=0;
	pitch=0;
	data_type=DATA_NONE;

	if (&m==NULL) return;

	const uint16_t c=m.columns,l=m.lines;

	if ((m.Coeff==NULL) || (c==0) || (l==0)) return;

	size_t coeff_size;

	switch(m.data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return;

	const ptrdiff_t p0=((((ptrdiff_t)c*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc((size_t)p0*(size_t)l,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return;

	size=(size_t)p0*(size_t)l;
	pitch=p0;
	columns=c; lines=l;
	data_type=m.data_type;

	const size_t n0=(size_t)c*coeff_size,n=(size_t)p0-n0;

	if (n>0) 
	{
		uint8_t *a=(uint8_t *)Coeff;

		for(uint16_t i=0; i<l; i++)
		{
			memset(a+n0,0,n);
			a+=p0;
		}
	}

	CopyStrict(m);
}


Matrix::~Matrix(void)
{
	Destroy();
}


bool Matrix::Create(void)
{
	if ((Coeff!=NULL) || (columns==0) || (lines==0)) return(false);

	size_t coeff_size;

	switch(data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const ptrdiff_t p0=((((ptrdiff_t)columns*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc((size_t)p0*(size_t)lines,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=(size_t)p0*(size_t)lines;
	pitch=p0;

	const size_t n0=(size_t)columns*coeff_size,n=(size_t)p0-n0;

	if (n>0) 
	{
		uint8_t *a=(uint8_t *)Coeff;

		for(uint16_t i=0; i<lines; i++)
		{
			memset(a+n0,0,n);
			a+=p0;
		}
	}

	return(true);
}


bool Matrix::Create(const Matrix &m)
{
	if ((Coeff!=NULL) || (&m==NULL)) return(false);

	const uint16_t c=m.columns,l=m.lines;

	if ((m.Coeff==NULL) || (c==0) || (l==0)) return(false);

	size_t coeff_size;

	switch(m.data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const ptrdiff_t p0=((((ptrdiff_t)c*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc((size_t)p0*(size_t)l,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=(size_t)p0*(size_t)l;
	pitch=p0;
	columns=c; lines=l;
	data_type=m.data_type;

	const size_t n0=(size_t)c*coeff_size,n=(size_t)p0-n0;

	if (n>0) 
	{
		uint8_t *a=(uint8_t *)Coeff;

		for(uint16_t i=0; i<l; i++)
		{
			memset(a+n0,0,n);
			a+=p0;
		}
	}

	return(CopyStrict(m));
}


bool Matrix::Create(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data)
{
	if ((Coeff!=NULL) || (c==0) || (l==0)) return(false);

	size_t coeff_size;

	switch(data)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const ptrdiff_t p0=((((ptrdiff_t)c*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc((size_t)p0*(size_t)l,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=(size_t)p0*(size_t)l;
	pitch=p0;
	columns=c; lines=l;

	const size_t n0=(size_t)c*coeff_size,n=(size_t)p0-n0;

	if (n>0) 
	{
		uint8_t *a=(uint8_t *)Coeff;

		for(uint16_t i=0; i<l; i++)
		{
			memset(a+n0,0,n);
			a+=p0;
		}
	}

	return(true);
}


void Matrix::Destroy(void)
{
	if (Coeff!=NULL)
	{
		_aligned_free(Coeff);
		Coeff=NULL;
	}
	columns=0; lines=0;
	size=0;
	pitch=0;
	data_type=DATA_NONE;
}


bool Matrix::FillD(const double data)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (c==0) || (l==0) || (data_type!=DATA_DOUBLE)) return(false);

	uint8_t *a=(uint8_t *)Coeff;

	for (uint16_t i=0; i<l; i++)
	{
		double *a0=(double *)a;

		for (uint16_t j=0; j<c; j++)
			*a0++=data;

		a+=pitch;
	}

	return(true);
}


bool Matrix::FillF(const float data)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (c==0) || (l==0) || (data_type!=DATA_FLOAT)) return(false);

	uint8_t *a=(uint8_t *)Coeff;

	for (uint16_t i=0; i<l; i++)
	{
		float *a0=(float *)a;

		for (uint16_t j=0; j<c; j++)
			*a0++=data;

		a+=pitch;
	}

	return(true);
}


bool Matrix::FillZero(void)
{
	if ((Coeff==NULL) || (columns==0) || (lines==0)) return(false);

	memset(Coeff,0,size);

	return(true);
}


bool Matrix::SetInfo(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data)
{
	if ((Coeff!=NULL) || (columns!=0) || (lines!=0) || (c==0) || (l==0) || (data_type==DATA_NONE)) return(false);

	columns=c; lines=l; data_type=data;

	return(true);
}


void Matrix::GetInfo(uint16_t &l,uint16_t &c,COEFF_DATA_TYPE &data) const
{
	c=columns; l=lines; data=data_type;
}


bool Matrix::GetSafeD(const uint16_t i,const uint16_t j,double &d) const
{
	if ((Coeff==NULL) || (columns==0) || (lines==0) || (i>=lines) || (j>=columns) || (data_type!=DATA_DOUBLE)) return(false);

	d=((double *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j];

	return(true);
}


bool Matrix::SetSafeD(const uint16_t i,const uint16_t j,const double d)
{
	if ((Coeff==NULL) || (columns==0) || (lines==0) || (i>=lines) || (j>=columns) || (data_type!=DATA_DOUBLE)) return(false);

	((double *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j]=d;

	return(true);
}


bool Matrix::GetSafeF(const uint16_t i,const uint16_t j,float &d) const
{
	if ((Coeff==NULL) || (columns==0) || (lines==0) || (i>=lines) || (j>=columns) || (data_type!=DATA_FLOAT)) return(false);

	d=((float *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j];

	return(true);
}


bool Matrix::SetSafeF(const uint16_t i,const uint16_t j,const float d)
{
	if ((Coeff==NULL) || (columns==0) || (lines==0) || (i>=lines) || (j>=columns) || (data_type!=DATA_FLOAT)) return(false);

	((float *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j]=d;

	return(true);
}


bool Matrix::CopyStrict(const Matrix &m)
{
	if ((Coeff==NULL) || (&m==NULL) || (columns==0) || (lines==0)) return(false);

	const uint16_t c=m.columns,l=m.lines;

	if ((m.Coeff==NULL) || (c==0) || (l==0) || (c!=columns) || (l!=lines) || (m.data_type!=data_type)) return(false);

	size_t coeff_size;

	switch(data_type)
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const ptrdiff_t pa=m.pitch,p=pitch;
	const uint8_t *a=(uint8_t *)m.Coeff;
	uint8_t *b=(uint8_t *)Coeff;
	const size_t size_line=(size_t)c*coeff_size;

	for (uint16_t i=0; i<l; i++)
	{
		memcpy(b,a,size_line);
		a+=pa;
		b+=p;
	}

	return(true);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Matrix_Compute::Matrix_Compute(void)
{
	zero_value=0.0;
	SSE2_Enable=false;
	AVX_Enable=false;
	AVX2_Enable=false;
}


Matrix_Compute::~Matrix_Compute(void)
{
}


Matrix_Compute::Matrix_Compute(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data):Matrix(l,c,data)
{
	zero_value=0.0;
	SSE2_Enable=false;
	AVX_Enable=false;
	AVX2_Enable=false;
}


Matrix_Compute::Matrix_Compute(const Matrix_Compute &m):Matrix(m)
{
	if (&m!=NULL)
	{
		zero_value=m.zero_value;
		SSE2_Enable=m.SSE2_Enable;
		AVX_Enable=m.AVX_Enable;
		AVX2_Enable=m.AVX2_Enable;
	}
	else
	{
		zero_value=0.0;
		SSE2_Enable=false;
		AVX_Enable=false;
		AVX2_Enable=false;
	}
}


bool Matrix_Compute::CopyStrict(const Matrix_Compute &m)
{
	if (!Matrix::CopyStrict(m)) return(false);

	zero_value=m.zero_value;

	return(true);
}


bool Matrix_Compute::CreateTranspose(const Matrix &m)
{
	if ((Coeff!=NULL) || (&m==NULL)) return(false);

	const uint16_t l=m.GetColumns(),c=m.GetLines();

	if ((m.GetPtrMatrix()==NULL) || (c==0) || (l==0)) return(false);

	size_t coeff_size;

	switch(m.GetDataType())
	{
		case DATA_FLOAT : coeff_size=sizeof(float); break;
		case DATA_DOUBLE : coeff_size=sizeof(double); break;
		case DATA_UINT64 : coeff_size=sizeof(uint64_t); break;
		case DATA_INT64 : coeff_size=sizeof(int64_t); break;
		case DATA_UINT32 : coeff_size=sizeof(uint32_t); break;
		case DATA_INT32 : coeff_size=sizeof(int32_t); break;
		case DATA_UINT16 : coeff_size=sizeof(uint16_t); break;
		case DATA_INT16 : coeff_size=sizeof(int16_t); break;
		case DATA_UINT8 : coeff_size=sizeof(uint8_t); break;
		case DATA_INT8 : coeff_size=sizeof(int8_t); break;
		default : coeff_size=0; break;
	}
	if (coeff_size==0) return(false);

	const ptrdiff_t p0=((((ptrdiff_t)c*coeff_size)+MATRIX_ALIGN_SIZE-1) >> MATRIX_ALIGN_SHIFT) << MATRIX_ALIGN_SHIFT;

	Coeff=(void *)_aligned_malloc((size_t)p0*(size_t)l,MATRIX_ALIGN_SIZE);
	if (Coeff==NULL) return(false);

	size=(size_t)p0*(size_t)l;
	pitch=p0;
	columns=c; lines=l;
	data_type=m.GetDataType();

	const size_t n0=(size_t)c*coeff_size,n=(size_t)p0-n0;

	if (n>0) 
	{
		uint8_t *a=(uint8_t *)Coeff;

		for(uint16_t i=0; i<l; i++)
		{
			memset(a+n0,0,n);
			a+=p0;
		}
	}

	switch(data_type)
	{
		case DATA_FLOAT : TransposeF(m); break;
		case DATA_DOUBLE : TransposeD(m); break;
		case DATA_UINT64 : TransposeU64(m); break;
		case DATA_INT64 : TransposeI64(m); break;
		case DATA_UINT32 : TransposeU32(m); break;
		case DATA_INT32 : TransposeI32(m); break;
		case DATA_UINT16 : TransposeU16(m); break;
		case DATA_INT16 : TransposeI16(m); break;
		case DATA_UINT8 : TransposeU8(m); break;
		case DATA_INT8 : TransposeI8(m); break;
		default : return(false);
	}

	return(true);
}


bool Matrix_Compute::Product_AB(const Matrix &ma, const Matrix &mb)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (&ma==NULL) || (&mb==NULL)) return(false);
	if (!ma.AllocCheck() || !mb.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns(),lb=mb.GetLines(),cb=mb.GetColumns();

	if ((ca!=lb) || (c!=cb) || (l!=la) || (ma.GetDataType()!=mb.GetDataType()) || (ma.GetDataType()!=data_type)) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : ProductF_AB(ma,mb); break;
		case DATA_DOUBLE : ProductD_AB(ma,mb); break;
		default : return(false);
	}

	return(true);
}


void Matrix_Compute::ProductF_AB(const Matrix &ma, const Matrix &mb)
{
	const uint16_t li=lines,co=columns;
	const uint16_t ca=ma.GetColumns();
	const uint8_t *a=(uint8_t *)ma.GetPtrMatrix();
	const uint8_t *b=(uint8_t *)mb.GetPtrMatrix();
	uint8_t *c=(uint8_t *)Coeff;
	const ptrdiff_t pa=ma.GetPitch(),pb=mb.GetPitch(),pc=pitch;
	const size_t size_line=(size_t)co*sizeof(float);

	const uint8_t *a0=a;
	uint8_t *c0=c;

	if (AVX_Enable)
	{
		const uint16_t n=(co+7)>>3;

		for (uint16_t i=0; i<li; i++)
		{
			const float a1=*(float *)a0;

			if (a1!=0.0) CoeffProductF_AVX(&a1,(float *)b,(float *)c0,n);
			else memset(c0,0,size_line);

			a0+=pa;
			c0+=pc;
		}
		b+=pb;

		for(uint16_t i=1; i<ca; i++)
		{
			const float *b1=(float *)b;

			a0=a;c0=c;
			for (uint16_t j=0; j<li; j++)
			{
				const float a1=((float *)a0)[i];

				if (a1!=0.0) CoeffAddProductF_AVX(&a1,b1,(float *)c0,n);
				a0+=pa;
				c0+=pc;
			}
			b+=pb;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(co+3)>>2;

			for (uint16_t i=0; i<li; i++)
			{
				const float a1=*(float *)a0;

				if (a1!=0.0) CoeffProductF_SSE2(&a1,(float *)b,(float *)c0,n);
				else memset(c0,0,size_line);

				a0+=pa;
				c0+=pc;
			}
			b+=pb;

			for(uint16_t i=1; i<ca; i++)
			{
				const float *b1=(float *)b;

				a0=a;c0=c;
				for (uint16_t j=0; j<li; j++)
				{
					const float a1=((float *)a0)[i];

					if (a1!=0.0) CoeffAddProductF_SSE2(&a1,b1,(float *)c0,n);
					a0+=pa;
					c0+=pc;
				}
				b+=pb;
			}
		}
		else
		{
			for (uint16_t i=0; i<li; i++)
			{
				const float *b1=(float *)b;
				const float a1=*(float *)a0;
				float *c1=(float *)c0;

				if (a1!=0.0)
				{
					for (uint16_t j=0; j<co; j++)
						c1[j]=a1*b1[j];
				}
				else memset(c1,0,size_line);

				a0+=pa;
				c0+=pc;
			}
			b+=pb;

			for(uint16_t i=1; i<ca; i++)
			{
				const float *b1=(float *)b;

				a0=a;c0=c;
				for (uint16_t j=0; j<li; j++)
				{
					float *c1=(float *)c0;
					const float a1=((float *)a0)[i];

					if (a1!=0.0)
					{
						for(uint16_t k=0; k<co; k++)
							c1[k]+=a1*b1[k];
					}
					a0+=pa;
					c0+=pc;
				}
				b+=pb;
			}
		}
	}
}


void Matrix_Compute::ProductD_AB(const Matrix &ma, const Matrix &mb)
{
	const uint16_t li=lines,co=columns;
	const uint16_t ca=ma.GetColumns();
	const uint8_t *a=(uint8_t *)ma.GetPtrMatrix();
	const uint8_t *b=(uint8_t *)mb.GetPtrMatrix();
	uint8_t *c=(uint8_t *)Coeff;
	const ptrdiff_t pa=ma.GetPitch(),pb=mb.GetPitch(),pc=pitch;
	const size_t size_line=(size_t)co*sizeof(double);

	const uint8_t *a0=a;
	uint8_t *c0=c;

	if (AVX_Enable)
	{
		const uint16_t n=(co+3)>>2;

		for (uint16_t i=0; i<li; i++)
		{
			const double a1=*(double *)a0;

			if (a1!=0.0) CoeffProductD_AVX(&a1,(double *)b,(double *)c0,n);
			else memset(c0,0,size_line);

			a0+=pa;
			c0+=pc;
		}
		b+=pb;

		for(uint16_t i=1; i<ca; i++)
		{
			const double *b1=(double *)b;

			a0=a;c0=c;
			for (uint16_t j=0; j<li; j++)
			{
				const double a1=((double *)a0)[i];

				if (a1!=0.0) CoeffAddProductD_AVX(&a1,b1,(double *)c0,n);
				a0+=pa;
				c0+=pc;
			}
			b+=pb;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(co+1)>>1;

			for (uint16_t i=0; i<li; i++)
			{
				const double a1=*(double *)a0;

				if (a1!=0.0) CoeffProductD_SSE2(&a1,(double *)b,(double *)c0,n);
				else memset(c0,0,size_line);

				a0+=pa;
				c0+=pc;
			}
			b+=pb;

			for(uint16_t i=1; i<ca; i++)
			{
				const double *b1=(double *)b;

				a0=a;c0=c;
				for (uint16_t j=0; j<li; j++)
				{
					const double a1=((double *)a0)[i];

					if (a1!=0.0) CoeffAddProductD_SSE2(&a1,b1,(double *)c0,n);
					a0+=pa;
					c0+=pc;
				}
				b+=pb;
			}
		}
		else
		{
			for (uint16_t i=0; i<li; i++)
			{
				const double *b1=(double *)b;
				const double a1=*(double *)a0;
				double *c1=(double *)c0;

				if (a1!=0.0)
				{
					for (uint16_t j=0; j<co; j++)
						c1[j]=a1*b1[j];
				}
				else memset(c1,0,size_line);

				a0+=pa;
				c0+=pc;
			}
			b+=pb;

			for(uint16_t i=1; i<ca; i++)
			{
				const double *b1=(double *)b;

				a0=a;c0=c;
				for (uint16_t j=0; j<li; j++)
				{
					double *c1=(double *)c0;
					const double a1=((double *)a0)[i];

					if (a1!=0.0)
					{
						for(uint16_t k=0; k<co; k++)
							c1[k]+=a1*b1[k];
					}
					a0+=pa;
					c0+=pc;
				}
				b+=pb;
			}
		}
	}
}


bool Matrix_Compute::Product_AtB(const Matrix &ma,const Matrix &mb)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (&ma==NULL) || (&mb==NULL)) return(false);
	if (!ma.AllocCheck() || !mb.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns(),lb=mb.GetLines(),cb=mb.GetColumns();

	if ((ca!=cb) || (c!=lb) || (l!=la) || (ma.GetDataType()!=mb.GetDataType()) || (ma.GetDataType()!=data_type)) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : ProductF_AtB(ma,mb); break;
		case DATA_DOUBLE : ProductD_AtB(ma,mb); break;
		default : return(false);
	}

	return(true);
}


void Matrix_Compute::ProductF_AtB(const Matrix &ma,const Matrix &mb)
{
	const uint16_t l=lines,co=columns;
	const uint16_t ca=ma.GetColumns();
	const uint8_t *a=(uint8_t *)ma.GetPtrMatrix();
	const uint8_t *b=(uint8_t *)mb.GetPtrMatrix();
	uint8_t *c=(uint8_t *)Coeff;
	const ptrdiff_t pa=ma.GetPitch(),pb=mb.GetPitch(),pc=pitch;

	if (AVX_Enable)
	{
		const uint16_t n=(ca+7)>>3;

		for (uint16_t i=0; i<l; i++)
		{
			const uint8_t *b0=b;
			float *c1=(float *)c;

			for (uint16_t j=0; j<co; j++)
			{
				VectorProductF_AVX((float *)a,(float *)b0,c1++,n);
				b0+=pb;
			}
			a+=pa;
			c+=pc;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(ca+3)>>2;

			for (uint16_t i=0; i<l; i++)
			{
				const uint8_t *b0=b;
				float *c1=(float *)c;

				for (uint16_t j=0; j<co; j++)
				{
					VectorProductF_SSE2((float *)a,(float *)b0,c1++,n);
					b0+=pb;
				}
				a+=pa;
				c+=pc;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				const float *a1=(float *)a;
				const uint8_t *b0=b;
				float *c1=(float *)c;

				for (uint16_t j=0; j<co; j++)
				{
					float s=0.0f;
					const float *b1=(float *)b0;

					for (uint16_t k=0; k<ca; k++)
						s+=a1[k]*b1[k];
					*c1++=s;
					b0+=pb;
				}
				a+=pa;
				c+=pc;
			}
		}
	}
}


void Matrix_Compute::ProductD_AtB(const Matrix &ma,const Matrix &mb)
{
	const uint16_t l=lines,co=columns;
	const uint16_t ca=ma.GetColumns();
	const uint8_t *a=(uint8_t *)ma.GetPtrMatrix();
	const uint8_t *b=(uint8_t *)mb.GetPtrMatrix();
	uint8_t *c=(uint8_t *)Coeff;
	const ptrdiff_t pa=ma.GetPitch(),pb=mb.GetPitch(),pc=pitch;

	if (AVX_Enable)
	{
		const uint16_t n=(ca+3)>>2;

		for (uint16_t i=0; i<l; i++)
		{
			const uint8_t *b0=b;
			double *c1=(double *)c;

			for (uint16_t j=0; j<co; j++)
			{
				VectorProductD_AVX((double *)a,(double *)b0,c1++,n);
				b0+=pb;
			}
			a+=pa;
			c+=pc;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(ca+1)>>1;

			for (uint16_t i=0; i<l; i++)
			{
				const uint8_t *b0=b;
				double *c1=(double *)c;

				for (uint16_t j=0; j<co; j++)
				{
					VectorProductD_SSE2((double *)a,(double *)b0,c1++,n);
					b0+=pb;
				}
				a+=pa;
				c+=pc;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				const double *a1=(double *)a;
				const uint8_t *b0=b;
				double *c1=(double *)c;

				for (uint16_t j=0; j<co; j++)
				{
					double s=0.0f;
					const double *b1=(double *)b0;

					for (uint16_t k=0; k<ca; k++)
						s+=a1[k]*b1[k];
					*c1++=s;
					b0+=pb;
				}
				a+=pa;
				c+=pc;
			}
		}
	}
}


bool Matrix_Compute::Product_tAA(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (&ma==NULL)) return(false);
	if (!ma.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns();

	if ((c!=ca) || (l!=ca) || (ma.GetDataType()!=data_type)) return(false);

	Matrix_Compute b;

	b.CreateTranspose(ma);

	if (!b.AllocCheck()) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : ProductF_AB(b,ma); break;
		case DATA_DOUBLE : ProductD_AB(b,ma); break;
		default : return(false);
	}

	return(true);
}


bool Matrix_Compute::Product_tAA(void)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (l!=c)) return(false);

	Matrix_Compute a(*this),b;

	if (!a.AllocCheck()) return(false);

	b.CreateTranspose(*this);

	if (!b.AllocCheck()) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : ProductF_AB(b,a); break;
		case DATA_DOUBLE : ProductF_AB(b,a); break;
		default : return(false);
	}

	return(true);
}


bool Matrix_Compute::Inverse(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (&ma==NULL)) return(false);
	if (!ma.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns();

	if ((ca!=la) || (l!=la) || (c!=ca) || (ma.GetDataType()!=data_type)) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : return(InverseF(ma)); break;
		case DATA_DOUBLE : return(InverseD(ma)); break;
		default : return(false);
	}

	return(false);
}


bool Matrix_Compute::Inverse(void)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (c!=l)) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : return(InverseF(*this)); break;
		case DATA_DOUBLE : return(InverseD(*this)); break;
		default : return(false);
	}

	return(false);
}


bool Matrix_Compute::InverseF(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	const uint16_t c2=c<<1;

	Matrix b(l,c2,data_type);

	if (!b.AllocCheck()) return(false);

	b.FillZero();

	const uint8_t *a0=(uint8_t *)ma.GetPtrMatrix();
	uint8_t *b_=(uint8_t *)b.GetPtrMatrix();
	const size_t size_line=(size_t)c*sizeof(float);
	const ptrdiff_t pa=ma.GetPitch(),pb=b.GetPitch(),pc=pitch;
	uint8_t *b0=b_;

	for (uint16_t i=0; i<l; i++)
	{
		memcpy(b0,a0,size_line);
		b0+=pb;
		a0+=pa;
		b.SetF(i,i+c,1.0f);
	}

	b0=b_;

	if (AVX_Enable)
	{
		const uint16_t n=(c2+7)>>3;

		for (uint16_t i=0; i<l; i++)
		{
			float *b2=(float *)b0;
			uint8_t *b1=b_;

			for (uint16_t j=0; j<c; j++)
			{
				float *b3=(float *)b1;

				if (i!=j)
				{
					const float ratio=-b3[i]/b2[i];

					if (ratio!=0.0) CoeffAddProductF_AVX(&ratio,b2,b3,n);
				}
				b1+=pb;
			}
			b0+=pb;
		}

		b0=b_;
		for (uint16_t i=0; i<l; i++)
		{
			float *b2=(float *)b0;
			const float a=1.0f/b2[i];

			CoeffProductF_AVX(&a,b2,b2,n);

			b0+=pb;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(c2+3)>>2;

			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					float *b3=(float *)b1;

					if (i!=j)
					{
						const float ratio=-b3[i]/b2[i];

						if (ratio!=0.0) CoeffAddProductF_SSE2(&ratio,b2,b3,n);
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;
				const float a=1.0f/b2[i];

				CoeffProductF_SSE2(&a,b2,b2,n);

				b0+=pb;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					float *b3=(float *)b1;

					if (i!=j)
					{
						const float ratio=-b3[i]/b2[i];

						if (ratio!=0.0)
						{
							for (uint16_t k=0; k<c2; k++)
								b3[k]+=ratio*b2[k];
						}
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;
				const float a=1.0f/b2[i];

				for (uint16_t j=0; j<c2; j++)
					b2[j]*=a;
				b0+=pb;
			}
		}
	}

	b0=b_;
	b0+=(ptrdiff_t)c*sizeof(float);
	uint8_t *c0=(uint8_t *)Coeff;
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(c0,b0,size_line);
		b0+=pb;
		c0+=pc;
	}

	return(true);
}


bool Matrix_Compute::InverseD(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	const uint16_t c2=c<<1;

	Matrix b(l,c2,data_type);

	if (!b.AllocCheck()) return(false);

	b.FillZero();

	const uint8_t *a0=(uint8_t *)ma.GetPtrMatrix();
	uint8_t *b_=(uint8_t *)b.GetPtrMatrix();
	const size_t size_line=(size_t)c*sizeof(double);
	const ptrdiff_t pa=ma.GetPitch(),pb=b.GetPitch(),pc=pitch;
	uint8_t *b0=b_;
	
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(b0,a0,size_line);
		b0+=pb;
		a0+=pa;
		b.SetF(i,i+c,1.0);
	}

	b0=b_;

	if (AVX_Enable)
	{
		const uint16_t n=(c2+3)>>2;

		for (uint16_t i=0; i<l; i++)
		{
			double *b2=(double *)b0;
			uint8_t *b1=b_;

			for (uint16_t j=0; j<c; j++)
			{
				double *b3=(double *)b1;

				if (i!=j)
				{
					const double ratio=-b3[i]/b2[i];

					if (ratio!=0.0) CoeffAddProductD_AVX(&ratio,b2,b3,n);
				}
				b1+=pb;
			}
			b0+=pb;
		}

		b0=b_;
		for (uint16_t i=0; i<l; i++)
		{
			double *b2=(double *)b0;
			const double a=1.0f/b2[i];

			CoeffProductD_AVX(&a,b2,b2,n);

			b0+=pb;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(c2+1)>>1;

			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					double *b3=(double *)b1;

					if (i!=j)
					{
						const double ratio=-b3[i]/b2[i];

						if (ratio!=0.0) CoeffAddProductD_SSE2(&ratio,b2,b3,n);
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;
				const double a=1.0f/b2[i];

				CoeffProductD_SSE2(&a,b2,b2,n);

				b0+=pb;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					double *b3=(double *)b1;

					if (i!=j)
					{
						const double ratio=-b3[i]/b2[i];

						if (ratio!=0.0)
						{
							for (uint16_t k=0; k<c2; k++)
								b3[k]+=ratio*b2[k];
						}
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;
				const double a=1.0f/b2[i];

				for (uint16_t j=0; j<c2; j++)
					b2[j]*=a;
				b0+=pb;
			}
		}
	}

	b0=b_;
	b0+=(ptrdiff_t)c*sizeof(double);
	uint8_t *c0=(uint8_t *)Coeff;
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(c0,b0,size_line);
		b0+=pb;
		c0+=pc;
	}

	return(true);
}


/*
Return :
 0 : Matrix is reversed.
 -1 : Allocation/Matrix configuration error.
 -2 : Matrix can't be reversed.
*/
int8_t Matrix_Compute::InverseSafe(const Matrix_Compute &ma)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (&ma==NULL)) return(-1);
	if (!ma.AllocCheck()) return(-1);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns();

	if ((ca!=la) || (l!=la) || (c!=ca) || (ma.GetDataType()!=data_type)) return(-1);

	switch(data_type)
	{
		case DATA_FLOAT : return(InverseSafeF(ma)); break;
		case DATA_DOUBLE : return(InverseSafeD(ma)); break;
		default : return(-1);
	}

	return(-1);
}


int8_t Matrix_Compute::InverseSafe(void)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (c!=l)) return(-1);

	switch(data_type)
	{
		case DATA_FLOAT : return(InverseSafeF(*this)); break;
		case DATA_DOUBLE : return(InverseSafeD(*this)); break;
		default : return(-1);
	}

	return(-1);
}


/*
Return :
 0 : Matrix is reversed.
 -1 : Allocation/Matrix configuration error.
 -2 : Matrix can't be reversed.
*/
int8_t Matrix_Compute::InverseSafeF(const Matrix_Compute &ma)
{
	const uint16_t l=lines,c=columns;
	const uint16_t c2=c<<1;

	Matrix b(l,c2,data_type);

	if (!b.AllocCheck()) return(-1);

	b.FillZero();

	const uint8_t *a0=(uint8_t *)ma.GetPtrMatrix();
	uint8_t *b_=(uint8_t *)b.GetPtrMatrix();
	const size_t size_line=(size_t)c*sizeof(float);
	const ptrdiff_t pa=ma.GetPitch(),pb=b.GetPitch(),pc=pitch;
	const float _zero=(float)ma.zero_value;
	uint8_t *b0=b_;
	
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(b0,a0,size_line);
		b0+=pb;
		a0+=pa;
		b.SetF(i,i+l,1.0f);
	}

	b0=b_;

	if (AVX_Enable)
	{
		const uint16_t n=(c2+7)>>3;

		for (uint16_t i=0; i<l; i++)
		{
			float *b2=(float *)b0;
			uint8_t *b1=b_;

			for (uint16_t j=0; j<c; j++)
			{
				float *b3=(float *)b1;

				if (i!=j)
				{
					if (fabs(b2[i])<=_zero) return(-2);

					const float ratio=-b3[i]/b2[i];

					if (ratio!=0.0) CoeffAddProductF_AVX(&ratio,b2,b3,n);
				}
				b1+=pb;
			}
			b0+=pb;
		}

		b0=b_;
		for (uint16_t i=0; i<l; i++)
		{
			float *b2=(float *)b0;

			if (fabs(b2[i])<=_zero) return(-2);

			const float a=1.0f/b2[i];

			CoeffProductF_AVX(&a,b2,b2,n);

			b0+=pb;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(c2+3)>>2;

			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					float *b3=(float *)b1;

					if (i!=j)
					{
						if (fabs(b2[i])<=_zero) return(-2);

						const float ratio=-b3[i]/b2[i];

						if (ratio!=0.0) CoeffAddProductF_SSE2(&ratio,b2,b3,n);
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;

				if (fabs(b2[i])<=_zero) return(-2);

				const float a=1.0f/b2[i];

				CoeffProductF_SSE2(&a,b2,b2,n);

				b0+=pb;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					float *b3=(float *)b1;

					if (i!=j)
					{
						if (fabs(b2[i])<=_zero) return(-2);

						const float ratio=-b3[i]/b2[i];

						if (ratio!=0.0)
						{
							for (uint16_t k=0; k<c2; k++)
								b3[k]+=ratio*b2[k];
						}
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				float *b2=(float *)b0;

				if (fabs(b2[i])<=_zero) return(-2);

				const float a=1.0f/b2[i];

				for (uint16_t j=0; j<c2; j++)
					b2[j]*=a;
				b0+=pb;
			}
		}
	}

	b0=b_;
	b0+=(ptrdiff_t)c*sizeof(float);
	uint8_t *c0=(uint8_t *)Coeff;
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(c0,b0,size_line);
		b0+=pb;
		c0+=pc;
	}

	return(0);
}


/*
Return :
 0 : Matrix is reversed.
 -1 : Allocation/Matrix configuration error.
 -2 : Matrix can't be reversed.
*/
int8_t Matrix_Compute::InverseSafeD(const Matrix_Compute &ma)
{
	const uint16_t l=lines,c=columns;
	const uint16_t c2=c<<1;

	Matrix b(l,c2,data_type);

	if (!b.AllocCheck()) return(-1);

	b.FillZero();

	const uint8_t *a0=(uint8_t *)ma.GetPtrMatrix();
	uint8_t *b_=(uint8_t *)b.GetPtrMatrix();
	const size_t size_line=(size_t)c*sizeof(double);
	const ptrdiff_t pa=ma.GetPitch(),pb=b.GetPitch(),pc=pitch;
	const double _zero=ma.zero_value;
	uint8_t *b0=b_;
	
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(b0,a0,size_line);
		b0+=pb;
		a0+=pa;
		b.SetF(i,i+l,1.0);
	}

	b0=b_;

	if (AVX_Enable)
	{
		const uint16_t n=(c2+3)>>2;

		for (uint16_t i=0; i<l; i++)
		{
			double *b2=(double *)b0;
			uint8_t *b1=b_;

			for (uint16_t j=0; j<c; j++)
			{
				double *b3=(double *)b1;

				if (i!=j)
				{
					if (fabs(b2[i])<=_zero) return(-2);

					const double ratio=-b3[i]/b2[i];

					if (ratio!=0.0) CoeffAddProductD_AVX(&ratio,b2,b3,n);
				}
				b1+=pb;
			}
			b0+=pb;
		}

		b0=b_;
		for (uint16_t i=0; i<l; i++)
		{
			double *b2=(double *)b0;

			if (fabs(b2[i])<=_zero) return(-2);

			const double a=1.0f/b2[i];

			CoeffProductD_AVX(&a,b2,b2,n);

			b0+=pb;
		}
	}
	else
	{
		if (SSE2_Enable)
		{
			const uint16_t n=(c2+1)>>1;

			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					double *b3=(double *)b1;

					if (i!=j)
					{
						if (fabs(b2[i])<=_zero) return(-2);

						const double ratio=-b3[i]/b2[i];

						if (ratio!=0.0) CoeffAddProductD_SSE2(&ratio,b2,b3,n);
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;

				if (fabs(b2[i])<=_zero) return(-2);

				const double a=1.0f/b2[i];

				CoeffProductD_SSE2(&a,b2,b2,n);

				b0+=pb;
			}
		}
		else
		{
			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;
				uint8_t *b1=b_;

				for (uint16_t j=0; j<c; j++)
				{
					double *b3=(double *)b1;

					if (i!=j)
					{
						if (fabs(b2[i])<=_zero) return(-2);

						const double ratio=-b3[i]/b2[i];

						if (ratio!=0.0)
						{
							for (uint16_t k=0; k<c2; k++)
								b3[k]+=ratio*b2[k];
						}
					}
					b1+=pb;
				}
				b0+=pb;
			}

			b0=b_;
			for (uint16_t i=0; i<l; i++)
			{
				double *b2=(double *)b0;

				if (fabs(b2[i])<=_zero) return(-2);

				const double a=1.0f/b2[i];

				for (uint16_t j=0; j<c2; j++)
					b2[j]*=a;
				b0+=pb;
			}
		}
	}

	b0=b_;
	b0+=(ptrdiff_t)c*sizeof(double);
	uint8_t *c0=(uint8_t *)Coeff;
	for (uint16_t i=0; i<l; i++)
	{
		memcpy(c0,b0,size_line);
		b0+=pb;
		c0+=pc;
	}

	return(0);
}


bool Matrix_Compute::Transpose(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (&ma==NULL)) return(false);
	if (!ma.AllocCheck()) return(false);

	const uint16_t la=ma.GetLines(),ca=ma.GetColumns();

	if ((l!=ca) || (c!=la) || (ma.GetDataType()!=data_type)) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : TransposeF(ma); break;
		case DATA_DOUBLE : TransposeD(ma); break;
		case DATA_UINT64 : TransposeU64(ma); break;
		case DATA_INT64 : TransposeI64(ma); break;
		case DATA_UINT32 : TransposeU32(ma); break;
		case DATA_INT32 : TransposeI32(ma); break;
		case DATA_UINT16 : TransposeU16(ma); break;
		case DATA_INT16 : TransposeI16(ma); break;
		case DATA_UINT8 : TransposeU8(ma); break;
		case DATA_INT8 : TransposeI8(ma); break;
		default : return(false);
	}

	return(true);
}

bool Matrix_Compute::Transpose(void)
{
	const uint16_t l=lines,c=columns;

	if ((Coeff==NULL) || (l==0) || (c==0) || (l!=c)) return(false);

	Matrix b(*this);

	if (!b.AllocCheck()) return(false);

	switch(data_type)
	{
		case DATA_FLOAT : TransposeF(b); break;
		case DATA_DOUBLE : TransposeD(b); break;
		case DATA_UINT64 : TransposeU64(b); break;
		case DATA_INT64 : TransposeI64(b); break;
		case DATA_UINT32 : TransposeU32(b); break;
		case DATA_INT32 : TransposeI32(b); break;
		case DATA_UINT16 : TransposeU16(b); break;
		case DATA_INT16 : TransposeI16(b); break;
		case DATA_UINT8 : TransposeU8(b); break;
		case DATA_INT8 : TransposeI8(b); break;
		default : return(false);
	}

	return(true);
}


void Matrix_Compute::TransposeF(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(float);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		float *a=(float *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(float *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeD(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(double);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		double *a=(double *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(double *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeU64(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(uint64_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		uint64_t *a=(uint64_t *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(uint64_t *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeI64(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(int64_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		int64_t *a=(int64_t *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(int64_t *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeU32(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(uint32_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		uint32_t *a=(uint32_t *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(uint32_t *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeI32(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(int32_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		int32_t *a=(int32_t *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(int32_t *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeU16(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(uint16_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		uint16_t *a=(uint16_t *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(uint16_t *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeI16(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(int16_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		int16_t *a=(int16_t *)a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*(int16_t *)b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeU8(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	uint8_t *a0=(uint8_t *)Coeff;
	const uint8_t *b0=(uint8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(uint8_t);

	for (uint16_t i=0; i<l; i++)
	{
		const uint8_t *b1=b0;
		uint8_t *a=a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}


void Matrix_Compute::TransposeI8(const Matrix &ma)
{
	const uint16_t l=lines,c=columns;
	int8_t *a0=(int8_t *)Coeff;
	const int8_t *b0=(int8_t *)ma.GetPtrMatrix();
	const ptrdiff_t pb=ma.GetPitch(),p=pitch,db0=sizeof(int8_t);

	for (uint16_t i=0; i<l; i++)
	{
		const int8_t *b1=b0;
		int8_t *a=a0;

		for (uint16_t j=0; j<c; j++)
		{
			*a++=*b1;
			b1+=pb;
		}
		a0+=p;
		b0+=db0;
	}
}
