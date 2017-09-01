#ifndef _MATRIX_CLASS_H
#define _MATRIX_CLASS_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef enum COEFF_DATA_TYPE_ {DATA_NONE,DATA_DOUBLE,DATA_FLOAT,DATA_UINT64,DATA_INT64,
	DATA_UINT32,DATA_INT32,DATA_UINT16,DATA_INT16,DATA_UINT8,DATA_INT8} COEFF_DATA_TYPE;


class Vector
{
public :
	Vector(void);
	Vector(const uint16_t l,const COEFF_DATA_TYPE data);
	Vector(const Vector &x);
	virtual ~Vector(void);

	bool AllocCheck(void) const {return(Coeff!=NULL);}
	bool Create(void);
	bool Create(const uint16_t l,const COEFF_DATA_TYPE data);
	bool Create(const Vector &x);
	bool CopyStrict(const Vector &x);
	void Destroy(void);
	bool FillD(const double data);
	bool FillF(const float data);
	bool FillZero(void);
	COEFF_DATA_TYPE GetDataType(void) const {return(data_type);}
	bool SetInfo(const uint16_t l,const COEFF_DATA_TYPE data);
	void GetInfo(uint16_t &l,COEFF_DATA_TYPE &data) const;
	uint16_t GetLength(void) const {return(length);}
	void* GetPtrVector(void) const {return(Coeff);}
	size_t GetDataSize(void) const {return(size);}
	double GetD(const uint16_t i) const {return(((double *)Coeff)[i]);}
	float GetF(const uint16_t i) const {return(((float *)Coeff)[i]);}
	void SetD(const uint16_t i,const double d) {((double *)Coeff)[i]=d;}
	void SetF(const uint16_t i,const float d) {((float *)Coeff)[i]=d;}
	bool GetSafeD(const uint16_t i,double &d) const ;
	bool SetSafeD(const uint16_t i,const double d);
	bool GetSafeF(const uint16_t i,float &d) const ;
	bool SetSafeF(const uint16_t i,const float d);

protected :
	void *Coeff;
	uint16_t length;
	size_t size;
	COEFF_DATA_TYPE data_type;
};

class Matrix;

class Vector_Compute : public Vector
{
public :
	Vector_Compute(void);
	Vector_Compute(const uint16_t l,const COEFF_DATA_TYPE data);
	Vector_Compute(const Vector_Compute &x);
	virtual ~Vector_Compute(void);

	void EnableSSE2(void) {SSE2_Enable=true;}
	void DisableSSE2(void) {SSE2_Enable=false;}
	void EnableAVX(void) {AVX_Enable=true;}
	void DisableAVX(void) {AVX_Enable=false;}
	void EnableAVX2(void) {AVX2_Enable=true;}
	void DisableAVX2(void) {AVX2_Enable=false;}

	bool Product_AX(const Matrix &ma,const Vector &x);
	bool Product_AX(const Matrix &ma);

protected :
	bool SSE2_Enable,AVX_Enable,AVX2_Enable;

	void ProductF_AX(const Matrix &ma,const Vector &x);
	void ProductD_AX(const Matrix &ma,const Vector &x);
};


class Matrix
{
public :
	Matrix(void);
	Matrix(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data);
	Matrix(const Matrix &m);
	virtual ~Matrix(void);

	bool AllocCheck(void) const {return(Coeff!=NULL);}
	bool Create(void);
	bool Create(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data);
	bool Create(const Matrix &m);
	virtual bool CopyStrict(const Matrix &m);
	void Destroy(void);
	bool FillD(const double data);
	bool FillF(const float data);
	bool FillZero(void);
	COEFF_DATA_TYPE GetDataType(void) const {return(data_type);}
	bool SetInfo(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data);
	void GetInfo(uint16_t &l,uint16_t &c,COEFF_DATA_TYPE &data) const;
	uint16_t GetLines(void) const {return(lines);}
	uint16_t GetColumns(void) const {return(columns);}
	void* GetPtrMatrix(void) const {return(Coeff);}
	void* GetPtrMatrixLine(const uint16_t i) const {return((void *)((uint8_t *)Coeff+i*pitch));}
	ptrdiff_t GetPitch(void) const {return(pitch);}
	size_t GetDataSize(void) const {return(size);}
	double GetD(const uint16_t i,const uint16_t j) const {return(((double *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j]);}
	float GetF(const uint16_t i,const uint16_t j) const {return(((float *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j]);}
	void SetD(const uint16_t i,const uint16_t j,const double d) {((double *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j]=d;}
	void SetF(const uint16_t i,const uint16_t j,const float d) {((float *)((uint8_t *)Coeff+(ptrdiff_t)i*pitch))[j]=d;}
	bool GetSafeD(const uint16_t i,const uint16_t j,double &d) const ;
	bool SetSafeD(const uint16_t i,const uint16_t j,const double d);
	bool GetSafeF(const uint16_t i,const uint16_t j,float &d) const ;
	bool SetSafeF(const uint16_t i,const uint16_t j,const float d);

protected :
	void *Coeff;
	uint16_t columns,lines;
	size_t size;
	ptrdiff_t pitch;
	COEFF_DATA_TYPE data_type;

	Matrix& operator=(const Matrix&){return(*this);}
};


class Matrix_Compute : public Matrix
{
public :
	Matrix_Compute(void);
	Matrix_Compute(const uint16_t l,const uint16_t c,const COEFF_DATA_TYPE data);
	Matrix_Compute(const Matrix_Compute &m);
	virtual ~Matrix_Compute(void);

	void EnableSSE2(void) {SSE2_Enable=true;}
	void DisableSSE2(void) {SSE2_Enable=false;}
	void EnableAVX(void) {AVX_Enable=true;}
	void DisableAVX(void) {AVX_Enable=false;}
	void EnableAVX2(void) {AVX2_Enable=true;}
	void DisableAVX2(void) {AVX2_Enable=false;}

	bool CreateTranspose(const Matrix &m);
	virtual bool CopyStrict(const Matrix_Compute &m);
	void SetZeroValue(const double z) {zero_value=fabs(z);}
	double GetZeroValue(void) const {return(zero_value);}
	bool Product_AB(const Matrix &ma,const Matrix &mb);
	bool Product_AtB(const Matrix &ma,const Matrix &mb);
	bool Product_tAA(const Matrix &ma);
	bool Product_tAA(void);
	bool Inverse(const Matrix &ma);
	bool Inverse(void);
	int8_t InverseSafe(const Matrix_Compute &ma);
	int8_t InverseSafe(void);
	bool Transpose(void);
	bool Transpose(const Matrix &ma);

protected :
	double zero_value;
	bool SSE2_Enable,AVX_Enable,AVX2_Enable;

	void ProductD_AB(const Matrix &ma,const Matrix &mb);
	void ProductD_AtB(const Matrix &ma,const Matrix &mb);
	bool InverseD(const Matrix &ma);
	int8_t InverseSafeD(const Matrix_Compute &ma);
	void TransposeD(const Matrix &ma);

	void ProductF_AB(const Matrix &ma,const Matrix &mb);
	void ProductF_AtB(const Matrix &ma,const Matrix &mb);
	bool InverseF(const Matrix &ma);
	int8_t InverseSafeF(const Matrix_Compute &ma);
	void TransposeF(const Matrix &ma);

	void TransposeU64(const Matrix &ma);

	void TransposeI64(const Matrix &ma);

	void TransposeU32(const Matrix &ma);

	void TransposeI32(const Matrix &ma);

	void TransposeU16(const Matrix &ma);

	void TransposeI16(const Matrix &ma);

	void TransposeU8(const Matrix &ma);

	void TransposeI8(const Matrix &ma);

	Matrix_Compute& operator=(const Matrix_Compute&){return(*this);}
};

#endif