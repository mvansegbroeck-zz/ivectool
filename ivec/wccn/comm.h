// description: env definition and utilities for fa project 
// author: xzhang

#ifndef _FA_COMM_H
#define _FA_COMM_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ipps.h"
#include "ippsr.h"
#include "ippm.h"

#define Relevance 12
#define ALIGN_4F(n) (((n)+3)&(~3))
#define FILE_HEADER_SIZE 512 
#define MIN_PROB 1.0e-5

#define PI 3.14159265428
#define log2pi    (log(2*PI))
#define LZERO (-1.0E10)

#define IPPUSE 1
#define VceSpaceFlag 1
#define ChlSpaceFlag 0

typedef float* PFLOAT; 
typedef			 char		char32[32];
typedef			 char		char128[128];

#include <map>
#include <string>
using namespace std;

typedef map<string,int> simap;

struct Feature_BaseInfo		// 特征文件头里的保存的信息
{
	char cFeatType[16];		// 例如：MFCC_A
	int  nFrameNum;			// 帧数
	int  nVecSize;			// 特征维数
	int  nVecSizeStore;		// 实际存储的特征维数
	int  nFeatKind;			// 特征类型
	int  nWinSize;			// 短时窗的窗长
	int  nFrameRate;		// 帧率

	int  nTotalParamSize;	// 特征总的参数量，单位：字节

	Feature_BaseInfo()
	{
		cFeatType[0]='\0';
		nFrameNum = nVecSizeStore = nVecSize = nFeatKind = -1;
		nFrameRate = nWinSize = -1;
		nTotalParamSize = 0;
	}
};

struct FMBLOCKMatrix			// 多个Matrix,type is float
{
	float   *pBuf;				// 所有数据放在Buf中
	int	   	nBufSize;			// Buf的大小，即有多少个float

	int     nBlockNum;		 	// 块的数目
	int		nRowDim;			// 每一行的维数
	int		nColDim;			// 每一列的维数
	PFLOAT *pBlockBuf;			// 每个Block的首地址

	FMBLOCKMatrix()
	{
		pBuf = NULL;
		nBufSize = -1;

		nBlockNum = nRowDim = nColDim = -1;
		pBlockBuf = NULL;
	}
};

struct INTMatrix
{
	int	    *pBuf;
	int		nRowDim;			// 行维数
	int		nColDim;			// 列维数

	INTMatrix()
	{
		pBuf=NULL;
		nRowDim = nColDim = -1;
	}
};

struct FMatrix
{
	PFLOAT  pBuf;
	int		nRowDim;			// 行数
	int		nColDim;			// 列数

	FMatrix()
	{
		pBuf=NULL;
		nRowDim = nColDim = -1;
	}
};

struct FVector
{
	PFLOAT  pBuf;
	int		nDim;			// 维数

	FVector()
	{
		pBuf=NULL;
		nDim = -1;
	}
};

inline void Free(void *pBuf);
inline float* Malloc_32f(int p_nSize);
inline void Zero_32f(float *p_pfBuf,int p_nSize);
inline bool AllocVector(FVector &p_Vector,int p_nDim);
inline void FreeVector(FVector &p_Vector);
inline void ZeroVector(FVector &p_Vector);
inline void FreeMatrix(FMatrix &p_Matrix);
inline bool AllocMatrix(FMatrix &p_Matrix, int p_nRowDim, int p_nColDim);
inline void ZeroMatrix(FMatrix &p_Matrix);

void Free(void *pBuf)
{
	if (IPPUSE)
	{
		if (pBuf!=NULL)	ippsFree(pBuf);
	}
	else
	{
		if (pBuf!=NULL)	free(pBuf);
	}
	pBuf=NULL;
}

float* Malloc_32f(int p_nSize)
{
	float *pfBuf;

	if (IPPUSE)
		pfBuf = ippsMalloc_32f(p_nSize);
	else
		pfBuf = (float *)malloc(p_nSize*sizeof(float));

	return pfBuf;
}

void Zero_32f(float *p_pfBuf,int p_nSize)
{
	if (p_pfBuf==NULL || p_nSize <=0 )	return;

	if (IPPUSE)
	{
		ippsZero_32f(p_pfBuf,p_nSize);
	}
	else
	{
		memset(p_pfBuf,sizeof(float),p_nSize);
	}
}

void FreeVector(FVector &p_Vector)
{
	if (p_Vector.pBuf!=NULL)
		Free(p_Vector.pBuf);
	p_Vector.pBuf = NULL;		

	p_Vector.nDim = -1;
}

bool AllocVector(FVector &p_Vector,int p_nDim)
{
	FreeVector(p_Vector);

	p_Vector.pBuf = Malloc_32f(p_nDim);
	if (p_Vector.pBuf==NULL)		return false;

	p_Vector.nDim = p_nDim;

	return true;
}

void ZeroVector(FVector &p_Vector)
{
	if (p_Vector.pBuf!=NULL)
	{
		Zero_32f(p_Vector.pBuf,p_Vector.nDim);
	}
}

void FreeMatrix(FMatrix &p_Matrix)
{
	if (p_Matrix.pBuf!=NULL)
		Free(p_Matrix.pBuf);
	p_Matrix.pBuf = NULL;		

	p_Matrix.nColDim = p_Matrix.nRowDim = -1;
}

bool AllocMatrix(FMatrix &p_Matrix, int p_nRowDim, int p_nColDim)
{
//	FreeMatrix(p_Matrix);

	p_Matrix.pBuf = Malloc_32f(p_nRowDim*p_nColDim);
	if (p_Matrix.pBuf==NULL)		return false;

	p_Matrix.nRowDim = p_nRowDim;
	p_Matrix.nColDim = p_nColDim;

	return true;
}

void ZeroMatrix(FMatrix &p_Matrix)
{
	if (p_Matrix.pBuf!=NULL)
	{
		Zero_32f(p_Matrix.pBuf,p_Matrix.nColDim*p_Matrix.nRowDim);
	}
}

inline float ComputeDotMulti(float *Vec1,float *Vec2,int Veclen)
{
	float totDotvalue=0.0;
	for(int i=0;i<Veclen;i++)
	{
		totDotvalue+=Vec1[i]*Vec2[i];
	}
	return(totDotvalue);
}

inline float ComputeDotMulti_IPP(float *Vec1, float *Vec2, int Veclen)
{
	Ipp32f pDst;
	IppStatus status = ippmDotProduct_vv_32f(Vec1, sizeof(float), Vec2, sizeof(float), &pDst, Veclen);
	if (status!=0) {
		printf("Error in ComputeDotMulti_IPP()!\n");
		exit(-1);
	}
	return(pDst);
}

inline bool matrixTranspos(FMatrix matrix, FMatrix &matrix_tans)
{
	if(!AllocMatrix(matrix_tans, matrix.nColDim, matrix.nRowDim))
		return false;

	for(int i=0; i<matrix.nRowDim; i++)
		for(int j=0; j<matrix.nColDim; j++)
			matrix_tans.pBuf[j*matrix.nRowDim+i]=matrix.pBuf[i*matrix.nColDim+j];

	return true;
}

inline bool matrix_invert(FMatrix matrix, FMatrix &matrix_inv)
{
	if (matrix.pBuf==NULL) 
		return false;
	int widthHeight = matrix.nRowDim;
	int srcStride2 = sizeof(float);
	int srcStride1 = widthHeight * sizeof(float);

	if(!AllocMatrix(matrix_inv, widthHeight, widthHeight))
		return false;
	int dstStride2 = sizeof(float);
	int dstStride1 = widthHeight * sizeof(float);

	Ipp32f *pBuffer;
	pBuffer = Malloc_32f(widthHeight*widthHeight+widthHeight);
	if (pBuffer==NULL) 
	 return false;

	IppStatus status = ippmInvert_m_32f((const Ipp32f*)matrix.pBuf, srcStride1, srcStride2, pBuffer, \
		(Ipp32f *)matrix_inv.pBuf, dstStride1, dstStride2, widthHeight);

	if (status==0) {
		printf("Invert zeta ok!\n");
	}
	else {
		printf("Invert zeta error!\n");
		return false;
	}
	Free(pBuffer);
	return true;
}

inline bool matrixAdd(FMatrix &matrix1, FMatrix matrix2)
{
	if (matrix1.nRowDim!=matrix2.nRowDim||matrix1.nColDim!=matrix2.nColDim) {
		return false;
	}
	if (matrix1.pBuf==NULL || matrix2.pBuf==NULL) {
		return false;
	}

	for (int i=0; i<matrix1.nRowDim; ++i) {
		for (int j=0; j<matrix1.nColDim; ++j) {
			matrix1.pBuf[i*matrix1.nColDim+j] += matrix2.pBuf[i*matrix1.nColDim+j];
		}
	}
	return true;
}

// multi内部开内存,外部释放
inline bool matrixMulti(FMatrix matrixA, FMatrix matrixB, FMatrix &multi)
{
	if (matrixA.nColDim != matrixB.nRowDim) {
		printf("Please check the width and height of the two matrixs.\n");
		return(false);
	}
	
	FMatrix transmatrixB;
	if (!matrixTranspos(matrixB, transmatrixB)) 
		return false;
//	if (!AllocMatrix(transmatrixB, matrixB.nColDim, matrixB.nRowDim)) 
//		return false;

//	ippmTranspose_m_32f((const Ipp32f*)matrixB.pBuf, matrixB.nColDim*sizeof(float), sizeof(float), matrixB.nColDim, \
//		matrixB.nRowDim, transmatrixB.pBuf, transmatrixB.nColDim*sizeof(float), sizeof(float));

	if (!AllocMatrix(multi, matrixA.nRowDim, matrixB.nColDim)) {
		return false;
	}
	
	for(int i=0;i<matrixA.nRowDim;i++)
		for(int j=0; j<matrixB.nColDim; j++) {
			multi.pBuf[i*matrixB.nColDim+j] = ComputeDotMulti_IPP(&matrixA.pBuf[i*matrixA.nColDim], \
				&transmatrixB.pBuf[j*matrixB.nRowDim], matrixA.nColDim);
		}
	
	FreeMatrix(transmatrixB);
	return(true);
}

inline bool matrixMulti_IPP(FMatrix matrixA, FMatrix matrixB, FMatrix &multi)
{
	float *pSrc1 = matrixA.pBuf;
	float *pSrc2 = matrixB.pBuf;
	int src1Width = matrixA.nColDim;
	int src1Height = matrixA.nRowDim;
	int src1Stride2 = sizeof(float);
	int src1Stride1 = src1Width * sizeof(float);

	int src2Width = matrixB.nColDim;
	int src2Height = matrixB.nRowDim;
	int src2Stride2 = sizeof(float);
	int src2Stride1 = src2Width * sizeof(float);

	if (!AllocMatrix(multi, src1Height, src2Width)) {
		return false;
	}

	float *pDst = multi.pBuf;
	int dstStride2 = sizeof(float);
	int dstStride1 = multi.nColDim * sizeof(float);
	
	IppStatus status = ippmMul_mm_32f(pSrc1, src1Stride1, src1Stride2, src1Width, src1Height, \
		pSrc2, src2Stride1, src2Stride2, src2Width, src2Height, pDst, dstStride1, dstStride2);

	if (status!=0) {
		return false;
	}

	return true;
}

inline bool matrixMulti_IPP_without_malloc(FMatrix matrixA, FMatrix matrixB, FMatrix &multi)
{
	float *pSrc1 = matrixA.pBuf;
	float *pSrc2 = matrixB.pBuf;
	int src1Width = matrixA.nColDim;
	int src1Height = matrixA.nRowDim;
	int src1Stride2 = sizeof(float);
	int src1Stride1 = src1Width * sizeof(float);

	int src2Width = matrixB.nColDim;
	int src2Height = matrixB.nRowDim;
	int src2Stride2 = sizeof(float);
	int src2Stride1 = src2Width * sizeof(float);

	if (multi.nRowDim!=src1Height || multi.nColDim!= src2Width) {
		return false;
	}

//	if (!AllocMatrix(multi, src1Height, src2Width)) {
//		return false;
//	}

	float *pDst = multi.pBuf;
	int dstStride2 = sizeof(float);
	int dstStride1 = multi.nColDim * sizeof(float);
	
	IppStatus status = ippmMul_mm_32f(pSrc1, src1Stride1, src1Stride2, src1Width, src1Height, \
		pSrc2, src2Stride1, src2Stride2, src2Width, src2Height, pDst, dstStride1, dstStride2);

	if (status!=0) {
		return false;
	}

	return true;
}

inline bool matrixMulti_Btrans(FMatrix matrixA, FMatrix matrixB, FMatrix &multi)
{
	//A*B^T
	if (matrixA.nColDim != matrixB.nColDim) {
		printf("Please check the width and height of the two matrixs.\n");
		return(false);
	}
	
//	FMatrix transmatrixB;
//	if (!matrixTranspos(matrixB, transmatrixB)) 
//		return false;

	if (!AllocMatrix(multi, matrixA.nRowDim, matrixB.nRowDim)) {
		return false;
	}
	
	for(int i=0;i<matrixA.nRowDim;i++)
		for(int j=0; j<matrixB.nRowDim; j++)
			multi.pBuf[i*matrixB.nRowDim+j] = ComputeDotMulti_IPP(&matrixA.pBuf[i*matrixA.nColDim], \
				&matrixB.pBuf[j*matrixB.nColDim], matrixB.nColDim);
	
//	FreeMatrix(transmatrixB);
	return(true);
}

//-------------------- float型的实对称矩阵运算 逆-----------------------------------
inline int mathRealMatrixSymmetryInversion (float *aa, int n)
{
	int		*is, *js, i,j,k,l,u,v;
	float	d,p;
	float	*a;
	
	a = new float[n*n];
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			a[i*n+j]=aa[i*n+j];
		
		is=new int[n];
		js=new int[n];
		
		for(k=0;k<n;k++)
		{
			d=0.0;
			for(i=k;i<n;i++)
				for(j=k;j<n;j++)
				{
					l=i*n+j;		p=(float)fabs(a[l]);
					
					if(p>d)
					{	d=p;	is[k]=i;	js[k]=j;			}
				}
				
				if( fabs(d) <= 1.0e-12)
				{
					delete is;		delete js;
					printf("error in dcinv!\n");			
					return false;		
				}
				
				if(is[k]!=k)
					for(j=0;j<n;j++)
					{
						u=k*n+j;	v=is[k]*n+j;		p=a[u];		a[u] = a[v];	a[v] = p;
					}
					
					if(js[k]!=k)
						for(i=0;i<n;i++)
						{
							u=i*n+k;	v=i*n+js[k];	p=a[u];		a[u]=a[v];		a[v]=p;
						}
						
						l=k*n+k;		a[l]=1.0f/a[l];
						for(j=0;j<n;j++)
							if(j!=k)
							{		u=k*n+j;		a[u]*=a[l];			}
							
							for(i=0;i<n;i++)
								if(i!=k)
									for(j=0;j<n;j++)
										if(j!=k)
										{		u=i*n+j;		a[u]-=a[i*n+k]*a[k*n+j];	}
										
										for(i=0;i<n;i++)
											if(i!=k)
											{		u=i*n+k;		a[u]*=-a[l];		}
		}// end of for(k...)
		
		
		for(k=n-1;k>=0;k--)
		{
			if(js[k]!=k)
				for(j=0;j<n;j++)
				{
					u=k*n+j;	v=js[k]*n+j;		p=a[u];		a[u]=a[v];		a[v]=p;
				}
				
				if(is[k]!=k)
					for(i=0;i<n;i++)
					{
						u=i*n+k;	v=i*n+is[k];		p=a[u];		a[u]=a[v];		a[v]=p;
					}
		}
		
		for(i=0;i<n;i++)
			for(j=0;j<n;j++)
				aa[i*n+j]=a[i*n+j];	
			
			delete is; delete js; delete a;
			return true;
			
}

inline FILE* fa_fopen(const char* file, const char* mode)
{
	FILE* fp = fopen(file, mode);
	if (NULL == fp)
	{
		printf("[ERROR]: cannot open file \'%s\' as \'%s\'\n", file ,mode);
		exit(1);
	}
	return fp;
}


inline bool ReadFeatFile(char *p_pcFeatFile,float *&p_pfFeatBuf,Feature_BaseInfo &p_sFeatInfo)
{
	if (p_pcFeatFile==NULL)	return false;

	////////////////////////////////////////////
	FILE *fp_feat;
	fp_feat = fopen(p_pcFeatFile,"rb");
	if (fp_feat==NULL)
	{
		printf("Error open %s for read!\n",p_pcFeatFile);
		return false;
	}

	// 读文件头
	fread(&p_sFeatInfo,sizeof(Feature_BaseInfo),1,fp_feat);
	if (p_sFeatInfo.nTotalParamSize<=0 || p_sFeatInfo.nTotalParamSize!=p_sFeatInfo.nFrameNum*p_sFeatInfo.nVecSizeStore*sizeof(float))
	{
		printf("Something is error in %s\n",p_pcFeatFile);
		fclose(fp_feat);
		return false;
	}

	fseek(fp_feat,FILE_HEADER_SIZE,SEEK_SET);

	if (p_pfFeatBuf!=NULL) {
		free(p_pfFeatBuf);
		p_pfFeatBuf=NULL;
	}
	p_pfFeatBuf = (float *)malloc(p_sFeatInfo.nTotalParamSize);
	if (p_pfFeatBuf==NULL)
	{
		printf("Error malloc(%d) in ReadFeatFile()\n",p_sFeatInfo.nTotalParamSize);
		fclose(fp_feat);
		return false;
	}
	int nNumRead = fread(p_pfFeatBuf,1,p_sFeatInfo.nTotalParamSize,fp_feat);
	if (nNumRead!=p_sFeatInfo.nTotalParamSize)
	{
		printf("Error in %s. Can't read enough param.\n",p_pcFeatFile);
		fclose(fp_feat);
		free(p_pfFeatBuf);
		p_pfFeatBuf = NULL;
		return false;
	}

	fclose(fp_feat);
	return true;	
}

inline const char* fa_path2file(const char* path)
{
	const char* p = path+strlen(path)-1;
	while (p>=path && *p!='\\' && *p!='/') --p;
	return p+1;

}
#endif
