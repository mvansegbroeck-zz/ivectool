/******************************DOCUMENT*COMMENT***********************************
*D
*D �ļ�����            : GMM_comm.h
*D
*D ��Ŀ����            : 
*D
*D �汾��              : 1.1.0001
*D
*D �ļ�����            :
*D
*D
*D �ļ��޸ļ�¼
*D ------------------------------------------------------------------------------ 
*D �汾��       �޸�����       �޸���     �Ķ�����
*D ------------------------------------------------------------------------------ 
*D 1.1.0001     2007.12.21     plu        �����ļ�
*D*******************************************************************************/
#ifndef _GMM_COMM_20071221_H_
#define _GMM_COMM_20071221_H_

#include <malloc.h>
#include <stdio.h>
#include <memory.h>

#define MAX_TOP_MIX_NUM		128		
#define MAX_ATNORM_NUM		1000	

#define ZERO				0.0000
#define MinuInf				-100000

typedef float* PFLOAT; 
typedef int    Frame_Gauss_Index[MAX_TOP_MIX_NUM];
typedef float  Frame_Gauss_Score[MAX_TOP_MIX_NUM];

//////////////////////////////////////////////////////////////////////////
struct GMMMapScore_BaseInfo
{
	int   nMixNum;			// ��ϸ�˹��
	int	  nVecSize;			// ����ά��

	float fOccFloor;		// ��˹���������µ�occ����
	float fMeanFactor;		// ��ֵ����Ӧ����
	float fIterFactor;		// ����map���������ӣ�[0,1]
	
	int   nTopMixNum;		// FastScore��

	bool  bDoZNorm;			// �Ƿ���ZNorm
	
	bool  bDoATNorm;		// �Ƿ���ATNorm
	int   nCohortMdlNum;	// Cohort����˵����ģ�͵���Ŀ
	int	  nATNormNum;		// ATNorm����Ŀ

	bool  bDoLFA;			// �Ƿ���LFA
	int   nLFA_RankNum;		// Rank��

	bool  bIPPUse;			// �Ƿ�ʹ��IPP

	GMMMapScore_BaseInfo()
	{
		nMixNum = nVecSize = -1;
		fOccFloor = fMeanFactor = fIterFactor = -1.f;

		nTopMixNum = -1;

		bDoZNorm = false;
		bDoATNorm = false;
		nCohortMdlNum = nATNormNum = -1;

		bDoLFA = false;
		nLFA_RankNum = -1;
		
		bIPPUse = false;
	}
};

///////////////////////////////////////////////////////////////////////////////
struct GMMFileHeader			// ��ϸ�˹ģ���ļ����ļ�ͷ
{
	int	  nMixNum;				// ��˹��������Ŀ
	int   nVecSize;				// ����ά��ά��
	int   nVecSizeStore;		// �ļ���ʵ�ʴ洢��ά��
	int   nFeatKind;			// ����������

	int   nTotalParamSize;		// �ļ����ܵĲ������ֽ���

	GMMFileHeader()
	{
		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;
		
		nTotalParamSize = 0;
	}
};

struct GaussMixModel			// ��ϸ�˹ģ��
{
    PFLOAT	pfWeightBuf;		// Ȩ��ʸ�� 
	PFLOAT	pfMeanBuf;			// ��ֵʸ��
	PFLOAT	pfDiagCovBuf;		// �Խ���Э���� 
	PFLOAT	pfMatBuf;			// Mat
	
	GaussMixModel()
	{
		pfWeightBuf = pfMeanBuf = pfDiagCovBuf = pfMatBuf =NULL;
	}
};

struct SpkMdlFileHeader			// ˵����ģ���ļ���ͷ�ļ���!!! ע�⣺Ŀǰ�汾ֻ����Ӧ�˾�ֵʸ��
{
	int   nMixNum;				// ��˹��������Ŀ
	int   nVecSize;				// ����ά��ά��
	int   nVecSizeStore;		// �ļ���ʵ�ʴ洢��ά��
	int   nFeatKind;			// ����������

	bool  bDoZNorm;				// �Ƿ�ZNorm
	bool  bDoATNorm;			// �Ƿ�ATNorm

	int   nATNormSelectNum;		// ATNorm�б�ѡ�е�CohortModel����Ŀ

	int   nMeanParamSize;		// �ļ��о�ֵʸ���Ĵ�С����λ�ֽ�
	int   nZNormParamSize;		// �ļ���ZNorm�����Ĵ�С����λ�ֽ�
	int   nATNormParamSize;		// �ļ���ATNorm�����Ĵ�С����λ�ֽ�

	int   nTotalParamSize;		// �ļ����ܲ����Ĵ�С����λ�ֽ�

	SpkMdlFileHeader()
	{
		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;
		bDoATNorm = bDoZNorm = false;

		nATNormSelectNum = 0;
		nMeanParamSize = nZNormParamSize = nATNormParamSize = nTotalParamSize = 0;
	}
};

struct GMMSpeakerModel			// ˵����ģ�Ͳ���
{
	PFLOAT  pfSpkMeanBuf;	    // ��ֵ����
	int    *pnCohortMdlIndex;	// Cohortģ�͵�index����ΪNULL
	PFLOAT  pfZNorm_Mean;		// ZNorm��Ӧ�ľ�ֵ����ΪNULL
	PFLOAT  pfZNorm_StdDev;		// ZNorm��Ӧ����ı�׼ƫ���ΪNULL

	GMMSpeakerModel()
	{
		pfSpkMeanBuf=NULL;
		pnCohortMdlIndex = NULL;
		pfZNorm_Mean = pfZNorm_StdDev = NULL;
	}
};

struct CohortModelFileHeader	// TNorm��ATNorm�õ�Cohortģ���ļ����ļ�ͷ��ע�⣺ĿǰCohortģ��ֻ�о�ֵʸ��
{
	int   nCohortMdlNum;
	
	int   nMixNum;				// ��˹��������Ŀ
	int   nVecSize;				// ����ά��ά��
	int   nVecSizeStore;		// �ļ���ʵ�ʴ洢��ά��
	int   nFeatKind;			// ����������

	bool  bDoZNorm;				// �Ƿ���ZNorm
		
	int   nTotalParamSize;		// �ļ������в����Ĵ�С���ֽ�

	CohortModelFileHeader()
	{
		nCohortMdlNum = -1;

		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;

		bDoZNorm = false;

		nTotalParamSize = 0;
	}
};

struct LFASpaceFileHeader		// LFA��V�ռ��ļ�ͷ
{
	int   nSuperDim;			// ��ʸ����ά��
	int   nRankNum;				// LFA �ӿռ��Rank��Ŀ
	int   nMixNum;				// ��ϸ�˹��Ŀ
	int   nDim;					// ʵ������ά��

	int   nTotalParamSize;		// �ļ������в����Ĵ�С���ֽ�

	LFASpaceFileHeader()
	{
		nRankNum = nSuperDim = nDim = nMixNum = -1;
		nTotalParamSize = 0;
	}
};

struct SuperVectorFileHeader		// SuperVector�ļ�ͷ //xxiao change -3.13
{
	int   nSuperVecNum;			// ��ʸ���ĸ���
	int   nSuperVecDim;				// ��ʸ������ά��
	short   nAddFloat;					// ͷ���Ƿ�����1.0, -1_No, +1_Yes
	short   nUseless;		// �ļ������в����Ĵ�С���ֽ�

	SuperVectorFileHeader()
	{
		nSuperVecNum = nSuperVecDim = -1;
		nAddFloat = 0;
	}
};

struct NGramBinFileHeader
{
	int nSmpNum;
	int nVecDim;
	short nUseless1;
	short nUseless2;
	NGramBinFileHeader()
	{
		nSmpNum = nVecDim = -1;
	}
};


//////////////////////////////////////////////////////////////////////////
struct GaussStatis				// ��˹ͳ������ѵ��ʱ�ã�
{
	PFLOAT pfAcc;				// һ��ͳ����
	float  fOcc;				// Pr(i|X) �۲��������ڵ�i����˹�����ĺ������

	GaussStatis()
	{
		pfAcc=NULL;
		fOcc=0.0;
	}
};

struct GMMStatistic				// GMMͳ����
{
	PFLOAT      pfBuf;			// �ܵ�buf 
	GaussStatis *pGauss;		// ��˹������ͳ��ֵ

	GMMStatistic()
	{
		pGauss = NULL;
		pfBuf = NULL;
	}
};

//////////////////////////////////////////////////////////////////////////
struct FMBLOCKMatrix			// ���Matrix,type is float
{
	float   *pBuf;				// �������ݷ���Buf��
	int	   	nBufSize;			// Buf�Ĵ�С�����ж��ٸ�float

	int     nBlockNum;		 	// �����Ŀ
	int		nRowDim;			// ÿһ�е�ά��
	int		nColDim;			// ÿһ�е�ά��
	PFLOAT *pBlockBuf;			// ÿ��Block���׵�ַ

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
	int		nRowDim;			// ÿһ�е�ά��
	int		nColDim;			// ÿһ�е�ά��

	INTMatrix()
	{
		pBuf=NULL;
		nRowDim = nColDim = -1;
	}
};

struct FMatrix
{
	PFLOAT  pBuf;
	int		nRowDim;			// ÿһ�е�ά��
	int		nColDim;			// ÿһ�е�ά��

	FMatrix()
	{
		pBuf=NULL;
		nRowDim = nColDim = -1;
	}
};

struct FVector
{
	PFLOAT  pBuf;
	int		nDim;			// ά��

	FVector()
	{
		pBuf=NULL;
		nDim = -1;
	}
};

//////////////////////////////////////////////////////////////////////////
struct COHORTSCORE_POOL					// TNorm��
{
	bool  bIsCalc;						// true��ʾ�Ѿ��������
	float fScore;	

	COHORTSCORE_POOL()
	{
		bIsCalc = false;
		fScore  = MinuInf;
	}
};
//////////////////////////////////////////////////////////////////////////
struct NGramModel 
{
	int nRowDim;	// smpNum
	int nColDim;	//Dim of each smp
	FMatrix p_fProbArray;
	INTMatrix p_nIndexArray;
	NGramModel()
	{
		p_fProbArray.nColDim=nColDim;
		p_fProbArray.nRowDim=nRowDim;
		p_nIndexArray.nColDim=nColDim;
		p_nIndexArray.nRowDim=nRowDim;
		p_fProbArray.pBuf=NULL;
		p_nIndexArray.pBuf=NULL;
	}
};

#endif // _GMM_COMM_H_