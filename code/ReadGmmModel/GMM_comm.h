/******************************DOCUMENT*COMMENT***********************************
*D
*D 文件名称            : GMM_comm.h
*D
*D 项目名称            : 
*D
*D 版本号              : 1.1.0001
*D
*D 文件描述            :
*D
*D
*D 文件修改记录
*D ------------------------------------------------------------------------------ 
*D 版本号       修改日期       修改人     改动内容
*D ------------------------------------------------------------------------------ 
*D 1.1.0001     2007.12.21     plu        创建文件
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
	int   nMixNum;			// 混合高斯数
	int	  nVecSize;			// 特征维数

	float fOccFloor;		// 高斯参数被更新的occ门限
	float fMeanFactor;		// 均值自适应因子
	float fIterFactor;		// 控制map次数的因子，[0,1]
	
	int   nTopMixNum;		// FastScore用

	bool  bDoZNorm;			// 是否做ZNorm
	
	bool  bDoATNorm;		// 是否做ATNorm
	int   nCohortMdlNum;	// Cohort集中说话人模型的数目
	int	  nATNormNum;		// ATNorm的数目

	bool  bDoLFA;			// 是否做LFA
	int   nLFA_RankNum;		// Rank数

	bool  bIPPUse;			// 是否使用IPP

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
struct GMMFileHeader			// 混合高斯模型文件的文件头
{
	int	  nMixNum;				// 高斯分量的数目
	int   nVecSize;				// 特征维数维数
	int   nVecSizeStore;		// 文件中实际存储的维数
	int   nFeatKind;			// 特征类型码

	int   nTotalParamSize;		// 文件中总的参数的字节数

	GMMFileHeader()
	{
		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;
		
		nTotalParamSize = 0;
	}
};

struct GaussMixModel			// 混合高斯模型
{
    PFLOAT	pfWeightBuf;		// 权重矢量 
	PFLOAT	pfMeanBuf;			// 均值矢量
	PFLOAT	pfDiagCovBuf;		// 对角逆协方差 
	PFLOAT	pfMatBuf;			// Mat
	
	GaussMixModel()
	{
		pfWeightBuf = pfMeanBuf = pfDiagCovBuf = pfMatBuf =NULL;
	}
};

struct SpkMdlFileHeader			// 说话人模型文件的头文件。!!! 注意：目前版本只自适应了均值矢量
{
	int   nMixNum;				// 高斯分量的数目
	int   nVecSize;				// 特征维数维数
	int   nVecSizeStore;		// 文件中实际存储的维数
	int   nFeatKind;			// 特征类型码

	bool  bDoZNorm;				// 是否ZNorm
	bool  bDoATNorm;			// 是否ATNorm

	int   nATNormSelectNum;		// ATNorm中被选中的CohortModel的数目

	int   nMeanParamSize;		// 文件中均值矢量的大小，单位字节
	int   nZNormParamSize;		// 文件中ZNorm参数的大小，单位字节
	int   nATNormParamSize;		// 文件中ATNorm参数的大小，单位字节

	int   nTotalParamSize;		// 文件中总参数的大小，单位字节

	SpkMdlFileHeader()
	{
		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;
		bDoATNorm = bDoZNorm = false;

		nATNormSelectNum = 0;
		nMeanParamSize = nZNormParamSize = nATNormParamSize = nTotalParamSize = 0;
	}
};

struct GMMSpeakerModel			// 说话人模型参数
{
	PFLOAT  pfSpkMeanBuf;	    // 均值参数
	int    *pnCohortMdlIndex;	// Cohort模型的index，可为NULL
	PFLOAT  pfZNorm_Mean;		// ZNorm对应的均值，可为NULL
	PFLOAT  pfZNorm_StdDev;		// ZNorm对应的逆的标准偏差，可为NULL

	GMMSpeakerModel()
	{
		pfSpkMeanBuf=NULL;
		pnCohortMdlIndex = NULL;
		pfZNorm_Mean = pfZNorm_StdDev = NULL;
	}
};

struct CohortModelFileHeader	// TNorm或ATNorm用的Cohort模型文件的文件头。注意：目前Cohort模型只有均值矢量
{
	int   nCohortMdlNum;
	
	int   nMixNum;				// 高斯分量的数目
	int   nVecSize;				// 特征维数维数
	int   nVecSizeStore;		// 文件中实际存储的维数
	int   nFeatKind;			// 特征类型码

	bool  bDoZNorm;				// 是否做ZNorm
		
	int   nTotalParamSize;		// 文件中所有参数的大小：字节

	CohortModelFileHeader()
	{
		nCohortMdlNum = -1;

		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;

		bDoZNorm = false;

		nTotalParamSize = 0;
	}
};

struct LFASpaceFileHeader		// LFA的V空间文件头
{
	int   nSuperDim;			// 超矢量的维数
	int   nRankNum;				// LFA 子空间的Rank数目
	int   nMixNum;				// 混合高斯数目
	int   nDim;					// 实际特征维数

	int   nTotalParamSize;		// 文件中所有参数的大小：字节

	LFASpaceFileHeader()
	{
		nRankNum = nSuperDim = nDim = nMixNum = -1;
		nTotalParamSize = 0;
	}
};

struct SuperVectorFileHeader		// SuperVector文件头 //xxiao change -3.13
{
	int   nSuperVecNum;			// 超矢量的个数
	int   nSuperVecDim;				// 超矢量的总维数
	short   nAddFloat;					// 头上是否多存有1.0, -1_No, +1_Yes
	short   nUseless;		// 文件中所有参数的大小：字节

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
struct GaussStatis				// 高斯统计量（训练时用）
{
	PFLOAT pfAcc;				// 一阶统计量
	float  fOcc;				// Pr(i|X) 观测特征属于第i个高斯分量的后验概率

	GaussStatis()
	{
		pfAcc=NULL;
		fOcc=0.0;
	}
};

struct GMMStatistic				// GMM统计量
{
	PFLOAT      pfBuf;			// 总的buf 
	GaussStatis *pGauss;		// 高斯分量的统计值

	GMMStatistic()
	{
		pGauss = NULL;
		pfBuf = NULL;
	}
};

//////////////////////////////////////////////////////////////////////////
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
	int		nRowDim;			// 每一行的维数
	int		nColDim;			// 每一列的维数

	INTMatrix()
	{
		pBuf=NULL;
		nRowDim = nColDim = -1;
	}
};

struct FMatrix
{
	PFLOAT  pBuf;
	int		nRowDim;			// 每一行的维数
	int		nColDim;			// 每一列的维数

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

//////////////////////////////////////////////////////////////////////////
struct COHORTSCORE_POOL					// TNorm用
{
	bool  bIsCalc;						// true表示已经计算过了
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