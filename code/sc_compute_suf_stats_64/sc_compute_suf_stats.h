// JFA Train Speaker-independent Parameters Define
// Author: xzhang

#ifndef _JFA_FATRAIN_H
#define _JFA_FATRAIN_H

#include "fa_comm.h"
//#define LOGFILE

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

struct Gaussstatis {           // 每个mixture上的统计量 
	float N_occ;
	float *p_F;
	float *p_S;

	Gaussstatis() {
		N_occ=0.0;
		p_F=p_S=NULL;
	}
};

struct GMMStatis				// GMM统计量
{
	float   *Facc;
	float   *Sacc;

	Gaussstatis *pGauss;		// 高斯分量的统计值

	GMMStatis()
	{
		pGauss = NULL;
		Facc=NULL;
		Sacc=NULL;
	}
};

class calc_suf_stats
{
public:
	calc_suf_stats();
	~calc_suf_stats();
	
protected:
	int				m_nMixNum;			// 混合高斯数目
	int				m_nVecSize;			// 特征维数
	int				m_nVecSizeStore;	// 内存空间中用的特征维数（主要是为了兼容IPP版本，4的倍数）
	int				m_nFeatKind;		// 特征类型码
	GaussMixModel	m_UBMGmmModel;		// UBM模型参数


	GMMStatis       m_Statistic;        // statistics for each utterance
	int             m_nUttNum;          // 总的句子数
	float           *m_pfIPPBuf_VecSizeStore;
//	FMatrix         N_utt;              // 每个说话人的total 0-order statistics
//	FMatrix         F_utt;				// 1-order statistics for each speaker
	FVector         N_utt;
	FVector         F_utt;

	char32             *m_cUttNames;
	int                *m_Utt_nFrmNums;


protected:
	bool  AllocGaussModel(GaussMixModel &p_GmmModel);	// 分配GMM模型空间
	void  FreeGaussModel(GaussMixModel  &p_GmmModel);	// 释放GMM模型空间

	bool  AllocBWStatisBuf(GMMStatis &p_Statistic);
	void  FreeBWStatisBuf(GMMStatis &p_Statistic);
	void  ResetBWStatisBuf(GMMStatis &p_Statistic);

	bool  BWStatisticEst(float *p_fFeatBuf, int p_nFrmNum);

public:
	bool  Initalize(int p_nMixNum, int p_nVecSize, int p_nUttNum);
	bool  LoadGMMModel(char *p_pcModelFile);

	bool  compute_stats(char *UttListfile, char *outstats);

	bool  save_stats(char *outfilename);
};

#endif