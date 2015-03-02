// JFA Train Speaker-independent Parameters Define
// Author: xzhang

#ifndef _JFA_FATRAIN_H
#define _JFA_FATRAIN_H

#include "fa_comm.h"
//#define LOGFILE

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

struct Gaussstatis {           // ÿ��mixture�ϵ�ͳ���� 
	float N_occ;
	float *p_F;
	float *p_S;

	Gaussstatis() {
		N_occ=0.0;
		p_F=p_S=NULL;
	}
};

struct GMMStatis				// GMMͳ����
{
	float   *Facc;
	float   *Sacc;

	Gaussstatis *pGauss;		// ��˹������ͳ��ֵ

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
	int				m_nMixNum;			// ��ϸ�˹��Ŀ
	int				m_nVecSize;			// ����ά��
	int				m_nVecSizeStore;	// �ڴ�ռ����õ�����ά������Ҫ��Ϊ�˼���IPP�汾��4�ı�����
	int				m_nFeatKind;		// ����������
	GaussMixModel	m_UBMGmmModel;		// UBMģ�Ͳ���


	GMMStatis       m_Statistic;        // statistics for each utterance
	int             m_nUttNum;          // �ܵľ�����
	float           *m_pfIPPBuf_VecSizeStore;
//	FMatrix         N_utt;              // ÿ��˵���˵�total 0-order statistics
//	FMatrix         F_utt;				// 1-order statistics for each speaker
	FVector         N_utt;
	FVector         F_utt;

	char32             *m_cUttNames;
	int                *m_Utt_nFrmNums;


protected:
	bool  AllocGaussModel(GaussMixModel &p_GmmModel);	// ����GMMģ�Ϳռ�
	void  FreeGaussModel(GaussMixModel  &p_GmmModel);	// �ͷ�GMMģ�Ϳռ�

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