#include "comm_srlr.h"
#include "memory_srlr.h"


struct GMMFileHeaderOld			// 混合高斯模型文件的文件头
{
	int nModelNum;				// 混合高斯模型的数目
	int nMixNum;				// 混合高斯的数目
	int nDim;					// 维数
	int nMfccKind;				// 2007.09.11 plu : 特征类型

	GMMFileHeaderOld()
	{
		nModelNum = nMixNum = nDim = nMfccKind = -1 ; // 表示无效
	}
};

struct GaussPDF					// 高斯概率密度
{	
	float	*pfMean;			// 均值矢量，不分配内存，只是保存每个高斯分量均值矢量的首地址
	float	*pfDiagCov;			// 对角逆斜方差阵（矢量），内存同上
	double	dMat;				// = log(weight) - 0.5*Dim*log(2pi) + 0.5*log|逆斜方差阵|
	int		nDim;				// 有效维数（不一定是4的倍数）

	GaussPDF()
	{
		pfMean=pfDiagCov=NULL;
		dMat = 0.0;
		nDim = 0;
	}
};

struct GaussMixModel			// 混合高斯模型
{
    float		*pfWeight;		// 权重矢量 
	float		*pfMeanBuf;		// 2007.09.17 plu : 所有的均值
	float		*pfDiagCovBuf;	// 2007.09.17 plu : 所有的方差
    GaussPDF	*pGauss;		// Gaussian component
	int			nMixNum;		// 混合高斯的数目

	GaussMixModel()
	{
		pfWeight= NULL;
		pGauss  = NULL;
		nMixNum = 0;

		pfMeanBuf = pfDiagCovBuf = NULL;	// 2007.09.17 plu : 
	}
};



struct GMMFileHeaderNew			// 混合高斯模型文件的文件头
{
	int	  nMixNum;				// 高斯分量的数目
	int   nVecSize;				// 特征维数维数
	int   nVecSizeStore;		// 文件中实际存储的维数
	int   nFeatKind;			// 特征类型码

	int   nTotalParamSize;		// 文件中总的参数的字节数

	GMMFileHeaderNew()
	{
		nMixNum = nVecSizeStore = nVecSize = nFeatKind = -1;
		
		nTotalParamSize = 0;
	}
};
void AllocGaussMixModel(GaussMixModel *p_pModel,int p_nMixNum,int p_nVecSize);
void FreeGaussMixModel(GaussMixModel *p_gmmParam,int p_nModelNum);

int main(int argc,char * argv[])
{

	if (argc!=3)
	{
		printf("OldUbm2NewUbm.exe oldubm newubm\n");
		exit(0);
	}
	
	GaussMixModel * m_pGmmModel;
	GMMFileHeaderOld  ModelHeaderOld;					// Header *header; 模型文件的文件头
	GMMFileHeaderNew  ModelHeaderNew;
	FILE *fpModel;
	ReadOpen(fpModel,argv[1]);

	fread(&ModelHeaderOld,sizeof(GMMFileHeaderOld),1,fpModel);
	ASSERT3(ModelHeaderOld.nDim>0, "Error in model file %s : nDim<0!",argv[1]);
	ASSERT3(ModelHeaderOld.nMixNum>0, "Error in model file %s : nMixNum<0!",argv[1]);
	ASSERT3(ModelHeaderOld.nModelNum>0, "Error in model file %s : nModelNum<0!",argv[1]);

	int m_nVecSize = ModelHeaderOld.nDim;
	int m_nVecSize4 = ALIGN_4F(m_nVecSize);
	int m_nMixNum  = ModelHeaderOld.nMixNum;		
	int m_nModelNum  = ModelHeaderOld.nModelNum;

	ASSERT3(m_nModelNum==1,"m_nModelNum=%d\n",m_nModelNum);

	m_pGmmModel = (GaussMixModel *)Malloc(m_nModelNum,sizeof(GaussMixModel),false);

	for(int i=0;i<m_nModelNum;i++)
	{
		// 分配 GMM 模型
		AllocGaussMixModel(&m_pGmmModel[i],m_nMixNum,m_nVecSize4);

		m_pGmmModel[i].nMixNum = m_nMixNum;

		// 读出weight
		fread(m_pGmmModel[i].pfWeight,sizeof(float),m_nMixNum,fpModel);

		// 读出均值，斜方差，mat
		for(int m=0;m<m_nMixNum;m++)
		{
			m_pGmmModel[i].pGauss[m].nDim = m_nVecSize;

			fread(m_pGmmModel[i].pGauss[m].pfMean,sizeof(float),m_nVecSize4,fpModel);
			fread(m_pGmmModel[i].pGauss[m].pfDiagCov,sizeof(float),m_nVecSize4,fpModel);
			fread(&m_pGmmModel[i].pGauss[m].dMat,sizeof(double),1,fpModel);
		
			// 如果mat值大于零，可能模型有问题，报警
			if (m_pGmmModel[i].pGauss[m].dMat>=0.0)
				printf("Warning : m_pGmmModel[%d].pGauss[%d].dMat=%.3f!\n",i,m,m_pGmmModel[i].pGauss[m].dMat);
		}
	}
	fclose(fpModel);


	ModelHeaderNew.nMixNum=m_nMixNum;
	ModelHeaderNew.nVecSize=m_nVecSize;
	ModelHeaderNew.nVecSizeStore=m_nVecSize4;
	ModelHeaderNew.nFeatKind=ModelHeaderOld.nMfccKind;
	ModelHeaderNew.nTotalParamSize=sizeof(float)*m_nMixNum*(m_nVecSize4*2+2);

	WriteOpen(fpModel,argv[2]);
	fwrite(&ModelHeaderNew,sizeof(GMMFileHeaderNew),1,fpModel);
	fseek(fpModel,FILE_HEADER_SIZE,SEEK_SET);
	int nNumRead;
	nNumRead=fwrite(m_pGmmModel[0].pfWeight,sizeof(float),m_nMixNum,fpModel);
	if (nNumRead!=m_nMixNum)
	{
		fclose(fpModel);	
		exit(1);
	}
	nNumRead=fwrite(m_pGmmModel[0].pfMeanBuf,sizeof(float),m_nVecSize4*m_nMixNum,fpModel);
	if (nNumRead!=m_nVecSize4*m_nMixNum)
	{
		fclose(fpModel);	
		exit(1);
	}
	nNumRead=fwrite(m_pGmmModel[0].pfDiagCovBuf,sizeof(float),m_nVecSize4*m_nMixNum,fpModel);
	if (nNumRead!=m_nVecSize4*m_nMixNum)
	{
		fclose(fpModel);	
		exit(1);
	}
	float temp;
	for (int i=0;i<m_nMixNum;i++){
		temp=float(m_pGmmModel[0].pGauss[i].dMat);
		fwrite(&temp,sizeof(float),1,fpModel);
	}
	fclose(fpModel);
	FreeGaussMixModel(m_pGmmModel,1);
	
}


void AllocGaussMixModel(GaussMixModel *p_pModel,int p_nMixNum,int p_nVecSize)
{
	ASSERT2(p_nMixNum,"Error call AllocGaussMixModel() : p_nMixNum<0!");
	ASSERT2(p_nVecSize,"Error call AllocGaussMixModel() : p_nVecSize<0!");
	
	// 分配
	p_pModel->pfWeight = (float *)Malloc(p_nMixNum*sizeof(float));
	p_pModel->pGauss = (GaussPDF*)Malloc (p_nMixNum*sizeof(GaussPDF));

	p_pModel->pfMeanBuf = (float *)Malloc(p_nMixNum*sizeof(float)*p_nVecSize,true);
	p_pModel->pfDiagCovBuf = (float *)Malloc(p_nMixNum*sizeof(float)*p_nVecSize,true);
	for(int m=0;m<p_nMixNum;m++)
	{
		p_pModel->pGauss[m].pfMean = p_pModel->pfMeanBuf + m*p_nVecSize;
		p_pModel->pGauss[m].pfDiagCov = p_pModel->pfDiagCovBuf + m*p_nVecSize;
	}

	// 初始化
	p_pModel->pfWeight[0] = 1.f;
}
void FreeGaussMixModel(GaussMixModel *p_gmmParam,int p_nModelNum)
{
	if (p_gmmParam)
	{
		ASSERT2(p_nModelNum,"Error call FreeGaussMixModel() : p_nModelNum<0!");
		for (int i=0;i<p_nModelNum;i++)
		{
			if (p_gmmParam[i].pGauss)		
			{
				Free(p_gmmParam[i].pfMeanBuf);
				Free(p_gmmParam[i].pfDiagCovBuf);

				Free(p_gmmParam[i].pGauss);

				p_gmmParam[i].pfMeanBuf = NULL;
				p_gmmParam[i].pfDiagCovBuf = NULL;
				p_gmmParam[i].pGauss = NULL;
			}
			Free(p_gmmParam[i].pfWeight);
			p_gmmParam[i].pfWeight = NULL;
		}
		Free(p_gmmParam);
		p_gmmParam = NULL;
	}
}
