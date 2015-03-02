#include "sc_compute_suf_stats.h"

calc_suf_stats::calc_suf_stats()
{
	m_pfIPPBuf_VecSizeStore=NULL;
	m_cUttNames=NULL;
}
calc_suf_stats::~calc_suf_stats()
{
	if (m_pfIPPBuf_VecSizeStore!=NULL) {
		Free(m_pfIPPBuf_VecSizeStore);
		m_pfIPPBuf_VecSizeStore=NULL;
	}
	
	FreeGaussModel(m_UBMGmmModel);
	FreeBWStatisBuf(m_Statistic);

	FreeVector(N_utt);
	FreeVector(F_utt);
	free(m_cUttNames);
	free(m_Utt_nFrmNums);
	m_Utt_nFrmNums=NULL;
}
bool calc_suf_stats::AllocGaussModel(GaussMixModel &p_GmmModel)
{
	p_GmmModel.pfMatBuf=Malloc_32f(m_nMixNum);
	if(!p_GmmModel.pfMatBuf)
		return false;

	p_GmmModel.pfWeightBuf=Malloc_32f(m_nMixNum);
	if(!p_GmmModel.pfWeightBuf)
		return false;

	p_GmmModel.pfMeanBuf=Malloc_32f(m_nMixNum*m_nVecSizeStore);
	if(!p_GmmModel.pfMeanBuf)
		return false;

	p_GmmModel.pfDiagCovBuf=Malloc_32f(m_nMixNum*m_nVecSizeStore);
	if(!p_GmmModel.pfDiagCovBuf)
		return false;

	return true;
}

void calc_suf_stats::FreeGaussModel(GaussMixModel &p_GmmModel)
{
	if(p_GmmModel.pfMeanBuf!=NULL) {
		Free(p_GmmModel.pfMeanBuf);
		p_GmmModel.pfMeanBuf=NULL;
	}

	if(p_GmmModel.pfDiagCovBuf!=NULL) {
		Free(p_GmmModel.pfDiagCovBuf);
		p_GmmModel.pfDiagCovBuf=NULL;
	}

	if(p_GmmModel.pfWeightBuf!=NULL) {
		Free(p_GmmModel.pfWeightBuf);
		p_GmmModel.pfWeightBuf=NULL;
	}
	
	if(p_GmmModel.pfMatBuf!=NULL) {
		Free(p_GmmModel.pfMatBuf);
		p_GmmModel.pfMatBuf=NULL;
	}
}




bool calc_suf_stats::Initalize(int p_nMixNum, int p_nVecSize, int p_nUttNum)
{	
	m_nMixNum     = p_nMixNum;
	m_nVecSize    = p_nVecSize;
	m_nVecSizeStore	 = ALIGN_4F(m_nVecSize);

	m_nUttNum=p_nUttNum;
	
	if(!AllocGaussModel(m_UBMGmmModel))
		return false;

//	if (!AllocMatrix(N_utt, m_nUttNum, m_nMixNum)) {
//		return false;
//	}

//	ZeroMatrix(N_utt);

	if(!AllocVector(N_utt, m_nMixNum)) {
		return false;
	}
	ZeroVector(N_utt);

	if (!AllocVector(F_utt, m_nMixNum*m_nVecSize)) {
		return false;
	}
	ZeroVector(F_utt);

	m_cUttNames = (char32 *)malloc(sizeof(char32)*m_nUttNum);
	if (m_cUttNames==NULL) {
		return false;
	}

	m_Utt_nFrmNums = (int *)malloc(sizeof(int)*m_nUttNum);
	if (m_Utt_nFrmNums==NULL) {
		return false;
	}

	if (!AllocBWStatisBuf(m_Statistic)) 
		return false;
	
	printf("Initializing for computing sufficient statistics successfully!\n");
	return true;
}



bool calc_suf_stats::LoadGMMModel(char *p_pcModelFile)
{
	if (NULL==p_pcModelFile)	return false;

	FILE *fpModel;
	fpModel=fopen(p_pcModelFile,"rb");
	if(fpModel==NULL) {
		printf("Error open %s for read.\n",p_pcModelFile);
		return false;
	}

	GMMFileHeader  ModelHeader;	
	fread(&ModelHeader,sizeof(GMMFileHeader),1,fpModel);
	if(m_nVecSize != ModelHeader.nVecSize) {
		fclose(fpModel);  
		return false;
	}
	if(m_nMixNum  != ModelHeader.nMixNum) { 
		fclose(fpModel);  
		return false;
	}
	if (m_nVecSizeStore != ModelHeader.nVecSizeStore) {
		fclose(fpModel);  
		return false;
	}

	if (m_nVecSizeStore%4 != 0) {
		fclose(fpModel);  
		return false;
	}

	m_nFeatKind = ModelHeader.nFeatKind;

	fseek(fpModel, FILE_HEADER_SIZE, SEEK_SET);	
	int nNumRead;
	nNumRead = fread(m_UBMGmmModel.pfWeightBuf,sizeof(float),m_nMixNum,fpModel);
	if (nNumRead!=m_nMixNum) {
		fclose(fpModel);	
		return false;
	}
	nNumRead = fread(m_UBMGmmModel.pfMeanBuf,sizeof(float),ModelHeader.nVecSizeStore*m_nMixNum,fpModel);
	if (nNumRead!=ModelHeader.nVecSizeStore*m_nMixNum) {
		fclose(fpModel);	
		return false;
	}
	nNumRead=fread(m_UBMGmmModel.pfDiagCovBuf,sizeof(float),ModelHeader.nVecSizeStore*m_nMixNum,fpModel);
	if (nNumRead!=ModelHeader.nVecSizeStore*m_nMixNum) {
		fclose(fpModel);	
		return false;
	}
	nNumRead=fread(m_UBMGmmModel.pfMatBuf,sizeof(float),m_nMixNum,fpModel);
	if (nNumRead!=m_nMixNum) {
		fclose(fpModel);	
		return false;
	}	
	fclose(fpModel);
	printf("Loading UBM parameters successfully!\n");	
	return true;
}

bool calc_suf_stats::AllocBWStatisBuf(GMMStatis &p_Statistic)
{
	FreeBWStatisBuf(p_Statistic);

	p_Statistic.pGauss = (Gaussstatis *)malloc(m_nMixNum*sizeof(Gaussstatis));
	if (!p_Statistic.pGauss) {
		return false;
	}

	p_Statistic.Facc = Malloc_32f(m_nMixNum*m_nVecSizeStore);
	if (!p_Statistic.Facc) {
		return false;
	}

	p_Statistic.Sacc = Malloc_32f(m_nMixNum*m_nVecSizeStore);
	if (!p_Statistic.Sacc) {
		return false;
	}

	float *pPos_F, *pPos_S;
	pPos_F = p_Statistic.Facc;
	pPos_S = p_Statistic.Sacc;
	for(int m=0;m<m_nMixNum;m++)
	{
		p_Statistic.pGauss[m].N_occ = 0.0;
		p_Statistic.pGauss[m].p_F = pPos_F;
		p_Statistic.pGauss[m].p_S = pPos_S;

		pPos_F += m_nVecSizeStore;
		pPos_S += m_nVecSizeStore;
	}
	
	m_pfIPPBuf_VecSizeStore = Malloc_32f(m_nVecSizeStore);

	return true;
}

void calc_suf_stats::FreeBWStatisBuf(GMMStatis &p_Statistic)
{
	if (p_Statistic.Facc!=NULL)
		Free(p_Statistic.Facc);
	if (p_Statistic.Sacc!=NULL)
		Free(p_Statistic.Sacc);
	if (p_Statistic.pGauss!=NULL)
		free(p_Statistic.pGauss);

	if (m_pfIPPBuf_VecSizeStore!=NULL) {
		Free(m_pfIPPBuf_VecSizeStore);
	}

	p_Statistic.Sacc = NULL;
	p_Statistic.Facc = NULL;
	p_Statistic.pGauss=NULL;
	m_pfIPPBuf_VecSizeStore=NULL;
}

void calc_suf_stats::ResetBWStatisBuf(GMMStatis &p_Statistic)
{
	if (p_Statistic.Facc!=NULL)
		Zero_32f(p_Statistic.Facc, m_nMixNum*m_nVecSizeStore);

	if (p_Statistic.Sacc!=NULL)
		Zero_32f(p_Statistic.Sacc, m_nMixNum*m_nVecSizeStore);

	if (p_Statistic.pGauss!=NULL) {
		for(int m=0;m<m_nMixNum;m++)
			p_Statistic.pGauss[m].N_occ = 0.0;	
	}
}

bool calc_suf_stats::BWStatisticEst(float *p_fFeatBuf, int p_nFrmNum)
{
	if (!p_fFeatBuf) {
		return false;
	}
	if (p_nFrmNum<=0) {
		return false;
	}

	ResetBWStatisBuf(m_Statistic);

	float*  pfProb_Buf;
	float** ppfProb_BufBuf;

	pfProb_Buf = (float *)malloc(sizeof(float)*p_nFrmNum);
	if (pfProb_Buf==NULL)
	{
		printf("Error malloc m_pfProb_Buf[]!\n");
		exit(-1);
	}
		
	ppfProb_BufBuf = (float **)malloc(m_nMixNum*sizeof(float *));
	if (ppfProb_BufBuf==NULL) 
	{
		printf("Error malloc m_ppfProb_BufBuf!\n");
		return false;
	}

	int i;
	for (i=0; i<m_nMixNum; i++)
	{
		ppfProb_BufBuf[i] = (float *)malloc(p_nFrmNum*sizeof(float));
		if (ppfProb_BufBuf[i]==NULL) 
		{
			printf("Error malloc m_ppfProb_BufBuf[%d]!\n",i);
			return false;
		}
	}

	// 计算第0个高斯分量与所有帧间的似然值

	int m;
	ippsLogGauss_32f_D2(p_fFeatBuf,
				m_nVecSizeStore,
				m_UBMGmmModel.pfMeanBuf,
				m_UBMGmmModel.pfDiagCovBuf,
				m_nVecSize,
				pfProb_Buf,
				p_nFrmNum,
				m_UBMGmmModel.pfMatBuf[0]);
        

	// 计算第m个高斯分量与所有帧间的似然值 --> m_ppfProb_BufBuf[m][t]
	// 计算当前模型与特征帧间的似然值 --> m_pfProb_Buf[t]

	ippsCopy_32f(pfProb_Buf,ppfProb_BufBuf[0],p_nFrmNum);
	
	for(m=1;m<m_nMixNum;m++)
	{
		ippsLogGauss_32f_D2(p_fFeatBuf,
				m_nVecSizeStore,
				&m_UBMGmmModel.pfMeanBuf[m*m_nVecSizeStore],
				&m_UBMGmmModel.pfDiagCovBuf[m*m_nVecSizeStore],
				m_nVecSize,				
				(float *)ppfProb_BufBuf[m],
				p_nFrmNum,
				m_UBMGmmModel.pfMatBuf[m]);

// /* DEBUGGING START */
//                 if (m==m_nMixNum-1)
//                 {
//                     printf("data=[\n");
//                     for(int k=0; k<36;k++) 
//                     {
//                       printf("%.5f \t %.5f \t %.5f\n",p_fFeatBuf[k],m_UBMGmmModel.pfMeanBuf[m*m_nVecSizeStore+k],m_UBMGmmModel.pfDiagCovBuf[m*m_nVecSizeStore+k]);
// 
//                     }
//                     printf("];\n");
//                     printf("constant=%.5f;\n",m_UBMGmmModel.pfMatBuf[m]);
//                     printf("prob=%.5f;\n",*ppfProb_BufBuf[m]);
//                     // x=data(:,1);mn=data(:,2);cv=data(:,3);prob=constant+log(exp(-0.5*sum(((x-mn).^2).*cv,1)))
//                 }
// /* DEBUGGING END */

		ippsLogAdd_32f(ppfProb_BufBuf[m],
				pfProb_Buf,
				p_nFrmNum,
				ippAlgHintNone);
	}

	float fProbSum = 0.f;
	ippsSum_32f(pfProb_Buf, p_nFrmNum, &fProbSum, ippAlgHintNone);

	printf("This file : Prob=%12.6f, fmNum=%5d, avgProb=%12.6f\n",
		fProbSum, p_nFrmNum, fProbSum/(float)p_nFrmNum);

	float *pCurFrame;
	float *featvector = Malloc_32f(m_nVecSizeStore);
	if (!featvector) {
		return false;
	}

	for(m=0;m<m_nMixNum;m++) {
		// ppfProb_BufBuf[m][i]-pfProb_Buf[i]  //m高斯,i代表帧的ID
		ippsSub_32f_I(pfProb_Buf, ppfProb_BufBuf[m], p_nFrmNum);
		// exp(ppfProb_BufBuf[m][i]-pfProb_Buf[i])
		ippsExp_32f_I(ppfProb_BufBuf[m],p_nFrmNum);  //第i帧数据在第m个高斯分量上的后验

		fProbSum = 0.0f;
		ippsSum_32f(ppfProb_BufBuf[m],p_nFrmNum,&fProbSum,ippAlgHintNone);
		m_Statistic.pGauss[m].N_occ = fProbSum;

		pCurFrame = p_fFeatBuf;
		for (int i=0;i<p_nFrmNum;i++) {
			if (ppfProb_BufBuf[m][i]<MIN_PROB) {
				pCurFrame += m_nVecSizeStore;
				continue;
			}
			memcpy(featvector, pCurFrame, sizeof(float)*m_nVecSizeStore);
			// meanAcc (Fc) first-order statistics
//			ippsSub_32f_I(&m_UBMGmmModel.pfMeanBuf[m*m_nVecSizeStore], featvector, m_nVecSizeStore);
			ippsMulC_32f(featvector, ppfProb_BufBuf[m][i], m_pfIPPBuf_VecSizeStore, m_nVecSizeStore);
			ippsAdd_32f_I(m_pfIPPBuf_VecSizeStore, m_Statistic.pGauss[m].p_F, m_nVecSizeStore);
//			ippsAdd_32f_I(featvector, m_Statistic.pGauss[m].p_F, m_nVecSizeStore);
	
			// varAcc (Sc) diagonal matrix, second-order statistics
//			ippsMul_32f(featvector, featvector, m_pfIPPBuf_VecSizeStore, m_nVecSizeStore);
//			ippsMulC_32f_I(ppfProb_BufBuf[m][i], m_pfIPPBuf_VecSizeStore, m_nVecSizeStore);
//			ippsAdd_32f_I(m_pfIPPBuf_VecSizeStore, m_Statistic.pGauss[m].p_S, m_nVecSizeStore);

			pCurFrame += m_nVecSizeStore;
		}
	}

	Free(featvector);
	//释放中间统计量
	free(pfProb_Buf);

	for(i=0; i<m_nMixNum; i++)
			free(ppfProb_BufBuf[i]);
	free(ppfProb_BufBuf);

	return true;	
}

bool calc_suf_stats::compute_stats(char *UttListfile, char *outstats)
{
//	FILE *fp = fa_fopen(UttListfile, "rb");

//	char buf[1024] = {0};

	printf("Total %d utterances!\n", m_nUttNum);
	
	FILE *fp_out = fopen(outstats, "wb");
	if (fp_out==NULL) {
		return false;
	}

	int flag=1;
	fwrite(&m_nMixNum, sizeof(int), 1, fp_out);
	fwrite(&m_nVecSize, sizeof(int), 1, fp_out);
	fwrite(&m_nUttNum, sizeof(int), 1, fp_out);
	fwrite(&flag, sizeof(int), 1, fp_out);
//	fwrite(N_utt.pBuf, sizeof(float), m_nUttNum*m_nMixNum, fp);
	fseek(fp_out, 4*sizeof(int)+m_nUttNum*m_nMixNum*sizeof(float), SEEK_SET);
	fwrite(&flag, sizeof(int), 1, fp_out);
//	fwrite(F_utt.pBuf, sizeof(float), m_nUttNum*m_nMixNum*m_nVecSize, fp);
//	fwrite(&flag, sizeof(int),1, fp);
//	fwrite(m_Utt_nFrmNums, sizeof(int), m_nUttNum, fp);
//	fwrite(&flag, sizeof(int),1, fp);

//	for (int i=0; i<m_nUttNum; ++i) {
//		fwrite(&m_cUttNames[i], sizeof(char32), 1, fp);
//	}


	FILE *fp = fa_fopen(UttListfile, "rb");
	char buf[1024] = {0};

	int index = 0;
	while (fgets(buf, 1024, fp)){
		char file[256] = {0};
		sscanf(buf, "%s", file);

		printf("Reading the %dth feature file %s\n", index+1, file);
		
		float *pFeatData=NULL;
		Feature_BaseInfo sHead;
		if(!ReadFeatFile(file, pFeatData, sHead))
			return false;
		if (sHead.nFeatKind!=m_nFeatKind) {
			printf("Error feat kind!\n");
			return false;
		}
		printf("Processing...\n");
		if(!BWStatisticEst(pFeatData, sHead.nFrameNum))
			return false;
		free(pFeatData);
		pFeatData=NULL;

		int i;
		for (i=0; i<m_nMixNum; i++) {
			N_utt.pBuf[i] = m_Statistic.pGauss[i].N_occ;
			memcpy(&F_utt.pBuf[i*m_nVecSize], m_Statistic.pGauss[i].p_F, \
				sizeof(float)*m_nVecSize);
		}

		fseek(fp_out, 4*sizeof(int)+index*m_nMixNum*sizeof(float), SEEK_SET);
		fwrite(N_utt.pBuf, sizeof(float), m_nMixNum, fp_out);
		fseek(fp_out, 4*sizeof(int)+m_nUttNum*m_nMixNum*sizeof(float)+sizeof(int)+index*m_nMixNum*m_nVecSize*sizeof(float), SEEK_SET);
		fwrite(F_utt.pBuf, sizeof(float), m_nMixNum*m_nVecSize, fp_out);

		char *pPos;
		char256 cFileName;
		pPos = strrchr(file,'\\');
		if (pPos!=NULL)
		{
			pPos++;
			strcpy(cFileName,pPos);
			pPos=strchr(cFileName,'.');
			if (pPos!=NULL)
				*pPos='\0';
		}
		else
		{
			strcpy(cFileName,file);
			pPos=strchr(cFileName,'.');
			if (pPos!=NULL)
				*pPos='\0';
		}

		m_Utt_nFrmNums[index]=sHead.nFrameNum;
		strcpy(m_cUttNames[index], cFileName);
		index++;
	}
	
	fwrite(&flag, sizeof(int), 1, fp_out);
	fwrite(m_Utt_nFrmNums, sizeof(int), m_nUttNum, fp_out);
	fwrite(&flag, sizeof(int),1, fp_out);

	for (int i=0; i<m_nUttNum; ++i) {
		fwrite(&m_cUttNames[i], sizeof(char32), 1, fp_out);
	}

	fclose(fp_out);
	fclose(fp);

	return true;
}

bool calc_suf_stats::save_stats(char *outfilename)
{
	FILE *fp = fopen(outfilename, "wb");
	if (fp==NULL) {
		return false;
	}

	int flag=1;
	fwrite(&m_nMixNum, sizeof(int), 1, fp);
	fwrite(&m_nVecSize, sizeof(int), 1, fp);
	fwrite(&m_nUttNum, sizeof(int), 1, fp);
	fwrite(&flag, sizeof(int), 1, fp);
	fwrite(N_utt.pBuf, sizeof(float), m_nUttNum*m_nMixNum, fp);
	fwrite(&flag, sizeof(int), 1, fp);
	fwrite(F_utt.pBuf, sizeof(float), m_nUttNum*m_nMixNum*m_nVecSize, fp);
	fwrite(&flag, sizeof(int), 1, fp);
	fwrite(m_Utt_nFrmNums, sizeof(int), m_nUttNum, fp);
	fwrite(&flag, sizeof(int),1, fp);

	for (int i=0; i<m_nUttNum; ++i) {
		fwrite(&m_cUttNames[i], sizeof(char32), 1, fp);
	}

	fclose(fp);

	return true;
}
