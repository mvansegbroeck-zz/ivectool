/******************************DOCUMENT*COMMENT***********************************
*D
*D 文件名称            : main.cpp
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
*D 1.1.0001     2007.12.25     plu        创建文件
*D*******************************************************************************/

#include "comm_srlr.h"
#include "GMM_comm.h"

int main(int argc,char *argv[])
{
	if (argc==1)
	{
		printf("Usage: %s -i GMMModelFile -o TxtFile\n",argv[0]);
		return 0;
	}

	char *pcModelFile,*pcTxtFile;
	pcTxtFile = pcModelFile = NULL;

	for (int i=1;i<argc;i++) 
	{
		if (!strcmp(argv[i],"-i"))
			pcModelFile = argv[++i];
		else if (!strcmp(argv[i],"-o"))
			pcTxtFile = argv[++i];
	}
	
	if (pcTxtFile==NULL)
	{
		printf("pls input param : -o TxtFile\n");	return 0;
	}
	if (pcModelFile==NULL)
	{
		printf("pls input param : -i CohortModelFile\n");	return 0;
	}
	//////////////////////////////////////////////////////
	FILE *fpMdl,*fpTxt;
	ReadOpen(fpMdl,pcModelFile);
	TextWriteOpen(fpTxt,pcTxtFile);

	GMMFileHeader sHeader;
	fread(&sHeader,sizeof(GMMFileHeader),1,fpMdl);
	
	float *pfParamBuf=(float *)malloc(sHeader.nTotalParamSize);
	ASSERT3(pfParamBuf,"Error malloc [%d]!",sHeader.nTotalParamSize);
	fseek(fpMdl,FILE_HEADER_SIZE,SEEK_SET);
	fread(pfParamBuf,1,sHeader.nTotalParamSize,fpMdl);
	fclose(fpMdl);
	
	printf("MixtureNum     = %d\n",sHeader.nMixNum);
	printf("VecSize        = %d\n",sHeader.nVecSize);
	printf("VecSizeStore   = %d\n",sHeader.nVecSizeStore);
	printf("FeatKind       = %d\n",sHeader.nFeatKind);
	printf("nTotalParamSize       = %d\n",sHeader.nTotalParamSize);

	fprintf(fpTxt,"MixtureNum     = %d\n",sHeader.nMixNum);
	fprintf(fpTxt,"VecSize        = %d\n",sHeader.nVecSize);
	fprintf(fpTxt,"VecSizeStore   = %d\n",sHeader.nVecSizeStore);
	fprintf(fpTxt,"nTotalParamSize       = %d\n",sHeader.nTotalParamSize);

	int nStep=sHeader.nVecSizeStore;
	float *pfCur;
	pfCur = pfParamBuf;

	fprintf(fpTxt,"*************** Weight ****************\n");
	for (int m=0;m<sHeader.nMixNum;m++)
	{
		fprintf(fpTxt,"[mix=%4d] fWeight=%f\n",m+1,*pfCur);
		pfCur++;
	}
	
	fprintf(fpTxt,"*************** Mean ****************\n");
	for (int m=0;m<sHeader.nMixNum;m++)
	{
		fprintf(fpTxt,"[mix=%4d] mean=",m+1);
		for (int nDim=0;nDim<sHeader.nVecSize;nDim++)
		{
			fprintf(fpTxt,"%12.6f",pfCur[nDim]);

		}
		fprintf(fpTxt,"\n");
		pfCur += nStep;
	}
	
	fprintf(fpTxt,"*************** DiagCov ****************\n");
	for (int m=0;m<sHeader.nMixNum;m++)
	{
		fprintf(fpTxt,"[mix=%4d] diagcov=",m+1);
		for (int nDim=0;nDim<sHeader.nVecSize;nDim++)
		{
			fprintf(fpTxt,"%12.6f",pfCur[nDim]);

		}
		fprintf(fpTxt,"\n");
		pfCur += nStep;
	}
	fprintf(fpTxt,"*************** Mat ****************\n");
	for (int m=0;m<sHeader.nMixNum;m++)
	{
		fprintf(fpTxt,"[mix=%4d] fMat=%f\n",m+1,*pfCur);
		if (*pfCur>=0.f)
		{
			printf("Warning : [mix=%4d] fMat=%f\a\n",m+1,*pfCur);
		}
		pfCur++;
	}
	
	fclose(fpTxt);

	//////////////////////////////////////////////////////
	free(pfParamBuf);
	printf("ok!\n");
	return 1;
}
