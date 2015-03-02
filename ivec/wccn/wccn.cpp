#include "comm.h"
#include "stdio.h"
struct SuperHead 
{
	int nSmpNum;
	int nSmpDim;

	short hh;
	short tt;
};

int main(int argc, char* argv[])
{
	if(argc!=4)
	{
		printf("Usage:\n");
		printf("%s normsdc.list wccn_matrix.txt superdim\n", argv[0]);
		return -1;
	}
	char *superflist=argv[1];
	char *outmatrix=argv[2];
	int  nDim = atoi(argv[3]); 

	FILE *fp=fopen(superflist, "rb");
	if (fp==NULL)
	{
		printf("% can not be opened!\n", superflist);
		return -1;
	}
	char buf[256] = {0};
	SuperHead sHead;
	
	int spknum=0;
	while (fgets(buf, 4096, fp)){
		char supername[256] = {0};
		sscanf(buf, "%s", supername);
		
		FILE *fp_super = fopen(supername, "rb");
		if(fp_super==NULL)
		{
			printf("%s is not exist!\n", supername);
			return -1;
		}
		fclose(fp_super);
		spknum++;
	}

	printf("Total %d speakers!\n", spknum);

	FMatrix *p_wi=(FMatrix *)malloc(sizeof(FMatrix)*spknum);
	if(p_wi==NULL)
	{
		return -1;
	}

	rewind(fp);
	int index=0;
	while (fgets(buf, 4096, fp)){
		char supername[256] = {0};
		sscanf(buf, "%s", supername);
		
		FILE *fp_super = fopen(supername, "rb");
		if(fp_super==NULL)
		{
			printf("%s is not exist!\n", supername);
			return -1;
		}

		fread(&sHead, sizeof(SuperHead), 1, fp_super);
		int nsuperNum=sHead.nSmpNum;
		if(sHead.nSmpDim!=nDim+1)
		{
			printf("error super dim!\n");
			return -1;
		}

		if(!AllocMatrix(p_wi[index], nsuperNum, nDim))
		{
			return -1;
		}
		for(int i=0; i<nsuperNum; i++)
		{
			float first_val=-1;
			fread(&first_val, sizeof(float), 1, fp_super);
			fread(&p_wi[index].pBuf[i*nDim], sizeof(float), nDim, fp_super);
		}
		fclose(fp_super);
		index++;
	}
	fclose(fp);

	FMatrix pp_mean_ws; // mean of each speaker
	if(!AllocMatrix(pp_mean_ws, spknum, nDim))
	{
		printf("Error malloc pp_mean_ws!\n");
		return -1;
	}
	ZeroMatrix(pp_mean_ws);

	int i, j,m,k;
	for(i=0; i<spknum; i++)
	{
		for(j=0; j<p_wi[i].nRowDim; j++)
		{
			for(k=0; k<nDim; k++)
			{
				pp_mean_ws.pBuf[i*nDim+k] += p_wi[i].pBuf[j*nDim+k];
			}
		}
		
		for(k=0; k<nDim; k++)
		{
			pp_mean_ws.pBuf[i*nDim+k] /= p_wi[i].nRowDim;
		}
	}


	FMatrix p_Sw;
	if(!AllocMatrix(p_Sw, nDim, nDim))
	{
		return -1;
	}
	ZeroMatrix(p_Sw);

	FMatrix p_Sw_sub;
	if(!AllocMatrix(p_Sw_sub, nDim, nDim))
		return -1;
	
	float *p_wi_temp_2=(float *)malloc(sizeof(float)*nDim);

	for(i=0; i<spknum; i++)
	{
		ZeroMatrix(p_Sw_sub);
		for(j=0; j<p_wi[i].nRowDim; j++)
		{
			for(k=0; k<nDim; k++)
			{
				p_wi_temp_2[k]=p_wi[i].pBuf[j*nDim+k]-pp_mean_ws.pBuf[i*nDim+k];
			}
			for(m=0; m<nDim; m++)
			{
				for(k=0; k<nDim; k++)
				{
					p_Sw_sub.pBuf[m*nDim+k] += p_wi_temp_2[m]*p_wi_temp_2[k];
				}
			}
		}
		for(j=0; j<nDim*nDim; j++)
			p_Sw.pBuf[j] += p_Sw_sub.pBuf[j]/p_wi[i].nRowDim;
	}
	for(j=0; j<nDim*nDim; j++)
		p_Sw.pBuf[j] /= spknum;

	FreeMatrix(p_Sw_sub);
	free(p_wi_temp_2);
	FreeMatrix(pp_mean_ws);
	for(i=0; i<spknum; i++)
	{
		FreeMatrix(p_wi[i]);
	}
	free(p_wi);

	if(!mathRealMatrixSymmetryInversion(p_Sw.pBuf, nDim))
			return false;

	printf("start Cholesky Decomp...\n");
	FMatrix f_A;
	if(!AllocMatrix(f_A, nDim, nDim))
		return false;
	ZeroMatrix(f_A);
	IppStatus status = ippStsNullPtrErr;
	status = ippmCholeskyDecomp_m_32f(p_Sw.pBuf, sizeof(float)*nDim, sizeof(float), \
		f_A.pBuf, sizeof(float)*nDim, sizeof(float), nDim);

	printf("end Cholesky Decomp.\n");

	FreeMatrix(p_Sw);
	
	FMatrix f_A_trans;
	if(!AllocMatrix(f_A_trans, nDim, nDim))
	{
		return -1;
	}
	if(!matrixTranspos(f_A, f_A_trans))
	{
		return -1;
	}

	for(i=0; i<nDim; i++)
	{
		f_A_trans.pBuf[i*nDim+i]=1/f_A_trans.pBuf[i*nDim+i];
	}
	
	
	FILE *fp_temp=fopen(outmatrix, "wt");
	
	
	
	for(int ii=0; ii<nDim; ii++)
	{
		for(int jj=0; jj<nDim; jj++)
		{
			fprintf(fp_temp, "%f ", f_A_trans.pBuf[ii*nDim+jj]);
		}
		fprintf(fp_temp, "\n");
	}
	fclose(fp_temp);
	
	
	FreeMatrix(f_A);
	FreeMatrix(f_A_trans);
	return 0;
}
