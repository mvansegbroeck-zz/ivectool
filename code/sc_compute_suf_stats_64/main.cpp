#include "sc_compute_suf_stats.h"

int main(int argc, char *argv[])
{
	if (argc==1) {
		printf("Usage:\n");
		printf("%s ubmfile flist.list mixnum vecsize outdir\n");
	}

	char * ubmfile = argv[1];
	char * filelist = argv[2];
	char * outdir = argv[5];
	int  p_nMixNum = atoi(argv[3]);
	int  p_nVecSize = atoi(argv[4]);

	calc_suf_stats *fa_stats = new calc_suf_stats();

	char buf[1024] = {0};
	int  nuttranceNum=0;
	FILE *fp = fa_fopen(filelist, "rb");
	while (fgets(buf, 1024, fp)){
		int referid = 0;
		char file[256] = {0};;
		char refername[256] = {0};
		sscanf(buf, "%s", file);
		nuttranceNum++;
	}
	fclose(fp);

	if (!fa_stats->Initalize(p_nMixNum, p_nVecSize, nuttranceNum)) {
		printf("Error in fa_train->Initalize()!\n");
		return -1;
	}

	if (!fa_stats->LoadGMMModel(ubmfile)) {
		printf("Error in fa_train->LoadGMMModel()!\n");
		return -1;
	}

	char *pPos;
	char256 cFileName;
	pPos = strrchr(filelist,'/');
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
		strcpy(cFileName,filelist);
		pPos=strchr(cFileName,'.');
		if (pPos!=NULL)
			*pPos='\0';
	}
	char outname[256];
	sprintf(outname, "%s/%s.stats", outdir, cFileName);

	if (!fa_stats->compute_stats(filelist, outname)) {
		printf("Error in fa_stats->compute_stats()!\n");
		return -1;
	}

	

//	if (!fa_stats->save_stats(outname)) {
//		printf("Error in save stats!\n");
//		return false;
//	}

//     Gave Segmentation Fault      
//	delete[] fa_stats;

	return 0;
}

