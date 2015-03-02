/******************************DOCUMENT*COMMENT***********************************
*D
*D �ļ�����            : memory_srlr.cpp
*D
*D ��Ŀ����            : 
*D
*D �汾��              : 1.1.0001
*D
*D �ļ�����            :Memory Management Modules
*D
*D
*D �ļ��޸ļ�¼
*D ------------------------------------------------------------------------------ 
*D �汾��       �޸�����       �޸���     �Ķ�����
*D ------------------------------------------------------------------------------ 
*D 1.1.0001     2007.01.01     plu        �޸��ļ�
*D*******************************************************************************/
#include <memory.h>
#include "comm_srlr.h"

// ����sz���ֽ�
void *Malloc(int sz,bool clear=false)
{
	ASSERT3(sz>0,"Error call Malloc(int sz,const bool clear) : sz=%d!",sz);	// 2007.02.12 plu : add

	void *pt=malloc(sz); 
	if (pt==NULL)
	{ 
	   printf ("Memory can't allocate %d bytes!\a\n", sz);					// 2007.02.12 plu : add
	   char256 str;
	   sprintf (str, "Memory allocate %d bytes", sz);
	   exit(-1);
	}

	if (clear)
	{
		memset(pt,0,sz);
	}
	
	return pt;
}

// ����p_nBlockNum*p_nSize���ֽ�
void *Malloc(int p_nBlockNum,int p_nBlockSize,bool clear=false)
{
	ASSERT3(p_nBlockNum*p_nBlockSize>0,
		"Error call Malloc(int p_nBlockNum,int p_nBlockSize, const bool clear) : p_nBlockNum*p_nBlockSize=%d!",
		p_nBlockNum*p_nBlockSize);	// 2007.02.12 plu : add
	
	return Malloc(p_nBlockNum*p_nBlockSize,clear);
}

// ����[p_nRow][p_nCol*p_nBlockSize]���ֽ�
void **Malloc(int p_nRow,int p_nCol,int p_nBlockSize,bool clear=false)
{
   void **buf=(void **)Malloc(sizeof(void *)*p_nRow,clear);

   for (int i=0;i<p_nRow;i++) 
	   buf[i]=Malloc(p_nCol,p_nBlockSize,clear);
   
   return buf;
}

// ����[p_nRow1][p_nRow2][p_nCol*p_nBlockSize]���ֽ�
void ***Malloc(int p_nRow1,int p_nRow2,int p_nCol,int p_nBlockSize,bool clear=false)
{
	void ***buf=(void ***)Malloc(sizeof(void **)*p_nRow1,clear);
	for (int i=0;i<p_nRow1;i++)
	{
		buf[i]=(void **)Malloc(sizeof(void *)*p_nRow2,clear);

		for (int j=0;j<p_nRow2;j++)
			buf[i][j]=Malloc(p_nCol*p_nBlockSize,clear);
	}
	return buf;
}


void *Malloc32(int sz,bool clear=false)
{
   char *ptr = (char*)Malloc(sz+32,clear);
   void *pt = (void*)(ptr + 32 - ((ptr - (char*)NULL) & 31));
   ((void**)pt)[-1]=ptr;

   return pt;
}

void *Malloc32(int n1,int itsz,bool clear=false)
{
   return Malloc32(n1*itsz,clear);
}

void Free(void *ptr)
{ 
	if (ptr!=NULL)					// 2007.02.12 plu : ����if(ptr) �ж�
	{
		free(ptr);
		ptr = NULL;
	}
}

void Free32(void *ptr) 
{ 
	free(((void**)ptr)[-1]); 
}

void Free(void **ptr,int n1)
{
	if (ptr)					// 2007.02.12 plu : ����if(ptr) �ж�
	{
		if (n1<=0)
			printf("Warning:  Free(void **ptr,int n1) n1=%d!\a\n",n1);		// 2007.02.12 plu : add

		for (int i=0;i<n1;i++) 
			Free(ptr[i]);

		Free(ptr);
	}
}

 void Free(void ***ptr,int n1,int n2) throw() 
{
	if (ptr)					// 2007.02.12 plu : ����if(ptr) �ж�
	{
		if (n1<=0)
			WARNING2("Error call Free(void **ptr,int n1,int n2) : n1=%d!",n1);		// 2007.02.12 plu : add
		if (n2<=0)
			WARNING2("Error call Free(void **ptr,int n1,int n2) : n2=%d!",n2);		// 2007.02.12 plu : add
		for (int i=0;i<n1;i++)
		{
			for (int j=0;j<n2;j++) Free(ptr[i][j]);
			Free(ptr[i]);
		}
		Free(ptr);
   }
}



