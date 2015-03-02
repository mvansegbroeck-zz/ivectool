/******************************DOCUMENT*COMMENT***********************************
*D
*D �ļ�����            : memory_srlr.h
*D
*D ��Ŀ����            : 
*D
*D �汾��              : 1.1.0001
*D
*D �ļ�����            : Memory Management Head File
*D
*D
*D �ļ��޸ļ�¼
*D ------------------------------------------------------------------------------ 
*D �汾��       �޸�����       �޸���     �Ķ�����
*D ------------------------------------------------------------------------------ 
*D 1.1.0001     2007.01.01                �޸��ļ�
*D*******************************************************************************/
#ifndef _MEMORY_SRLR_H_
#define _MEMORY_SRLR_H_

extern void *Malloc(int sz,bool clear=false);
extern void *Malloc(int p_nNum,int p_nSize,bool clear=false);
extern void **Malloc(int p_nRow,int p_nCol,int p_nBlockSize,bool clear=false);
extern void ***Malloc(int p_nRow1,int p_nRow2,int p_nCol,int p_nBlockSize,bool clear=false);

extern void Free(void *ptr); 
extern void Free(void **ptr,int n1);
extern void Free(void ***ptr,int n1,int n2);

extern void *Malloc32(int sz,bool clear=false);
extern void *Malloc32(int n1,int itsz,bool clear=false);
extern void Free32(void *ptr);

#endif
