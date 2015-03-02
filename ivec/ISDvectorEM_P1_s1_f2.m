% P1. E STEP OF EM TRAINING

function ISDvectorEM_P1_s1_f2(TrainMatName,TrainLabelName,StasticsMatName,N_Rank,Idx_iteration,nbClassesTrain,n_MeansuperDim,T_init,InvSigma1,matrix1,nMixNum,nVecSize,ProjectMatrixPool,log_N_s_table,log_N_s_table_length,log_N_s_table_interval,W_init,InvSigma2,matrix2)

%load(LastModelMatName);
load(TrainMatName);
label=load(TrainLabelName);

n_TrainSample=size(F,2);
assert(n_MeansuperDim==size(F,1));

N_s=sum(N);

L_train=zeros(nbClassesTrain,n_TrainSample);

for i=1:nbClassesTrain
   idx=find(label==i);
   L_train(i,idx)=1; 
end

X=zeros(N_Rank,n_TrainSample);
t1=zeros(sum(nVecSize*nMixNum),n_TrainSample);
t2=zeros(nbClassesTrain,n_TrainSample);
Matrix_N_E_xx=single(zeros(N_Rank,N_Rank));
Matrix_N_F_E_x=single(zeros(sum(nVecSize*nMixNum),N_Rank));
Matrix_N_L_E_x=single(zeros(nbClassesTrain,N_Rank));

InvSigmaVector1=reshape(InvSigma1,sum(nVecSize*nMixNum),1);

min_log_N_s_table=min(log_N_s_table);
n_trainingSampleIdx=[1:n_TrainSample];    
n_ProjMatrixPoolIdx=[1:n_TrainSample];    
for i=1:n_TrainSample
   ProjMatrixPoolIdx=round((log(N_s(i))-min_log_N_s_table)/log_N_s_table_interval);
   ProjMatrixPoolIdx=min(ProjMatrixPoolIdx,log_N_s_table_length);
   ProjMatrixPoolIdx=max(ProjMatrixPoolIdx,1);

   n_ProjMatrixPoolIdx(i)=ProjMatrixPoolIdx;

   N_s(i)=exp(min(log_N_s_table)+log_N_s_table_interval*(ProjMatrixPoolIdx-1));
   
   t1(:,i)=InvSigmaVector1.*(F(:,i).*N_s(i));
   t2(:,i)=InvSigma2.*(L_train(:,i).*N_s(i));
   
end
t1=T_init'*t1;
t2=W_init'*t2;
x=t1+t2;
% unsupervised
%x=t1;

if length(unique(N_s))==1,
  X=ProjectMatrixPool{n_ProjMatrixPoolIdx(1)}*x;
  for i=1:n_TrainSample
   if((sum(isnan(X(:,i)))~=0)&&(sum(isnan(F(:,i)))~=0))
     now=sprintf('NaN Iter %d SplitName %s Sample %d\n',Idx_iteration,TrainMatName,i)
   end
  end
  Matrix_N_E_xx=ProjectMatrixPool{n_ProjMatrixPoolIdx(1)}.*N_s(1).*n_TrainSample;
else   
  for i=1:n_TrainSample
     X(:,i)=ProjectMatrixPool{n_ProjMatrixPoolIdx(i)}*x(:,i);
     if((sum(isnan(X(:,i)))==0)&&(sum(isnan(F(:,i)))==0))
         Matrix_N_E_xx=Matrix_N_E_xx+ProjectMatrixPool{n_ProjMatrixPoolIdx(i)}.*N_s(i);
     else
         n_trainingSampleIdx(find(n_trainingSampleIdx==i))=[];
         now=sprintf('NaN Iter %d SplitName %s Sample %d\n',Idx_iteration,TrainMatName,i)
     end
  end
end

N_samplereal=length(n_trainingSampleIdx);
X_temp=zeros(N_Rank,N_samplereal);
for i=1:N_samplereal
    X_temp(:,i)=X(:,n_trainingSampleIdx(i)).*N_s(n_trainingSampleIdx(i));
end
Matrix_N_E_xx=Matrix_N_E_xx+X(:,n_trainingSampleIdx)*X_temp';
Matrix_N_F_E_x=F(:,n_trainingSampleIdx)*X_temp';
Matrix_N_L_E_x=L_train(:,n_trainingSampleIdx)*X_temp';
    
save(StasticsMatName,'X','Matrix_N_E_xx','Matrix_N_F_E_x','n_trainingSampleIdx','L_train','Matrix_N_L_E_x');
