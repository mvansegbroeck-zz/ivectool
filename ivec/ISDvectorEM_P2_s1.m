% P2. M STEP OF EM TRAINING, ACCUMULATION OVER SPLIT FILES

function ISDvectorEM_P2_s1(LastModelMatName,trainingset,nb_splits,NewModelMatName,ivec_dim,odir)


load (LastModelMatName);

ProjectMatrixPool_=cell(size(ProjectMatrixPool,3));
for i=1:size(ProjectMatrixPool,3)
   ProjectMatrixPool_{i}=ProjectMatrixPool(:,:,i);
end
ProjectMatrixPool=ProjectMatrixPool_;

statsdir=sprintf('%s/split_files_stats/',odir);
N_StaticsMat=nb_splits;
N_TrainMat=nb_splits;

Matrix_N_E_xx_total=zeros(N_Rank,N_Rank);
Matrix_N_F_E_x_total=zeros(sum(nVecSize*nMixNum),N_Rank);
Matrix_N_L_E_x_total=zeros(nbClassesTrain,N_Rank);

VectorEnergyError1=zeros(n_MeansuperDim,1);
VectorEnergyError2=zeros(nbClassesTrain,1);


min_log_N_s_table=min(log_N_s_table);

for i=1:N_StaticsMat
   fprintf('part 2 sec 1 Iter %d Split %d\n',Idx_iteration,i);
   load(fullfile(statsdir,sprintf('%s_split_%i_iVdim%i_prenormalized_tmpstats.mat',trainingset,i,ivec_dim)));
   Matrix_N_E_xx_total=Matrix_N_E_xx_total+Matrix_N_E_xx;
   Matrix_N_F_E_x_total=Matrix_N_F_E_x_total+Matrix_N_F_E_x;
   Matrix_N_L_E_x_total=Matrix_N_L_E_x_total+Matrix_N_L_E_x;
end
T_init=Matrix_N_F_E_x_total*inv(Matrix_N_E_xx_total);
W_init=Matrix_N_L_E_x_total*inv(Matrix_N_E_xx_total);

N_trainingsampleTotal=0;
for ii=1:N_TrainMat
    load(fullfile(statsdir,sprintf('%s_split_%i_iVdim%i_prenormalized_tmpstats.mat',trainingset,ii,ivec_dim)));
    load(fullfile(statsdir,sprintf('%s_split_%i_prenormalized.mat',trainingset,ii)));
    fprintf('part 2 sec 2 Iter %d Split %d\n',Idx_iteration,ii);

    N_s=sum(N);
    
     for jj=1:length(N_s)
        
        ProjMatrixPoolIdx=round((log(N_s(jj))-min_log_N_s_table)/log_N_s_table_interval);
        ProjMatrixPoolIdx=min(ProjMatrixPoolIdx,log_N_s_table_length);
        ProjMatrixPoolIdx=max(ProjMatrixPoolIdx,1); 
        N_s(jj)=exp(min(log_N_s_table)+log_N_s_table_interval*(ProjMatrixPoolIdx-1));
     end
    
    n_TrainSample=size(F,2);
    F_recon=T_init*X;
    L_recon=W_init*X;

    for k=1:length(n_trainingSampleIdx)
    i=n_trainingSampleIdx(k);
     t6=(F(:,i)-F_recon(:,i)).*F(:,i);
     VectorEnergyError1=VectorEnergyError1+t6.*N_s(i);
     t6=(L_train(:,i)-L_recon(:,i)).*L_train(:,i);
     VectorEnergyError2=VectorEnergyError2+t6.*N_s(i);
    end
    N_trainingsampleTotal=N_trainingsampleTotal+length(n_trainingSampleIdx);
end
fprintf('part 2 sec 3 Iter %d\n',Idx_iteration);
VectorEnergyError1=max(VectorEnergyError1,1);
InvSigma1=N_trainingsampleTotal./(VectorEnergyError1);
InvSigma1=reshape(InvSigma1,sum(nVecSize),nMixNum);
InvSigma2=N_trainingsampleTotal./(VectorEnergyError2);

Idx_iteration=Idx_iteration+1;   
matrix1=zeros(N_Rank,N_Rank);
for i=1:nMixNum
    matrix1=matrix1+T_init((i-1)*sum(nVecSize)+1:i*sum(nVecSize),:)'*diag(InvSigma1(:,i))*T_init((i-1)*sum(nVecSize)+1:i*sum(nVecSize),:);    
end
matrix2=zeros(N_Rank,N_Rank);
matrix2=W_init'*diag(InvSigma2)*W_init;    

%fprintf('%.3f %.3f\n',min(min(VectorEnergyError1)), max(max(VectorEnergyError1)));
%fprintf('%.3f %.3f %.3f %.3f\n',min(min(InvSigma1)),max(max(InvSigma1)),min(min(T_init)),max(max(T_init)));
%fprintf('%.3f %.3f %.3f %.3f\n',min(min(matrix1)),max(max(matrix1)),min(min(matrix2)),max(max(matrix2)));

for i=1:log_N_s_table_length
   N_s_app= exp(min(log_N_s_table)+log_N_s_table_interval*(i-1));
   ProjectMatrixPool{i}=inv(eye(N_Rank)+matrix1.*N_s_app+matrix2.*N_s_app);
   %ProjectMatrixPool{i}=inv(eye(N_Rank)+matrix1.*N_s_app);
end
save (NewModelMatName,'N_Rank','Idx_iteration','nbClassesTrain','n_MeansuperDim','T_init','InvSigma1','matrix1','nMixNum','nVecSize','ProjectMatrixPool','log_N_s_table','log_N_s_table_length','log_N_s_table_interval','W_init','InvSigma2','matrix2');

