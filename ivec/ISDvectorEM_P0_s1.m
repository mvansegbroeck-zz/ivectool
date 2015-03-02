% P0. INITIALIZE PARAMETERS FOR EM TRAINING
function ISDvectorEM_P0_s1(trainingdir, N_Rank, nMixNum, featType, nVecSize, nbClassesTrain)

%clear all;
%N_Rank=400; 
% 600(for RATS) = dimension of i-vector, related to the #Gaussians, in the GMM and the amount of the data
%nMixNum=128; 
% 2048 (for RATS)

Idx_iteration=0;
n_MeansuperDim=sum(nMixNum*nVecSize);

T_init=single(randn(n_MeansuperDim,N_Rank));

for i=1:N_Rank
    T_init(:,i)=(T_init(:,i)./norm(T_init(:,i),2)).*0.07;
end

W_init=single(randn(nbClassesTrain,N_Rank));

for i=1:N_Rank
    W_init(:,i)=(W_init(:,i)./norm(W_init(:,i),2)).*0.07;
end

InvSigma1=ones(sum(nVecSize),nMixNum);
InvSigma2=ones(nbClassesTrain,1);

matrix1=single(zeros(N_Rank,N_Rank));   
matrix1=T_init'*T_init;    
matrix2=single(zeros(N_Rank,N_Rank));
matrix2=W_init'*W_init;

N_s_table=[301:30000];
log_N_s_table=log(N_s_table);
log_N_s_table_length=300;
log_N_s_table_interval=(max(log_N_s_table)-min(log_N_s_table))/(log_N_s_table_length-1);

%ProjectMatrixPool=zeros(N_Rank,N_Rank,log_N_s_table_length,'single');
ProjectMatrixPool=cell(log_N_s_table_length,1);

for i=1:log_N_s_table_length
   N_s_app= exp(min(log_N_s_table)+log_N_s_table_interval*(i-1));
   ProjectMatrixPool{i}=inv(eye(N_Rank)+matrix1.*N_s_app+matrix2.*N_s_app);
   %ProjectMatrixPool{i}=inv(eye(N_Rank)+matrix1.*N_s_app);
end
    
save(fullfile(trainingdir,sprintf('T_init_isdvecfa_%s_%i_p0_s1_it0.mat',featType,N_Rank)),'-v7.3');


