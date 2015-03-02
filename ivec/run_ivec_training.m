function run_ivec_training(featType,trainset_id,NbIter,ivec_dim,ubmSize,odir)

% 3 seconds full training set: run_ivec_training('mfcc',10,600,2048)

switch featType
 case {'gfcc'}
   feat_dim = 69;
 otherwise
   disp('Unknown feature type.'); 
   return
end
trainingdir = fullfile(odir,'ivec_training',trainset_id);
if ~exist(trainingdir), mkdir(trainingdir); end

% define splits
splits=dir(sprintf('%s/split_files/%s_split_*.list',odir,trainset_id));
nb_splits=length(splits);

% define nb of labels
lbls=load(sprintf('%s/split_files/%s.lst',odir,trainset_id));
nbClassesTrain=length(unique(lbls));

% initialize
fprintf('Initialize EM training \n');
ISDvectorEM_P0_s1(trainingdir,ivec_dim,ubmSize,featType, feat_dim, nbClassesTrain);

% EM training
fprintf('EM training \n');

for j=1:NbIter
   lastmodelfile=fullfile(trainingdir,sprintf('T_init_isdvecfa_%s_%i_p0_s1_it%d.mat',featType,ivec_dim,j-1));
   load(lastmodelfile);
   for i=1:nb_splits
       fprintf('Step E: Iter %d Split %d\n',j,i);
       train_mat_name=sprintf('%s/split_files_stats/%s_split_%d_prenormalized.mat',odir,trainset_id,i);
       train_lst_name=sprintf('%s/split_files/%s_split_%d.lst',odir,trainset_id,i);
       stats_mat_name=sprintf('%s/split_files_stats/%s_split_%d_iVdim%i_prenormalized_tmpstats.mat',odir,trainset_id,i,ivec_dim);
       ISDvectorEM_P1_s1_f2(train_mat_name,train_lst_name,stats_mat_name,N_Rank,Idx_iteration,nbClassesTrain,n_MeansuperDim,T_init,InvSigma1,matrix1,nMixNum,nVecSize,ProjectMatrixPool,log_N_s_table,log_N_s_table_length,log_N_s_table_interval,W_init,InvSigma2,matrix2);
   end
   
   fprintf('Step M: Iter %d \n',j);
   newmodelfile=fullfile(trainingdir,sprintf('T_init_isdvecfa_%s_%i_p0_s1_it%d.mat',featType,ivec_dim,j));
   ISDvectorEM_P2_s1(lastmodelfile,trainset_id,nb_splits,newmodelfile,ivec_dim,odir);
   unix(sprintf('\\rm %s',lastmodelfile));
end

% clean up
for i=1:nb_splits
   stats_mat_name=sprintf('%s/split_files_stats/%s_split_%d_iVdim%i_prenormalized_tmpstats.mat',odir,trainset_id,i,ivec_dim);
   unix(sprintf('\\rm %s',stats_mat_name));
end
exit;
