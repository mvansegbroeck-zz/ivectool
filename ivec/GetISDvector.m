function GetISDvector(trainset_id,datasets,featType,ivec_dim,NbIter, odir)

% 3 seconds full training set: GetISDvector('mfcc_3s_train_addlen','mfcc_3s_train_addlen',10)

trainingdir = fullfile(odir,'ivec_training',trainset_id);
statsdir= sprintf('%s/split_files_stats',odir);

j=NbIter; % final iteration of i-vector model
ivec_modelfile=fullfile(trainingdir,sprintf('T_init_isdvecfa_%s_%i_p0_s1_it%d.mat',featType,ivec_dim,NbIter));
fprintf('i-vector model: %s\n',ivec_modelfile);

% BUILD IVECTOR DATA SETS AND LABELS
for setID = 1 : length(datasets),
    train_mat_names = dir(sprintf('%s/%s*prenormalized.mat',statsdir,datasets{setID}));
    for k = 1 : length(train_mat_names),
        
        train_mat_name = fullfile(statsdir,train_mat_names(k).name);
        X_mat_name = strrep(train_mat_name,'_prenormalized.mat',sprintf('_iVdim%i_prenormalized.gvector.mat',ivec_dim));
    
        fprintf('extracting i-vectors for %s (split %i) \n',datasets{setID},k);
        GetISDvector_single(train_mat_name,X_mat_name,ivec_modelfile);
    end
end
exit;
