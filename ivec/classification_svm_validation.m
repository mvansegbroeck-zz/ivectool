function classification_svm_task_loo_train(trainset_id,testset_id,featType,ivec_dim,odir)
% classification_svm5_allsamples('gtfccmcc2spchpca80plpnr_3s_train_80hrs_ubm2048_iVdim400','gtfccmcc2spchpca80plpnr_3s_dev2_80hrs_ubm2048_iVdim400','gtfccmcc2spchpca80plpnr',400)

nbClasses = 3;  % NUMBER OF LANGUAGES

substrmatch = @(x,y) ~cellfun(@isempty,strfind(y,x));
findmatching = @(x,y) substrmatch(x,y);

statsdir= sprintf('%s/split_files_stats',odir);
labeldir= sprintf('%s/split_files',odir);
svmdir='svm.validation';
if ~exist(svmdir), mkdir(svmdir); end

testset={trainset_id,testset_id};
%tasks={'readingspanSentence','stroopdualtask','strooptimepressure'};
tasks={'r'}; % global model
%tasks={'readingspanSentence'};
%tasks={'stroopdualtask'};
%tasks={'stroop'};
addpath(genpath('~/matlab/mRMR_0.9_compiled/'));

% ---speakers---
spkrs=[1 2 7 13 14 15 18 22 24 25 27 6 9 11 12 20 21 26 90:99];
spkr_norm={'no','MN','MVN'};
%spkr_norm={'MVN'};

for spkr_norm_id = 1:length(spkr_norm) 

  test_label_real_tot=[];
  l_tot=[];
  s_tot=[];
  for taskID = 1 : length(tasks),
  
    train_data={[];[];[]};
    train_label={[];[];[]};
    comb_data={[]};
    comb_label={[]};
    test_data={[]};
    test_label={[]};
    train_spkr={};
    test_spkr={};
  
    % BUILD IVECTOR DATA SETS AND LABELS
    for setID = 1 : length(testset),
    
        gvectorfiles = dir(sprintf('%s/%s*_iVdim%i*.gvector.mat',statsdir,testset{setID},ivec_dim));
        for durID = 1,
    
            for k = 1 : length(gvectorfiles),
    
                datafile = fullfile(statsdir,gvectorfiles(k).name);
                [~,filename] = fileparts(datafile);
    
                labelfile = fullfile(labeldir,strrep(filename,'_prenormalized.gvector','.lst'));
                labelfile = strrep(labelfile,sprintf('_iVdim%i',ivec_dim),'');
                t1=load(labelfile);
                load(datafile);
  
                % label and ivector extraction per task
                listfile = strrep(labelfile,'.lst','.list');
                fid=fopen(listfile); featfiles=textscan(fid,'%s\n'); featfiles=featfiles{1}; fclose(fid);
                taskndx=find(findmatching(tasks{taskID},featfiles));;
                t1=t1(taskndx);
                X=X(:,taskndx);
               
                if strcmp(testset{setID},trainset_id)   
                    if length(t1)~=size(X,2)
                       fprintf('label files for %s do not match.\n', labelfile); return;
                    end
                   comb_label{durID}=[comb_label{durID};t1];
                   comb_data{durID}=[comb_data{durID} X];
                   for spkrID=spkrs
                      train_spkr{end+1}=find(substrmatch(sprintf('subj%.2i',spkrID),featfiles));
                   end
                   train_spkr=train_spkr(~cellfun('isempty',train_spkr)); 
                elseif strcmp(testset{setID},testset_id)
                   if length(t1)~=size(X,2)
                      fprintf('label files for %s do not match. \n', testset{setID}); return;
                   end
                   test_label{durID}=[test_label{durID};t1];
                   test_data{durID}=[test_data{durID} X];
                   for spkrID=spkrs
                      test_spkr{end+1}=find(substrmatch(sprintf('subj%.2i',spkrID),featfiles));
                   end
                   test_spkr=test_spkr(~cellfun('isempty',test_spkr)); 
                end
    
                clear t1;
                clear X;
    
            end
        end
    end
    
    train_data_real = comb_data{1};
    train_label_real = comb_label{1};
    
    test_data_real=test_data{1};
    test_label_real=test_label{1};

    switch spkr_norm{spkr_norm_id}
     case { 'MN' } 
       for k=1:length(train_spkr)
          train_data_real(:,train_spkr{k})=train_data_real(:,train_spkr{k})-repmat(mean(train_data_real(:,train_spkr{k}),2),1,length(train_spkr{k}));
       end
       for k=1:length(test_spkr)
          test_data_real(:,test_spkr{k})=test_data_real(:,test_spkr{k})-repmat(mean(test_data_real(:,test_spkr{k}),2),1,length(test_spkr{k}));
       end
     case { 'MVN' } 
       for k=1:length(train_spkr)
          train_data_real(:,train_spkr{k})=mvn(train_data_real(:,train_spkr{k}));
       end
       for k=1:length(test_spkr)
          test_data_real(:,test_spkr{k})=mvn(test_data_real(:,test_spkr{k}));
       end
    end

    % REMOVE OUTLIERS (CORRUPTED FEATURES), NORMALIZE WITH 2NORM
    
    idx=find(sum(isnan(train_data_real))>0);
    train_data_real(:,idx)=[];
    train_label_real(idx)=[];
    idx=find(sum(isnan(test_data_real))>0);
    if(length(idx)>0)
    test_data_real(:,idx)=test_data_real(:,1); % OUTLIERS TEST ARE JUST SET TO THE FIRST IVECTOR
    %test_label_real(idx)=[];
    end
    for i=1:length(train_label_real)
        train_data_real(:,i)=train_data_real(:,i)./norm(train_data_real(:,i),2); 
    end
    
    for i=1:length(test_label_real)
        test_data_real(:,i)=test_data_real(:,i)./norm(test_data_real(:,i),2); 
    end
    
    nMaxSampPerClassReal=zeros(nbClasses,1);
    
    train_data_real_order=[];
    train_label_real_order=[];
    
    for i=1:nbClasses
       idx=find(train_label_real==i);
       
       t1=rand(1,length(idx));
       [t2 t3]=sort(t1);
       idx=idx(t3);
       
       clear  t1;
       clear  t2;
       clear  t3;
       
       nMaxSampPerClassReal(i)=length(idx);
       idx=idx(1:nMaxSampPerClassReal(i));%only 500 samples
       train_data_real_order=[train_data_real_order train_data_real(:,idx)];
       train_label_real_order=[train_label_real_order;train_label_real(idx)];
    end
    
    fprintf('nb training samples %i\n', size(train_data_real_order,2));
    clear train_data_real;
    clear train_label_real;
    
    featsel=false;
    if featsel == true
      featsel_frac=0.9;
      fprintf('selected comp %i \n',fix(featsel_frac*ivec_dim));
      train_data_real_order=fix(train_data_real_order*1000);
      test_data_real=fix(test_data_real*1000);
      %featselndx = [1:ivec_dim]; 
      featselndx = mrmr_mid_d(train_data_real_order', train_label_real_order, fix(featsel_frac*ivec_dim));
      %train_data_real_order(featselndx,:)=[];
      %test_data_real(featselndx,:)=[];
      %train_data_real_order=train_data_real_order./1000;
      %test_data_real=test_data_real./1000;
    
      train_data_real_order=train_data_real_order(featselndx,:)./1000;
      test_data_real=test_data_real(featselndx,:)./1000;
    end

    wccnFlag = false;
    
    if wccnFlag
        ivec_dim=size(train_data_real_order,1);
        
        for i=1:nbClasses
           idx=find(train_label_real_order==i);
           data_class_i=train_data_real_order(:,idx);
           ivecdir=sprintf('%s/ivec_training/%s',odir,trainset_id);
           if ~exist(ivecdir), mkdir(ivecdir); end
           name2=sprintf('%s/X.isdvector.train.%s.%i.%d.normsdc',ivecdir,featType,ivec_dim,i);
           fid=fopen(name2,'wb');
           fwrite(fid,size(data_class_i,2),'int32');
           fwrite(fid,size(data_class_i,1)+1,'int32');
           fwrite(fid,1,'int16');
           fwrite(fid,0,'int16');
           data_class_i=[ones(1,size(data_class_i,2));data_class_i];
           data_class_i=reshape(data_class_i,(size(data_class_i,1))*size(data_class_i,2),1);
           fwrite(fid,data_class_i,'float32');
           fclose(fid);
    
           clear data_class_i;
        end
        unix(sprintf('ls %s/ivec_training/%s/X.isdvector.train.%s.%i.*.normsdc > %s/ivec_training/%s/X.isdvector.train.%s.%i.normsdc.list',odir,trainset_id,featType,ivec_dim,odir,trainset_id,featType,ivec_dim));
        
        % DO WCCN
        unix(sprintf('./wccn/WCCN %s/ivec_training/%s/X.isdvector.train.%s.%i.normsdc.list %s/ivec_training/%s/train.%s.normsdc.wccn.txt %i',odir,trainset_id,featType,ivec_dim,odir,trainset_id,featType,ivec_dim));
        
        WCCNMatrix=load(sprintf('%s/ivec_training/%s/train.%s.normsdc.wccn.txt',odir,trainset_id,featType));
        train_data_real_order=WCCNMatrix*train_data_real_order;
        test_data_real=WCCNMatrix*test_data_real;
    
        for i=1:length(train_label_real_order)
            train_data_real_order(:,i)=train_data_real_order(:,i)./norm(train_data_real_order(:,i),2); 
        end
    
        for i=1:length(test_label_real)
            test_data_real(:,i)=test_data_real(:,i)./norm(test_data_real(:,i),2); 
        end
        
        
    end

    %pca_dim=300;
    %[~,pca_mat]=pca(train_data_real_order',pca_dim); 
    %train_data_real_order=pca_mat.M'*train_data_real_order;
    %test_data_real=pca_mat.M'*test_data_real;

    trainls=fullfile(svmdir,sprintf('%s_iVdim%i_%s_spkrNorm%s.ls',trainset_id,ivec_dim,tasks{taskID},spkr_norm{spkr_norm_id}));
    testls=fullfile(svmdir,sprintf('%s_iVdim%i_%s_spkrNorm%s.ls',testset_id,ivec_dim,tasks{taskID},spkr_norm{spkr_norm_id}));
    libsvmwrite(trainls, train_label_real_order, sparse(train_data_real_order'));
    libsvmwrite(testls, test_label_real, sparse(test_data_real'));
    
    tmp=strfind(testset_id,'_');
    
    fprintf('start SVM training - spkr normalisation: %s \n', spkr_norm{spkr_norm_id});
   
    for svm_order = [3 5 7] 
    
     % TRAIN SVM
     miss_class_cost=1;
     param=sprintf('-t 1 -d %i -r 1 -g 1 -b 1 -c %.3f -m 20 -w1 %f -w2 %f -w3 1',svm_order, miss_class_cost, nMaxSampPerClassReal(3)/nMaxSampPerClassReal(1),nMaxSampPerClassReal(3)/nMaxSampPerClassReal(2));
     %fprintf('%s \n', param);
     trainmdl=fullfile(svmdir,sprintf('%s_iVdim%i_svm%i_%s_spkrNorm%s.mdl',trainset_id,ivec_dim,svm_order,tasks{taskID},spkr_norm{spkr_norm_id}));
     testscore=fullfile(svmdir,sprintf('%s_iVdim%i_svm%i_%s_spkrNorm%s.score',testset_id,ivec_dim,svm_order,tasks{taskID},spkr_norm{spkr_norm_id}));
     [~,~]=  unix(sprintf('./svm-train-mp %s %s %s',param,trainls,trainmdl));
     
     % PREDICT SVM
     [~,~]=unix(sprintf('./svm-predict-mp -b 1 %s %s %s',testls,trainmdl,testscore));
     
     % SCORING
     [l,s1,s2,s3] = textread(testscore, '%s%s%s%s');
     s=zeros(length(s1)-1,nbClasses);
     l=str2double(l(2:end));
     s(:,1)=str2double(s1(2:end)); s(:,2)=str2double(s2(2:end)); s(:,3)=str2double(s3(2:end));
         
     Accuracy_total=mean(test_label_real==l);
     
     accuracy=zeros(nbClasses,1);
     for i=1:nbClasses
        idx=find(test_label_real==i);
        accuracy(i)=mean(l(idx)==i);
     end
     
     [WeightedPercentage UnweightedPercentage ConfusionMatrix] = ComputeUnweightedPercentageAndConfusionMatrix(test_label_real,l); 
     fprintf('norm %s - svm %i: %.3f %.3f \t',spkr_norm{spkr_norm_id}, svm_order,WeightedPercentage,UnweightedPercentage);
     fprintf('%i %i %i\n', sum(l==1), sum(l==2), sum(l==3))

     fid=fopen(sprintf('%s.labels',testscore),'w');
     for kk = 1 : length(featfiles),
         [~,filename,~]=fileparts(featfiles{kk});
         [~,filename,~]=fileparts(filename);
         fprintf(fid,'%s.wav,L%i\n', filename, l(kk));
     end
     fclose(fid);
     
     test_label_real_tot=[test_label_real_tot;test_label_real];
     l_tot=[l_tot;l];
     s_tot=[s_tot;s];
 
   end
 
  end
end
exit;
