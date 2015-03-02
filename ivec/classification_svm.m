function classification_svm(trainset_id,testset_id,featType,ivec_dim,odir)
% e.g. classification_svm('mfcc_3s_train_80hrs_ubm1024_iVdim200','mfcc_3s_dev2_80hrs_ubm1024_iVdim200','mfcc',200)

addpath('./DETware_v2.1') 

train_data={[];[];[]};
train_label={[];[];[]};
comb_data={[]};
comb_label={[]};
test_data={[]};
test_label={[]};

statsdir= sprintf('%s/split_files_stats',odir);
labeldir= sprintf('%s/split_files',odir);

testset={trainset_id,testset_id};

% BUILD IVECTOR DATA SETS AND LABELS
for setID = 1 : length(testset),

    gvectorfiles = dir(sprintf('%s/%s*.gvector.mat',statsdir,testset{setID}));
    if length(gvectorfiles)==0, statsdir=statsdir2; gvectorfiles = dir(sprintf('%s/%s*.gvector.mat',statsdir2,testset{setID})); end
    for durID = 1,

        for k = 1 : length(gvectorfiles),

            datafile = fullfile(statsdir,gvectorfiles(k).name);
            [~,filename] = fileparts(datafile);

            labelfile = fullfile(labeldir,strrep(filename,'_prenormalized.gvector','.lst'));
            labelfile = strrep(labelfile,sprintf('_iVdim%i',ivec_dim),'');
            t1=load(labelfile);
            load(datafile);
            if strcmp(testset{setID},trainset_id)   
               comb_label{durID}=[comb_label{durID};t1];
               comb_data{durID}=[comb_data{durID} X];
            elseif strcmp(testset{setID},testset_id)
                if length(t1)~=size(X,2)
                   fprintf('label files for %s do not match. Rebuilding. \n', testset{setID});
                   [~,~]=unix(sprintf('bash list2lst.bash %s %s/split_files; ',testset{setID},odir));
                   t1=load(labelfile);
                   if length(t1)~=size(X,2), fprintf('error in labels files\n'); return; end;
                end
                test_label{durID}=[test_label{durID};t1];
                test_data{durID}=[test_data{durID} X];
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

% NUMBER OF LANGUAGES
nbClasses = length(unique(train_label_real));  

% ORDERING and COMPUTING NUMBER OF CLASS SAMPLES 

train_data_real_order=[];
train_label_real_order=[];
nMaxSampPerClassReal=zeros(1,nbClasses);

for i=1:nbClasses
   idx=find(train_label_real==i);
   
   t1=rand(1,length(idx));
   [t2 t3]=sort(t1);
   idx=idx(t3);
   
   clear  t1;
   clear  t2;
   clear  t3;
   
   nMaxSampPerClassReal(i)=length(idx);
   train_data_real_order=[train_data_real_order train_data_real(:,idx)];
   train_label_real_order=[train_label_real_order;train_label_real(idx)];
end

fprintf('nb training samples %i\n', size(train_data_real_order,2));
clear train_data_real;
clear train_label_real;

wccnFlag = true;

if wccnFlag
    
    for i=1:nbClasses
       idx=find(train_label_real_order==i);
       data_class_i=train_data_real_order(:,idx);
       name2=sprintf('%s/ivec_training/%s/X.isdvector.train.%s.%i.%d.normsdc',odir,trainset_id,featType,ivec_dim,i);
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
    
modeldir='model';
if ~exist(modeldir), mkdir(modeldir); end
modelfile=fullfile(modeldir,sprintf('%s.mat',trainset_id));
disp(modelfile)

if ~exist(modelfile)
   % TRAIN SVM
   param=sprintf('-s 1 -c 5 -r 1 -g 5 -w1 %f -w2 %f -w3 1 -w4 %f -w5 %f -w6 %f',nMaxSampPerClassReal(3)/nMaxSampPerClassReal(1),nMaxSampPerClassReal(3)/nMaxSampPerClassReal(2),nMaxSampPerClassReal(3)/nMaxSampPerClassReal(4),nMaxSampPerClassReal(3)/nMaxSampPerClassReal(5),nMaxSampPerClassReal(3)/nMaxSampPerClassReal(6));
   fprintf('%s \n', param);
   model=train_libpoly(train_label_real_order,sparse(train_data_real_order)',param);
   save(modelfile,'model');
else
   load(modelfile,'model');
end

% PREDICT SVM
[l,a,s]=predict_libpoly(test_label_real,sparse(test_data_real)',model,'-b 1');

scoredir='score';
if ~exist(scoredir), mkdir(scoredir); end
save(fullfile(scoredir,sprintf('%s.mat',strrep(trainset_id,'_3s_train',''))),'test_label_real','s','l');

% SCORING
Accuracy_total=mean(test_label_real==l);

accuracy=zeros(6,1);
for i=1:6
   idx=find(test_label_real==i);
   accuracy(i)=mean(l(idx)==i);
end

for i=1:size(s,1)
    tt=s(i,:);
    tt=exp(tt);
    score(i,:)=tt./sum(tt);
end

truescore=[];
falsescore=[];
for i=1:length(test_label_real)
   for j=1:5
      if( test_label_real(i)==j)
          truescore=[truescore;score(i,j)];
      else
          falsescore=[falsescore;score(i,j)];
      end
   end
end

[P_miss,P_fa] = Compute_DET (truescore, falsescore);
Pmiss_min = 0.001;
Pmiss_max = 0.2;
Pfa_min = 0.001;
Pfa_max = 0.2;
Set_DET_limits(Pmiss_min,Pmiss_max,Pfa_min,Pfa_max);

C_miss = 1;
C_fa = 1;
P_target = 0.5;
Set_DCF(C_miss,C_fa,P_target);

[DCF_opt Popt_miss Popt_fa] = Min_DCF(P_miss,P_fa);

[tt idx]=min(abs(P_miss-P_fa));
EER=P_miss(idx);

fprintf('ACC = %.4f\n',Accuracy_total);
fprintf('EER = %.4f\n',EER);
fprintf('DCF = %.4f\n',DCF_opt);
