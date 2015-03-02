function stats2mat_feat(featType,ubmfile,statsfile,odir)

% stats
statsdir=strcat(odir,'/split_files_stats/'); 

matfile = strrep(statsfile,'.stats','_prenormalized.mat');
fprintf('processing: %s \n', matfile);

fid=fopen(ubmfile,'rb');
nMixNum=fread(fid,1,'int32');
nVecSize=fread(fid,1,'int32');
nVecSizeStore=fread(fid,1,'int32');
nFeatKind=fread(fid,1,'int32');
nTotalParamSize=fread(fid,1,'int32');
fseek(fid,512,-1);
pfWeightBuf=fread(fid,nMixNum,'float32');
pfMeanBuf=fread(fid,nMixNum*nVecSizeStore,'float32');
pfDiagCovBufInv=fread(fid,nMixNum*nVecSizeStore,'float32');
pfDiagCovBufInv=reshape(pfDiagCovBufInv,nVecSizeStore,nMixNum);
pfMatBuf=fread(fid,nMixNum,'float32');
fclose(fid);

pfDiagCovBufInv=pfDiagCovBufInv(1:nVecSize,:);
pfMeanBuf=reshape(pfMeanBuf,nVecSizeStore,nMixNum);
pfMeanBuf=pfMeanBuf(1:nVecSize,:);
pfMeanBuf=reshape(pfMeanBuf,nVecSize*nMixNum,1);

fid=fopen(statsfile,'rb');
nMixNum=fread(fid,1,'int32');
nVecSize=fread(fid,1,'int32');
nUttNum=fread(fid,1,'int32');
flag1=fread(fid,1,'int32');
N=fread(fid,nMixNum*nUttNum,'float32');
flag2=fread(fid,1,'int32');
F=fread(fid,nMixNum*nVecSize*nUttNum,'float32');
flag3=fread(fid,1,'int32');
UtterFrameNum=fread(fid,nUttNum,'int32');
flag4=fread(fid,1,'int32');
UtterNames=fread(fid,32*nUttNum,'char');
fclose(fid);

N=reshape(N,nMixNum,nUttNum);
F=reshape(F,nMixNum*nVecSize,nUttNum);
assert(size(N,1)==nMixNum);

N=N+1e-20;


for i=1:nMixNum
    F((i-1)*nVecSize+1:i*nVecSize,:)=F((i-1)*nVecSize+1:i*nVecSize,:)./(repmat(N(i,:),nVecSize,1));
end
F=single(F-repmat(pfMeanBuf,1,nUttNum));



for k=1:size(F,2)
   W=(sqrt(pfDiagCovBufInv).*[sqrt(N(:,k)/sum(N(:,k)))*ones(1,nVecSize)]');
   W=W(:);
   %W=zeros(nMixNum*nVecSize,1);
   %for i=1:nMixNum
   %     W((i-1)*nVecSize+1:i*nVecSize)=sqrt(pfDiagCovBufInv(:,i)).*sqrt(N(i,k)/sum(N(:,k)));      
   %end

   if(sum(isnan(W))>0)
      now=sprintf('W isnan error! k=%d',k)
   end
   F(:,k)= F(:,k).*W;
end

save(matfile,'F','N','-v7.3');

exit;
