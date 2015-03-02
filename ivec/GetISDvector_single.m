function GetISDvector_single(TrainMatName,XMatName,ModelName)

load(ModelName);
load(TrainMatName);

N_s=sum(N);

n_TrainSample=size(F,2);
n_MeansuperDim=size(F,1);

X=zeros(N_Rank,n_TrainSample);
InvSigmaVector=reshape(InvSigma1,sum(nVecSize*nMixNum),1);

for i=1:n_TrainSample
    %ProjMatrix=inv(eye(N_Rank)+matrix1.*N_s(i));
    %t1=InvSigmaVector.*(F(:,i).*N_s(i));
    %X(:,i)=ProjMatrix*(T_init'*t1);
    ProjMatrix=(eye(N_Rank)+matrix1.*N_s(i));
    t1=InvSigmaVector.*(F(:,i).*N_s(i));
    X(:,i)=ProjMatrix\(T_init'*t1);
    clear t1;
end

save(XMatName,'X');
