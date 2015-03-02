function make_gfcc(wav_list)
% Generate MFCC features in HTK format (*.mfc)
% for the files given in the database file
% wav_list: list of wav files to be processed

 if ~exist('wav_list', 'var')
     fprintf(1, 'Input file list is missing!\n');
     exit;
 end
 
 fd=fopen(wav_list,'r');
 list_files=textscan(fd,'%s\n');
 fclose(fd);
 
 for fileID = 1:length(list_files{1})
  
     entry=list_files{1}{fileID};
     [filedir,filename,ext]=fileparts(entry);
     destdir=strrep(filedir,'/wavs','/feats'); if ~exist(destdir); mkdir(destdir); end
     gfccfile=fullfile(destdir,strcat(filename,'.gfcc'));
     if exist(strcat(gfccfile,'.mfc')),
        fprintf('features for %s already exist. Skipping...',filename); fprintf('\n');
        continue
     end
 
     % Training files are stored in the WAV format
     [sig, fs] = wavread(entry);
     
     if length(sig)==0, continue; end
 
     fprintf('file read %s',filename); fprintf('\n');
    
     % extract GFCC features + first and second order delta's 
     gfcc=FE_GFCC(sig, fs);

     % apply mean variance normalization
     gfcc = mvn(gfcc);  
 
     writeHTK(gfccfile, gfcc);  % gtfcc mvn sva
     [~,~]=unix(sprintf('./Feat2Mfc %s %s.mfc',gfccfile,gfccfile));
     [~,~]=unix(sprintf('\\rm %s',gfccfile));
 
     clear gfcc sig
         
 end

end

function [gfcc,gf,logE]=FE_GFCC(sam, fs, FrameLen, numChannel, numGFCCs )
% Generate gammatone features (GF) and gammatone frequency cepstral
% coefficients (GFCC).

 if ~exist('sampFreq', 'var')
     sampFreq = 16000;
 end
 
 if ~exist('numChannel', 'var')
     numChannel = 64;
 end
 
 if ~exist('numGFCCs', 'var')
     numGFCCs = 23;
 end
 
 if ~exist('FrameLen', 'var')
     FrameLen = 0.025*fs;
 end
 
 [M,Mstatic]=vec2featmat2(numChannel,numGFCCs-1);
 [gf,logE]=FE_GT(sam,fs,numChannel,FrameLen);
 gfaug=stacklp(gf,4);
 gfcc=M*gfaug;

end

function [gf,logE,X]=FE_GT(sam,fs,NbCh,FrameLen,tGtMat)
% make Gammatone Frequency Filtered Spectra

 Efloor=exp(-50);
 Gfloor=exp(-50);
 PreEmp=0.97;
 
 if ~exist('fs', 'var')
    fs=16000; % sampling frequency
 end
 if ~exist('FrameLen', 'var')
    FrameLen=0.025*fs;
 end
 FrameShift=0.01*fs;
 Nfft=2^ceil(log(FrameLen)/log(2));
 
 if ~exist('NbCh', 'var')
     NbCh=64; % number of GT filter banks
 end
 if ~exist('tGtMat', 'var')
     tGtMat=gammatone_matrix(Nfft,fs,NbCh)';
 end
 % truncate as in ReadWave
 NbFr=floor( (length(sam)-FrameLen+FrameShift)/FrameShift);
 sam=sam(1:NbFr*FrameShift+FrameLen-FrameShift);
 
 % DC removal
 sam=filter([1 -1],[1 -0.999],[0 reshape(sam,1,length(sam))]);
 
 % framing
 ind1=1:FrameShift:length(sam)-1-FrameLen+FrameShift;
 ind2=(1:FrameLen)';
 x=sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr)));
 
 % logE
 logE=log(max(sum(x.^2,1),Efloor));
 
 % preemphasis & windowing
 T=length(sam);
 sam=[0 sam(2:T)-PreEmp*sam(1:T-1)];
 win=hamming(FrameLen);
 % Fourier Transform
 X=fft(win(:,ones(1,NbFr)).*sam(ind1(ones(FrameLen,1),:)+ind2(:,ones(1,NbFr))),Nfft);
 X(1,:)=0;
 X=abs(X(1:Nfft/2,:));
 gf=max(tGtMat'*X,Gfloor).^(1/3);

end

function [wts,cfreqs] = gammatone_matrix(nfft, sr, nfilts, width, minfreq, maxfreq, maxlen)
% [wts,cfreqa] = fft2gammatonemx(nfft, sr, nfilts, width, minfreq, maxfreq, maxlen)
%      Generate a matrix of weights to combine FFT bins into
%      Gammatone bins.  nfft defines the source FFT size at
%      sampling rate sr.  Optional nfilts specifies the number of
%      output bands required (default 64), and width is the
%      constant width of each band in Bark (default 1).
%      minfreq, maxfreq specify range covered in Hz (100, sr/2).
%      While wts has nfft columns, the second half are all zero.
%      Hence, aud spectrum is
%      fft2gammatonemx(nfft,sr)*abs(fft(xincols,nfft));
%      maxlen truncates the rows to this many bins.
%      cfreqs returns the actual center frequencies of each
%      gammatone band in Hz.
%
% 2004-09-05  Dan Ellis dpwe@ee.columbia.edu  based on rastamat/audspec.m
% Last updated: $Date: 2009/02/22 02:29:25 $

 if nargin < 2;    sr = 16000; end
 if nargin < 3;    nfilts = 64; end
 if nargin < 4;    width = 0.5; end
 if nargin < 5;    minfreq = 50; end
 if nargin < 6;    maxfreq = sr/2; end
 if nargin < 7;    maxlen = nfft/2; end
 
 wts = zeros(nfilts, nfft);
 
 % after Slaney's MakeERBFilters
 EarQ = 9.26449;
 minBW = 24.7;
 order = 1;
 
 cfreqs = -(EarQ*minBW) + exp((1:nfilts)'*(-log(maxfreq + EarQ*minBW) + ...
                 log(minfreq + EarQ*minBW))/nfilts) * (maxfreq + EarQ*minBW);
 cfreqs = flipud(cfreqs);
 
 GTord = 4;
 
 ucirc = exp(j*2*pi*[0:(nfft/2)]/nfft);
 
 justpoles = 0;
 
 for i = 1:nfilts
   cf = cfreqs(i);
   ERB = width*((cf/EarQ).^order + minBW^order).^(1/order);
   B = 1.019*2*pi*ERB;
   r = exp(-B/sr);
   theta = 2*pi*cf/sr;
   pole = r*exp(j*theta);
 
   if justpoles == 1
     % point on unit circle of maximum gain, from differentiating magnitude
     cosomegamax = (1+r*r)/(2*r)*cos(theta);
     if abs(cosomegamax) > 1
       if theta < pi/2;  omegamax = 0;
       else              omegamax = pi;   end
     else
       omegamax = acos(cosomegamax);
     end
     center = exp(j*omegamax);
     gain = abs((pole-center).*(pole'-center)).^GTord;
     wts(i,1:(nfft/2+1)) = gain * (abs((pole-ucirc).*(pole'- ...
                                                      ucirc)).^-GTord);
   else
     % poles and zeros, following Malcolm's MakeERBFilter
     T = 1/sr;
     A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2* ...
                                                       cf*pi*T)./exp(B*T))/2;
     A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2* ...
                                                       cf*pi*T)./exp(B*T))/2;
     A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2* ...
                                                       cf*pi*T)./exp(B*T))/2;
     A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2* ...
                                                       cf*pi*T)./exp(B*T))/2;
     zros = -[A11 A12 A13 A14]/T;
 
     r100=2*exp(4*j*cf*pi*T);
     r101=-T*r100;
     r102=exp(B*T);
     r11=2*exp(-(B*T)).*exp(2*j*cf*pi*T).*T;
     r12=cos(2*cf*pi*T);
     r13=sqrt(3 - 2^(3/2))* sin(2*cf*pi*T);
     r14=sqrt(3 + 2^(3/2))* sin(2*cf*pi*T);
 
     %gain(i) = abs( (r101 + r11.*(r12 - r13)) .*  (r101 + r11.*(r12 + r13))  .* ...
     %               (r101 + r11.*(r12 - r14)) .*  (r101 + r11.*(r12 + r14))  ./ ...
     %               (-2*(r102*r102)^(-1) - r100 + 2*(r102)^(-1)+r100*(r102)^(-1)).^4);
     gain(i) =  abs((-2*exp(4*j*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                 (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                  sin(2*cf*pi*T))) .* ...
                (-2*exp(4*j*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                 (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
                  sin(2*cf*pi*T))).* ...
                (-2*exp(4*j*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                 (cos(2*cf*pi*T) - ...
                  sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
                (-2*exp(4*j*cf*pi*T)*T + 2*exp(-(B*T) + 2*j*cf*pi*T).*T.* ...
                 (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
                (-2 ./ exp(2*B*T) - 2*exp(4*j*cf*pi*T) +  ...
                 2*(1 + exp(4*j*cf*pi*T))./exp(B*T)).^4);
     wts(i,1:(nfft/2+1)) = ((T^4)/gain(i)) ...
         * abs(ucirc-zros(1)).*abs(ucirc-zros(2))...
         .*abs(ucirc-zros(3)).*abs(ucirc-zros(4))...
         .*(abs((pole-ucirc).*(pole'-ucirc)).^-GTord);
   end
 end
 
 wts = wts(:,1:maxlen);

end

function [M,Mdct]=vec2featmat2(D,cep_ord)
% Make DCT transform
 
 b0=[0 0 0 0 1  0  0  0 0];
 b1=[0 0 2 1 0 -1 -2 0 0]/10;
 b2=conv(b1(3:7),[2 1 0 -1 -2]/10);

 Mdct=[cos(pi*(1:cep_ord)'*((1:D)-0.5)/D);ones(1,D)];

 M=kron([b0;-b1;b2],Mdct);

end

function y=stacklp(x,N)
% Stack features
 if nargin<2,
   N=5;
 end
 
 [D,T]=size(x);
 xext=x(:,[ones(1,N) 1:T T*ones(1,N)]);
 
 y=zeros((2*N+1)*D,T);
 if islogical(x), y=logical(y);end % copy logical property
 for k=0:2*N,
   y(k*D+(1:D),:)=xext(:,k+1:k+T);
 end

end
