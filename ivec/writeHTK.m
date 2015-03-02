function writeHTK(outfile,e)
% Write data into HTK format
% outfile: file name
% e: data matrix
% Written by Yang Shao, and adapted by Xiaojia Zhao in Oct'11

dtfile=fopen(outfile,'wb');

% header
fwrite(dtfile,size(e,2),'int32');       % nSamples
fwrite(dtfile,size(e,1),'int32');       % sampSize

% write data
ee=reshape(e,size(e,1)*size(e,2),1);

%for i=1:size(e,2)
%for j=1:size(e,1)    
%    fwrite(dtfile,e(j,i),'real*4');
%end
%end

fwrite(dtfile,e,'float32');

fclose(dtfile);