#!/bin/bash

featType=$1
featdim=$2
ubmsize=$3
wavlist=$4

trainset='train'

cd ubm

mkdir -p cfg/ models/ out/ feature_list/
cfgfile='ubm.'${ubmsize}'.'${featType}'.'${trainset}'.cfg'

featlist=`basename $wavlist`
featlist=${featlist%.*}.txt
echo $featlist
cp -f $wavlist feature_list/$featlist
sed -i 's/\.wav/\.'$featType'\.mfc/' feature_list/$featlist
sed -i 's/wavs\//feats\//' feature_list/$featlist

ubmdir=ubm${ubmsize}.${featType}.${trainset}
cp TEMPLATE.cfg cfg/$cfgfile

sed -i 's/featlist/'${featlist}'/g' cfg/$cfgfile
sed -i 's/ubmsize/'$ubmsize'/g' cfg/$cfgfile
sed -i 's/featdimension/'$featdim'/g' cfg/$cfgfile
sed -i 's/ubmdir/'${ubmdir}'/g' cfg/$cfgfile

mkdir -p models/$ubmdir/UBM_1

./GMM_UBM_train_linux_em64t_MultiThread16 cfg/$cfgfile > out/log.${cfgfile%.*} 

# convert ubm file to new format
ubm_file_old=models/$ubmdir/UBM_$ubmsize/gmm_train12.param
ubm_file=models/$ubmdir.param
./oldUBM2newUBM $ubm_file_old $ubm_file &
wait
  
# convert ubm file to txt
./ReadGmmModel -i $ubm_file -o ${ubm_file}.txt > out/ReadGmmModel_${ubmdir} &
wait
echo "ubm conversion finished"

\rm -rf models/$ubmdir
echo "ubm dirs cleaned up"
