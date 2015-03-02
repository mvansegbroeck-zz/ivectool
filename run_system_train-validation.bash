#!/bin/bash

# process flags
initialize=true
feat_extract=true
ubm_train=true
ubm_stats=true
prenorm_stats=true
ivec_train=true
ivec_extract=true
ivec_class=true

# define parameters
feat_type=gfcc
feat_dim=69;
# --
nb_splits=4
ubm_size=1024
# --
ivec_iter=10
ivec_dim=400

# db files for train and test
dbfile_train=train_wavs_12k.db
dbfile_test=all_testset_10s.db

# path setting 
wav_dir=$PWD/wavs/
feat_dir=$PWD/feats/
dbfile_dir=$PWD/dbfile/
odir=$PWD/ivec_train-validation/
classdir=$odir/classification

# useful shortcuts
train_set_id=${feat_type}_ubm${ubm_size}_training
val_set_id=${feat_type}_ubm${ubm_size}_validation
wav_list_train=$odir/split_files/${train_set_id}.wavlist
wav_list_val=$odir/split_files/${val_set_id}.wavlist

# UBM 
ubm_file_dir=$PWD/ivec/ubm
ubm_file=$ubm_file_dir/models/ubm${ubm_size}.${feat_type}.train.param
echo $ubm_file

# create dirs
mkdir -p $odir/split_files $odir/split_files_stats $odir/ubm $odir/ivec_training $odir/out $odir/classification 

# change directory
cd ivec

if $initialize
  then

  # build file list for train and test set
  cp $dbfile_dir/${dbfile_train} $wav_list_train
  sed -i -e 's/^/'$(echo $wav_dir | sed 's/\//\\\//g')'/' $wav_list_train
  cp  $wav_list_train $odir/split_files/${train_set_id}.list
  sed -i -e 's/\/wavs/\/feats/' $odir/split_files/${train_set_id}.list
  sed -i -e 's/.wav/.'${feat_type}'.mfc/' $odir/split_files/${train_set_id}.list

  cp $dbfile_dir/${dbfile_test} $wav_list_val
  sed -i -e 's/^/'$(echo $wav_dir | sed 's/\//\\\//g')'/' $wav_list_val
  cp  $wav_list_val $odir/split_files/${val_set_id}.list
  sed -i -e 's/\/wavs/\/feats/' $odir/split_files/${val_set_id}.list
  sed -i -e 's/.wav/.'${feat_type}'.mfc/' $odir/split_files/${val_set_id}.list
  
  echo "file lists made"

  # split the training file list
  python Splitlist.py $odir/split_files/${train_set_id}.list $nb_splits $odir/split_files/${train_set_id}

  echo "file list train splitted in "$nb_splits
  
  # make label files
  python list2lst.py $odir/split_files/${train_set_id}.list $odir/split_files/${train_set_id}.lst
  for k in `seq 1 $nb_splits`; do
    python list2lst.py $odir/split_files/${train_set_id}_split_${k}.list $odir/split_files/${train_set_id}_split_${k}.lst
  done
  python list2lst.py $odir/split_files/${val_set_id}.list $odir/split_files/${val_set_id}.lst

fi

if $feat_extract
  then
  # extract features training set
  matlab_cmd="make_gfcc('"$wav_list_train"');"
  echo $matlab_cmd
  matlab -nojvm -r "$matlab_cmd" > $odir/out/feature_extraction_${train_set_id} 2>&1 &
  wait

  # extract features validation set
  matlab_cmd="make_gfcc('"$wav_list_val"');"
  echo $matlab_cmd
  matlab -nojvm -r "$matlab_cmd" > $odir/out/feature_extraction_${val_set_id} 2>&1 &
  wait

  echo "feature extraction finished"
fi

if $ubm_train
  then
  if [ ! -e $ubm_file ] 
  then
    echo "bash train_ubm.bash $feat_type $feat_dim $ubm_size $wav_list_train > $odir/ubm/ubm_training_${train_set_id}"
    bash train_ubm.bash $feat_type $feat_dim $ubm_size $wav_list_train > $odir/ubm/ubm_training_${train_set_id} &
    wait
  else
    echo "ubm file "$ubm_file" already exist"
  fi

  echo "UBM training finished"
fi

if $ubm_stats
  then
  # compute UBM stats for training set
  for k in `seq 1 $nb_splits`; do
  ./sc_compute_suf_stats_64 $ubm_file $odir/split_files/${train_set_id}_split_${k}.list $ubm_size $feat_dim $odir/split_files_stats/  > $odir/out/sc_compute_suf_stats_${train_set_id}_split_${k} & 
   wait
  done

  # compute UBM stats for validation set
  ./sc_compute_suf_stats_64 $ubm_file $odir/split_files/${val_set_id}.list $ubm_size $feat_dim $odir/split_files_stats/ > $odir/out/sc_compute_suf_stats_${val_set_id} & 
  wait

  echo "UBM stats calculation finished"
fi

if $prenorm_stats
  then
  # prenormalize UBM stats for training set
  for statsfile in `ls ${odir}/split_files_stats/${train_set_id}*.stats`; do
      matlab_cmd="stats2mat_feat('"$feat_type"','"$ubm_file"','"$statsfile"','"$odir"');"
      statsfile_id=`basename $statsfile`
      echo $matlab_cmd
      matlab -nojvm -r "$matlab_cmd" > $odir/out/stats2mat_feat_${statsfile_id} 2>&1 &
      wait
  done

  # prenormalize UBM stats for validation set
  for statsfile in `ls ${odir}/split_files_stats/${val_set_id}*.stats`; do
      matlab_cmd="stats2mat_feat('"$feat_type"','"$ubm_file"','"$statsfile"','"$odir"');"
      statsfile_id=`basename $statsfile`
      echo $matlab_cmd
      matlab -nojvm -r "$matlab_cmd" > $odir/out/stats2mat_feat_${statsfile_id} 2>&1 &
      wait
  done

  echo "UBM stats prenormalization finished"
fi

if $ivec_train
  then  
  # do ivec training
  matlab_cmd="run_ivec_training('"$feat_type"','"$train_set_id"',"$ivec_iter","$ivec_dim","$ubm_size",'"$odir"');"
  echo $matlab_cmd
  matlab -nojvm -r "$matlab_cmd" > $odir/out/run_ivec_training_${train_set_id}_iVitr${ivec_iter}_iVdim${ivec_dim} 2>&1 &
  wait

  echo "ivec training finished"
fi

if $ivec_extract
  then
  # extract ivec from training
  datasets="{'"$train_set_id"','"$val_set_id"'}"
  matlab_cmd="GetISDvector('"$train_set_id"',"$datasets",'"$feat_type"',"$ivec_dim","$ivec_iter",'"$odir"')"
  echo $matlab_cmd
  matlab -nojvm -r "$matlab_cmd" > $odir/out/GetISDvector_${train_set_id}_iVitr${ivec_iter}_iVdim${ivec_dim} 2>&1 &
  wait
  
  echo "ivec extraction finished"
fi

if $ivec_class
  then
  matlab_cmd="classification_svm('"$train_set_id"','"$val_set_id"','"$feat_type"',"$ivec_dim",'"$odir"'); "
  echo $matlab_cmd
  matlab -nojvm -r "$matlab_cmd" > $classdir/classification_svm_task_${val_set_id}_iVitr${ivec_iter}_iVdim${ivec_dim} 2>&1 &
  wait

  echo "ivec classification finished"
fi
