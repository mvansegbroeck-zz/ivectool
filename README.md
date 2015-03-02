# ivectool
## 1. extract audio: 
tar -xvzf audio.tar.gz

## 2. convert flac to wav: 

mkdir wavs; 

for line in flac/*.flac; 

 do sox $line wavs/$(basename ${line%.*}).wav ;
 
done

## 3. run pipeline script

bash run_system_train-validation.bash
