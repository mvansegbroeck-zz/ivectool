# ivectool
extract audio: tar -xvzf audio.tar.gz
convert flac to wav: mkdir wavs; for line in flac/*.flac; do sox $line wavs/$(basename ${line%.*}).wav ; done
