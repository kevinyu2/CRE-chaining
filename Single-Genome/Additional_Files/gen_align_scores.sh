#!/bin/bash

rand_dir="/home/mwarr/Data/alignment/rand_temp"
setb_dir="/home/mwarr/Data/alignment/setb_temp" 
align_dir="/home/mwarr/Data/alignment/local/align_temp"
script_dir="/home/mwarr/CRE-chaining/Single-Genome/alignment_rand_vs_acr.py"

for i in {0..49};
do 
    nohup python $script_dir ${rand_dir}/temp_${i} ${align_dir}/temp_rand_${i} &
    nohup python $script_dir ${setb_dir}/temp_${i} ${align_dir}/temp_${i} &
done;




