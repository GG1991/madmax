#!/bin/bash

./run_1.sh 1 > out_1.dat
mv mic_info.dat mic_info_1.dat
./run_1.sh 2 > out_2.dat
mv mic_info.dat mic_info_2.dat
vim mic_info_1.dat mic_info_2.dat
