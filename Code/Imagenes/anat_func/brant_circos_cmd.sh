#!/usr/bin/env bash

chmod u+x "/Users/hose/Downloads/circos-0.69-5/circos"
circosbin="/Users/hose/Downloads/circos-0.69-5/circos"
circosconf="/Applications/MATLAB_R2015b.app/toolbox/brant-stable-master/circos/circos.conf"
labelfile="/Users/hose/Desktop/TFM_TECI/MINT_trabajo/brant_labels.txt"
bandfile="/Users/hose/Desktop/TFM_TECI/MINT_trabajo/brant_band.txt"
linkfile="/Users/hose/Desktop/TFM_TECI/MINT_trabajo/brant_links.txt"
outdir="/Users/hose/Desktop/TFM_TECI/MINT_trabajo"
${circosbin} -conf ${circosconf} -param image/background=white -param karyotype=${bandfile} -param plots/plot/file=${labelfile} -param links/link/file=${linkfile} -outputdir ${outdir}