#!/bin/bash

HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality
while getopts 's:' option
do
    case "$option" in
        s) species=${OPTARG};;
    esac
done
mkdir ${HOME}/gene_position/${species}/combined
for gene in IGH IGK IGL
do
    cat ${HOME}/gene_position/${species}/pri/${species}_${gene}_pos.sorted.bed ${HOME}/gene_position/${species}/alt/${species}_${gene}_pos.sorted.bed > ${HOME}/gene_position/${species}/combined/${species}_${gene}_pos.sorted.bed
done