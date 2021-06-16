#!/usr/env bash

definitions=(bedtools sambamba picard)
declare -A containers=( 
    [bwa]=https://depot.galaxyproject.org/singularity/bwa%3A0.7.17--pl5.22.0_2
    [kraken]=https://depot.galaxyproject.org/singularity/kraken:1.1.1--pl5262h7d875b9_5
    [kraken2]=https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5262h7d875b9_0
    [bracken]=https://depot.galaxyproject.org/singularity/bracken:2.6.1--py39h7cff6ad_2
    [samtools]=https://depot.galaxyproject.org/singularity/samtools%3A1.12--hd5e65b6_0
    [ariba]=https://depot.galaxyproject.org/singularity/ariba:2.14.6--py38h6ed170a_0
    [mlst]=https://depot.galaxyproject.org/singularity/mlst:2.19.0--hdfd78af_1
    [spades]=https://depot.galaxyproject.org/singularity/spades:3.15.2--h95f258a_1
    [quast]=https://depot.galaxyproject.org/singularity/quast%3A5.0.2--py37pl5262h190e900_4
    [freebayes]=https://depot.galaxyproject.org/singularity/freebayes%3A1.3.5--py39hba5d119_3
)

for tool in "${definitions[@]}"; do
    output_file="${tool}_$(date +%Y-%m-%d).sif";
    if [[ ! -f $output_file ]]; then
        echo "Building tool ${tool}";
        singularity build --force "${output_file}" "${tool}";
    else
        echo "Tool ${tool} already exist, skipping...";
    fi;
done;
echo "Download pre built singularity containers";
for tool in ${!containers[@]}; do
    output_file="${tool}.sif";
    if [[ ! -f $output_file ]]; then
        echo "Downloading ${tool}";
        wget -O "${tool}.sif" "${containers[$tool]}"
    else
        echo "Tool ${tool} already exist, skipping...";
    fi;
done
