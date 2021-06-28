#! /usr/bin/env bash
input_files="${1}"
chewbacca_result="${2}"
missing_loci="${3}"

# validate inputs
[ ! -f ${input_files} ] && echo "no file ${input_files}" && exit 1
[ ! -f ${chewbacca_result} ] && echo "no file ${chewbacca_result}" && exit 1
[ -z ${missing_loci} ] && echo "no missing loci output file set (arg nr 3)" && exit 1

# iterate over results and calculate the number of missing loci
wc -l "${input_files}" | cut -d ' ' -f 1 | xargs -I '{}' tail -n '{}' "${chewbacca_result}" | while read line || [ -n "${line}" ]; do
    n_missing=`echo ${line} | sed -e "s/NIPHEM/-/g" -e "s/NIPH/-/g" -e "s/LNF/-/g" -e "s/INF-*//g" -e "s/PLOT[^\t]*/-/g" -e "s/ALM/-/g" -e "s/ASM/-/g" | fmt -w 1 | tail -n +2 | grep '-' | wc -l`
    echo -e "$(echo ${line} | cut -d ' ' -f 1)\t${n_missing}" >> "${missing_loci}"
done;
