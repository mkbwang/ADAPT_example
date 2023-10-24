#!/bin/bash

declare -a setting1=("N100_P5" "N100_P10" "N100_P20" "N50_P5" "N50_P10" "N50_P20")
declare -a setting2=("enrich_abundant" "deplete_abundant" "mix_abundant" 
                    "enrich_rare" "deplete_rare" "mix_rare")
                    
for s1 in ${setting1[@]}; do
  mkdir $s1
  for s2 in ${setting2[@]}; do
    nested_dir="${s1}/${s2}"
    mkdir $nested_dir
  done
done
