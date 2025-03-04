#!/bin/bash

for k in {2..10}; do
  python assemble_genome.py -i input/reads.fna -k $k -t 1 --debug --seed 42 --clean_graph --output_dir output_dir/ > logs/k_${k}.txt
done
