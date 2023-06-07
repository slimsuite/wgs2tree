#!/bin/sh

for dir in "abyss" "genomes" "BUSCO" "BUSCOMP" "msa" "parsed_genes" "trees"; do
    if [ -d ./$dir ]; then
        rm -r ./$dir;
    fi
done

rm in.tree
