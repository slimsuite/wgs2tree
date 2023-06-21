# wgs2tree - rapid phylogenomics from whole genome sequencing data

_Stephanie H. Chen, Caitlin Ramsay, Jade Liang, Weilin Wu, Ying Xu, Kavitha Krishna & Richard J. Edwards_

`wgs2tree` is a snakemake workflow for generating consensus phylogenomic trees from whole genome sequencing (WGS) data. It was originally coded as a UNSW BINF6112 project under the supervision of [Stephanie Chen](https://github.com/stephanie-h-chen) and [Richard Edwards](https://github.com/cabbagesofdoom). The students were 
Caitlin Ramsay, Jade Liang, Weilin Wu, Ying Xu, and Kavitha Krishna. Limited documentation can be found in the provided [Technical Specification](https://github.com/slimsuite/wgs2tree/blob/main/wgs2tree%20Tech%20Spec.pdf).

The `wgs2tree` workflow is as follows:

1. Generate low coverage WGS data (e.g. Illumina) per sample.
2. Generate a number of quick draft assemblies per sample.
3. Use BUSCO to identify single-copy orthologues across assemblies.
4. Use BUSCOMP to compile the maximal set of BUSCO orthologues for each sample.
5. Generate multiple sequence alignments and phylogenetic trees per BUSCO gene across samples.
6. Generate a consensus phylogenomic tree with ASTRAL.

If you would like more information, or for help adapting to different data types, please raise an Issue or contact the authors.

## License

`wgs2tree` is freely available for use under a [GNU GPL v3 license](https://github.com/slimsuite/wgs2tree/blob/main/LICENSE).

## Citation

`wgs2tree` has not yet been published. Please cite this repository and the tools used by the `wgs2tree` workflow.



