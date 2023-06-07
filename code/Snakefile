import os
import re

# for ABySS
READS = {}  # sample names and the corresponding read files location
SAMPLES = []  # sample names = READS.keys()

# for BUSCO
AUGUSTUS_CONFIG_FILE = ""
BUSCO_CONFIG_FILE = ""
lineage = ""

# read through the configuration file
with open("snakemake.config", "r") as config:
    configs = list(config)
for line in configs:
    line = line.rstrip()
    if (re.match("^#", line)):
        continue
    elif (re.match("^sample: ", line)):
        name = re.search(": (\w+)$", line).group(1)
        if name not in SAMPLES:
            SAMPLES.append(name)
    elif (re.match("^read: ", line)):
        species = re.search("/(\w+)_[12]\.fastq", line).group(1)
        file = re.search("^read: (.+)$", line).group(1)
        if species in READS: 
            READS[species].append(file)
        else:
            READS[species] = []
            READS[species].append(file)
    elif (re.match("^k: ", line)):
        start = re.search("^k: ([0-9]+) ([0-9]+) ([0-9]+)$", line).group(1)
        end = re.search("^k: ([0-9]+) ([0-9]+) ([0-9]+)$", line).group(2)
        step = re.search("^k: ([0-9]+) ([0-9]+) ([0-9]+)$", line).group(3)
        K_values = range(int(start), int(end), int(step))
    elif (re.match("^AUGUSTUS_CONFIG: ", line)):
        AUGUSTUS_CONFIG_FILE = re.search("^AUGUSTUS_CONFIG: (.*)$", line).group(1)
        # print("AUGUSTUS_CONFIG: ", AUGUSTUS_CONFIG_FILE)
        shell("export AUGUSTUS_CONFIG_PATH='{}'".format(AUGUSTUS_CONFIG_FILE))
    elif (re.match("^BUSCO_CONFIG: ", line)):
        BUSCO_CONFIG_FILE = re.search("^BUSCO_CONFIG: (.+)$", line).group(1)
        # print("BUSCO_CONFIG:", BUSCO_CONFIG_FILE)
    elif (re.match("^lineage: ", line)):
        lineage = re.search("^lineage: (.+)$", line).group(1)
        GENES = list(map(lambda x: x.split(".")[0], os.listdir(os.path.join(lineage, "hmms"))))

rule all:
    input:
        "final.tree"

rule tree_all:
    input:
        "in.tree" 

rule extract_all:
    input:
        expand("BUSCOMP/{sample}", sample=SAMPLES) 

rule assess_all:
    input:
        expand("BUSCO/{sample}/{k}-{sample}", sample=SAMPLES, k=K_values) 

rule assemble_all:
    input:
        expand("genomes/{sample}/{k}-{sample}.fa", sample=SAMPLES, k=K_values) 

rule align_all:
    input:
        expand("msa/{gene}.fasta", gene=GENES)

rule genome_assembly:
    output: 
        "genomes/{sample}/{k}-{sample}.fa"
    run:
        cur_dir = os.getcwd()
        for k in K_values:
            if not os.path.isdir('{}/abyss/{}'.format(cur_dir, k)):
                path = os.path.join(cur_dir, "abyss/{}".format(k)) 
                os.makedirs(path)

        read_files = READS
        reads = read_files[wildcards.sample]
        shell("module load abyss/2.2.0")
        if len(reads) == 2:
            species = wildcards.sample
            shell("abyss-pe k={} --directory './abyss/{}' in='{} {}' name={}".format(wildcards.k, wildcards.k, reads[0], reads[1], species))
            # Cleaning up:
            # Our target file species-scaffolds.fa is a soft link, so we convert it to be the original file
            shell("cp --remove-destination \"./abyss/{}/$(readlink ./abyss/{}/{}-scaffolds.fa)\" ./abyss/{}/{}-scaffolds.fa".format(wildcards.k, wildcards.k, species, wildcards.k, species))
            # Remove files that are "{species}-*" but not "{species}-scaffolds.fa"
            # shell("find ./abyss/ -type f ! -name '{}-scaffolds.fa' -delete".format(species))
            shell("find ./abyss/{} | grep '{}-*' | grep -v '{}-scaffolds.fa' | xargs rm -f".format(wildcards.k, species, species))
            # # Remove 'coverage.hist' file 
            shell("find ./abyss/{} -type f -name 'coverage.hist' -delete".format(wildcards.k))
        elif len(reads) == 1:
            species = wildcards.sample
            shell("ABYSS -k{} {} -o ./abyss/{}/{}-scaffolds.fa".format(wildcards.k, reads[0], wildcards.k, species))
        shell("mv ./abyss/{}/{}-scaffolds.fa ./genomes/{}/{}-{}.fa".format(wildcards.k, wildcards.sample, wildcards.sample, wildcards.k, wildcards.sample))

rule BUSCO_assessment:
    input:
        "genomes/{sample}/{k}-{sample}.fa"
    output:
        directory("BUSCO/{sample}/{k}-{sample}")
    shell:
        """
        module add python/3.8.2 sepp/4.3.10 blast+/2.10.1 hmmer/3.3 augustus/3.3.2 busco/4.1.0;
        export BUSCO_CONFIG_FILE='{BUSCO_CONFIG_FILE}';
        export AUGUSTUS_CONFIG_PATH="{AUGUSTUS_CONFIG_FILE}"
        busco -q -m genome -i {input} -o {wildcards.k}-{wildcards.sample} -l {lineage};
        mv {wildcards.k}-{wildcards.sample} {output};
        """

rule BUSCO_plot:
    input:
        expand("BUSCO/{sample}/{k}-{sample}", sample=SAMPLES, k=K_values)
    output:
        "busco_figure.png"
    shell:
        """
        mkdir busco_summary;
        for assessment in BUSCO/*; do
            for d in $assessment/*; do
                cp $d/short_summary* ./busco_summary
            done
            
        done
        touch busco_figure.png;
        module unload busco/4.1.0;
        module add python/3.8.2 sepp/4.3.10 blast+/2.10.1 hmmer/3.3 augustus/3.3.2 emboss/6.6.0 busco/3.0.2b;
        export BUSCO_CONFIG_FILE="./busco_plot.config.ini";
        python3 generate_plot.py -wd busco_summary;
        mv busco_summary/busco_figure.png .;
        rm -r busco_summary;
        module unload busco/3.0.2b;
        module load busco/4.1.0;
        """

rule BUSCOMP_extraction:
    input:
        busco = expand("BUSCO/{{sample}}/{k}-{{sample}}", k=K_values),
        genomes = expand("genomes/{{sample}}/{k}-{{sample}}.fa", k=K_values) 
    output:
        directory("BUSCOMP/{sample}")
    shell:
        """
        export RSTUDIO_PANDOC=/apps/rstudio/1.1.456/bin/pandoc;
        mkdir BUSCOMP/{wildcards.sample};
        module load minimap2/2.17;
        python /home/z3452659/slimsuitedev/tools/buscomp.py basefile={wildcards.sample} runs=./BUSCO/{wildcards.sample}/* fastadir=./genomes/{wildcards.sample} i=-1;
        mv {wildcards.sample}* BUSCOMP/{wildcards.sample};
        """

rule group_genes:
    input:
        expand("BUSCOMP/{sample}", sample=SAMPLES)
    output:
        directory("parsed_genes")
    shell:
        """
        python3 parse_genes.py;
        """

rule alignment: 
    input: 
        "parsed_genes",
    output:
        "msa/{gene}.fasta"
    shell: 
        '''
        module load muscle/3.8.31;
        if ! test -d "./msa"; then 
         mkdir msa; 
        fi 
        infile=./parsed_genes/{wildcards.gene}.fasta;
        out=msa/{wildcards.gene}.fasta;
        if ! test -f $infile; then
            touch $out;
        else 
            muscle -quiet -in $infile -out $out -maxiters 2;
        fi
        '''

rule build_tree:
    input:
        "msa/{gene}.fasta"
    output:
        "trees/{gene}.treefile"
    shell:
        """
        module load iqtree/2.0.4;
        if ! test -d "./trees"; then 
            mkdir trees; 
        fi
        if [ ! -s ./msa/{wildcards.gene}.fasta ]; then
            touch ./trees/{wildcards.gene}.treefile;
        else
            iqtree2 -s ./msa/{wildcards.gene}.fasta -pre {wildcards.gene} -quiet;
            mv {wildcards.gene}.treefile ./trees;
            ls | egrep "{wildcards.gene}" | xargs rm -f;
        fi
        
        """

rule collect_trees:
    input:
        expand("trees/{gene}.treefile", gene=GENES)
    output:
        "in.tree"
    shell:
        """
        if ! test -f "./in.tree"; then 
            touch in.tree; 
        fi
        for file in ./trees/*; do
            cat $file >> in.tree;
        done
        """

rule super_tree:
    input:
        "in.tree"
    output:
        "final.tree"
    shell:
        """
        java -jar ./ASTRAL/Astral/astral.5.7.4.jar -i in.tree -o final.tree
        """
