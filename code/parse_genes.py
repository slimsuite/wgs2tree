#!/apps/python/3.8.3/bin/python3
import sys, os, re
from Bio import SeqIO

# working_dir = /srv/scratch/z5206348/

working_dir = os.getcwd() + "/"
busco_genes = []

os.mkdir(working_dir + "parsed_genes/")

# open readfiles config file and store contents as text
with open("snakemake.config", "r") as config:
    file_contents = config.read()

# store sample names
samples = re.findall(r"sample: (\w+)", file_contents)
print(samples)

for sample in samples:
    sample_dir = working_dir + "BUSCOMP/" + sample + "/" + sample + ".buscomp.fasta"
    with open(sample_dir, "r") as buscomp_output:
        for record in SeqIO.parse(buscomp_output, "fasta"):
            output_filename = working_dir + "parsed_genes/" + record.id + ".fasta"
            record.id = sample
            mode = 'a' if os.path.exists(output_filename) else 'w'
            with open(output_filename, mode) as output_file:
                SeqIO.write(record, output_file, "fasta")
