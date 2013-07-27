## Accurate estimates of allele-specific expression with a single genome sequence and RNA-seq data

### Quinn, A.M.; Juneja, P.; Jiggins, F.M.

For recreating the analysis described in the paper

---

These scripts are not meant to stand alone, but rather to serve as a record of our work.  They could also be used as a blueprint for future improvements to our pipeline.  Our intent is not to provide a polished piece of software, but to provide an example of how our methods might be implemented.

All the scripts have been tested, however, there may be path or software dependencies that need to be set up before they can be used.  I have done my best to indicate where this is needed. If you find that you are unable to fix a broken script, feel free to [contact us](http://www.gen.cam.ac.uk/research/jiggins/).

---

Running these scripts requires the Quinn2013_supporting_files directory available for download [here](https://www.gen.cam.ac.uk/research/jiggins/data/).

You will need to set the path to this directory in mglobals.py. You will also need to set the path to a new working directory in mglobals.py each time you run the pipeline.

Scripts are organized into folders which correspond to headers in the paper. These folders are described below:


## Alignment to a single reference

### single_alignment

Aligned RNA-seq data from F1 progeny of cross between two DGRP lines to published reference sequence.

Run:

    python pipeline-singleend_altfastas.py

Output is: "16_A12_pUn_down_INTER_py.csv" in working directory

### tophat_vary

Aligned to published reference with different tophat parameters

* 2 mismatches
* 3 mismatches
* 5 mismatches
* 10 mismatches

Each set of scripts is in the corresponding folder (except 2, this is the default, and is the same as that produced by "single_alignment."

Run:

    python <number>mismatches/pipeline-singleend_altfastas.py

Output is: "16_A12_pUn_down_INTER_py.csv"

### parental_alignment

Aligned samples to reference sequences created by substituting SNPs present in parental genomic sequences.

To generate reference sequences, set path to supporting files and run:

    build_fastas.py

Output is prefixed with "freeze2_362_noINDEL_1-1" and "freeze2_765_noINDEL_1-1."

You must manually add the output of build_fastas.py to what will be the working directory in mglobals.py.

Run:

    mkdir parental_alignment
    cd parental_alignment
    mkdir -p original/16_A12_pUn_down/
    mkdir -p alternate/16_A12_pUn_down/
    # move reference fastas and indices in 16_A12_pUn_down directories
    mv /path/to/<build_fastas.py output>/freeze2_362_noINDEL_1-1* original/16_A12_pUn_down/
    mv /path/to/<build_fastas.py output>/freeze2_765_noINDEL_1-1* alternate/16_A12_pUn_down/

Run:

    python pipeline-singleend_altfastas.py

Output is: "parental_alignment.csv"

## Alignment to multiple references

### multiple_reference/unphased

Alignment to the published reference sequence as well as one generated from SNPs called from the RNA-seq data. SNPs are unphased.

Run:

    python pipeline-singleend_altfastas.py

Output is: "16_A12_pUn_down_snps.csv"

### multiple_reference/phased

Input to hapCUT:

1. Set of filtered variants from multiple_reference_2genomes (file: "16_A12_pUn_down_freeze2.vcf")
2. Sorted bam file from multiple_reference_2genomes (file: "accepted_hits.sorted.bam")
3. Reference fasta from Quinn2013_supporting_files ("genome.fa")

To sort run:

    samtools sort accepted_hits.bam accepted_hits.sorted

For hapCUT run:

    /path/to/Quinn2013_supporting_files/HAPCUT-v0.5/x86_64/extractHAIRS --VCF 16_A12_pUn_down_freeze2.vcf --bam accepted_hits.sorted.bam --maxIS 600 --ref genome.fa --indels 1 > fragment_matrix_file
    /path/to/Quinn2013_supporting_files/HAPCUT-v0.5/x86_64/HAPCUT --fragments fragment_matrix_file --VCF 16_A12_pUn_down_freeze2.vcf --output output_haplotype_file > hapcut.log

To make phased vcfs run (there are some path dependencies to files produced in the above steps at the bottom of this file):

    python3 /path/to/phasing_SNPs/make_two_vcfs_from_phased.py

To build alternate fastas run (there are some path dependencies to files produced by make_to_vcfs_from_phased.py):

    python /path/to/phasing_SNPs/build_fastas.py

Finally, run (move to new fastas into place, run pipeline):

    # There are path dependencies you must fill in here
    bash /path/to/phasing_SNPs/prepare_env.sh
    python pipeline-singleend_altfastas.py

Output is: "phased_alignment.csv"

## SNP calling: using the RNA-seq data

### SNP_calling/no_filter

Same as multiple_reference_2genomes but with NO filtering by known SNPs.

Run:

    python pipeline-singleend_altfastas.py

Output is: "16_A12_pUn_down_snps.csv"

### SNP_calling/DGRP_all

Same as multiple_reference_2genomes but with filtering by ALL DGRP SNPs.

Run:

    python pipeline-singleend_altfastas.py

Output is: "16_A12_pUn_down_snps.csv"

### SNP_calling/DGRP_2homo

Same as multiple_reference_2genomes but with filtering by SNPs that were called homozygous in a minimum of 2 DGRP lines

Run:

    python pipeline-singleend_altfastas.py

Output is: "16_A12_pUn_down_snps.csv"

## Per-gene ASE estimates

### gene_estimates

This is the pipeline for getting per-gene ASE estimates using a single reference sequence and SNPs called from paired-end RNA-seq data.

This pipeline also incorporates the step where conflicting snps and reads are removed (see snp2gene.py)

Run:

    python pipeline-pairedend_altfastas.py

Output files end with "_genes.csv."

## Software Requirements (must be in PATH)

* TopHat (version 2.0.8 with Bowtie 2 version 2.1.0)
* SAMtools (version 0.1.18) or greater
* Bedtools intersect (version 2.17.0)
* Other software is in Quinn2013_supporting_files



