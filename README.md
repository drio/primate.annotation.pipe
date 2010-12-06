primate annotation pipeline
===========================

#### Introduction ####

These set of tools implement a pipeline to perform annotation of primate data 
The original input is SNP calls from our Sanger pipeline (SNP detector). The
typical data flow would look like this:

[![Diagram](http://drio.github.com/primate.annotation.pipe/images/annotation.primate.diagram.png)](http://github.com/drio/primate.annotation.pipe)

#### Requirements: ####

1. ruby interpreter
2. perl distribution with ensembl APi packages and its dependencies. 

#### Installation: ####

1. Download the tools. You can clone with git or copy over the directory 
to your destination machine.

2. Add an entry in your PATH to the bin directory.

At this point you should be able to use any of the tools in bin.

#### Input files ####

You should have 4 files as input files:

1. snps.txt          : Alleles at the snps position for ALL the samples
2. snps.summary.txt  : Genotype call at the snp positions
3. indels.txt        : Alleles at the indel positions.
4. indels.summary.txt: Genotype call at the indel positions.
