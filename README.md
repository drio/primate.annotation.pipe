primate annotation pipeline
===========================

#### Introduction ####

These set of tools implement a pipeline to perform annotation of primate data 
(Check diagram in docs/).  The original input is SNP calls from our Sanger 
pipeline (SNP detector).

[![Diagram](http://drio.github.com/primate.annotation.pipe/images/annotation.primate.diagram.png)](http://github.com/drio/primate.annotation.pipe)

#### Requirements: ####

  1. ruby interpreter
  2. perl distribution with ensembl APi packages and its dependencies. 

#### Installation: ####

1. Download the tools. You can clone with git or copy over the directory 
to your destination machine.

2. Add an entry in your PATH to the bin directory.

At this point you should be able to use any of the tools:

<pre>
$ ls bin/*
-rwxr-xr-x  1 drio  staff   1.1K Oct 15 15:30 ensembl2lff.rb
-rwxr-xr-x  1 drio  staff   3.1K Oct 14 10:55 primate_bridge.rb
-rwxr-xr-x  1 drio  staff   1.1K Oct 14 11:17 run_snp_effect_predictor.sh
</pre>

