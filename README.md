# GRASS
========

[![Join the chat at https://gitter.im/COMBINE-lab/GRASS](https://badges.gitter.im/COMBINE-lab/GRASS.svg)](https://gitter.im/COMBINE-lab/GRASS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

GRASS (Graph Regularized Annotation via Semi-Supervised learning) is a tool for annotating *de novo* transcriptome assemblies using data from closely related species with previously annotated genomes. GRASS efficiently utilizes information from three major sources:
  1. A graph generated by [RapClust](https://github.com/COMBINE-lab/RapClust), which is a tool for clustering contigs in a *de novo* assembly. The graph is constructed using previously computed *fragment equivalence classes*. 
  2. Quantification results obtained using [Sailfish](https://github.com/kingsfordgroup/sailfish) or [Salmon](https://github.com/COMBINE-lab/salmon) which are used to improve the graph and filter out spurious contigs from the assembly. (*fragment equivalence classes*) are also generated by these.
  3. Annotation information used as *seed* to label the graph. This is obtained by a two-way nucleotide BLAST by default but any other labeling method can be used as well.

## Dependencies
----------------

GRASS depends on the following external programs (to be available in the environment where it runs):

  1. The [MCL](http://micans.org/mcl/) clustering tool
  2. The [Sailfish](https://github.com/kingsfordgroup/sailfish) (or [Salmon](https://github.com/COMBINE-lab/salmon)) quantification tool.
  3. The [junto](https://github.com/parthatalukdar/junto) library for label propagation.

Further, it depends on the following Python packages:
  
  1. [RapClust](https://pypi.python.org/pypi/rapclust/)
  2. [Click](http://click.pocoo.org/5/)
  3. [PyYAML](https://pypi.python.org/pypi/PyYAML)
  4. [Pandas](http://pandas.pydata.org/)
  5. [NumPy](http://www.numpy.org/)

If you wish to use the default labeling method in GRASS and allow it to run BLAST on the input FASTA files, ensure that you have the BLAST utilities installed (ncbi-blast+). 

Also ensure that junto and MCL are installed. For junto, download the library source code from the github page and set JUNTO_DIR to be the top level directory of junto. Next, add the directory JUNTO_DIR/bin to your path "export PATH="$PATH:$JUNTO_DIR/bin". A ready-to-install tarball for MCL can be downloaded from their webpage.

To install GRASS via pip (along with the dependencies getting installed automatically), you can use:

```
> pip install grass
```

You should now have a `grass` executable in your path.  You can test this with the following command:

```
> grass --help
```

You should see the following output:

```
Usage: grass [OPTIONS]

Options:
  --config TEXT  Config file describing the experimental setup
  --help         Show this message and exit.
```

## Using GRASS
---------------

The code for GRASS automatically runs RapClust to generate the mapping ambiguity graph. However, Sailfish or Salmon need to be run on the input samples to provide the information required in order to generate this graph. Hence, there are two main steps involved:

  1. Run Sailfish on each sample in your experiment, passing it the `--dumpEq` option.  This will tell Sailfish to dump a representation of the fragment equivalence classes that it computed during quasi-mapping of each sample.  Apart from this additional option, Sailfish should be run normally (i.e. passing in whatever other options are appropriate for your samples).
  2. Run GRASS, providing it with a configuration file that describes the experimental setup of your samples, and where the Sailfish quantification results have been written.

### Run Sailfish/Salmon
Let's illustrate this pipeline with a particular example, the following experimental data from *[Trapnell et al.](http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html)*:

Accession | Condition | Replicate
----------|-----------|----------
SRR493366 | scramble  | 1
SRR493367	| scramble  | 2
SRR493368	| scramble  | 3
SRR493369	| HOXA1KD	  | 1
SRR493370	| HOXA1KD	  | 2
SRR493371 | HOXA1KD   | 3

We'll assume that the raw read files reside in the directory `reads`.  Assuming that you've already built the index on the transcriptome you wish to quantify, a typical run of Sailfish on this data would look something like.

```
> parallel -j 6 "samp={}; sailfish quant -i index -l IU -1 <(gunzip -c reads/{$samp}_1.fq.gz) -2 <(gunzip -c reads/{$samp}_2.fq.gz) -o {$samp}_quant --dumpEq -p 4" ::: SRR493366 SRR493367 SRR493368 SRR4933669 SRR493370 SRR493371
```

This will quantify each sample, and write the result to the directory `samplename_quant`. 

### Setting up config file for GRASS
Our test species in this example is human and the closely related annotated species is chimp. Given this setup, we're now ready to run GRASS.  First, we have to make an appropriate config file, an example config file looks like following:
```
conditions:
    - Control
    - HOXA1 Knockdown
samples:
    Control:
        - SRR493366_quant
        - SRR493367_quant
        - SRR493368_quant
    HOXA1 Knockdown:
        - SRR493369_quant
        - SRR493370_quant
        - SRR493371_quant
outdir: human_rapclust
```

Along with the above information GRASS requires related species information in *any* **one** of the following formats:

  1. You can pass the FASTA files to GRASS in the following way and it will run a two-way BLAST assigning seed labels to contigs. Ensure that the FASTA files are passed in the following order (the first is from the test species, second from the annotated species)

  ```
  fasta:
      - human.transcripts.fa
      - chimp.transcripts.fa
  ```

  2. If you have already run BLAST, you can pass the output files (in BLAST outfmt 6). Again, ensure that the first one is BLAST of contigs from test species against the annotated species and the second is BLAST of contigs from annotated species against the test species. 

  ```
  labels:
      - human.chimpdb.txt
      - chimp.humandb.txt
  ```

  3. If you wish to use a pre-processed label file, you can pass a two-column file where the first is the set of contigs from the test species and second the label. If a contig has multiple labels in the input file, one will be chosen arbitrarily as seed.

  ```
  labels:
      - human.labels.txt
  ```

So a sample config file provided with the FASTA files (example 1) would look something like this:
```
conditions:
    - Control
    - HOXA1 Knockdown
samples:
    Control:
        - SRR493366_quant
        - SRR493367_quant
        - SRR493368_quant
    HOXA1 Knockdown:
        - SRR493369_quant
        - SRR493370_quant
        - SRR493371_quant
fasta:
    - human.transcripts.fa
    - chimp.transcripts.fa
outdir: human_rapclust
```

### Run GRASS
Once we have our `config.yaml` file ready with the above information we can run GRASS.  GRASS uses [YAML](http://yaml.org/) to specify its configuration files.  The configuration file must contain the following three entries; `conditions`, `samples`, and `outdir` and only one of the following two: `fasta` or `labels`.  The `conditions` entry lists the conditions present in the sample. The `samples` entry is a nested dictionary of lists; there is a key corrseponding to each condition listed in the `conditions` entry, and the value associated with this key is a list of quantification directories of the samples for this condition.  Finally, the `outdir` entry specifies where the output and intermediate files should be stored. The last entry has been explained above. Given these, you can run GRASS as:

```
> grass --config config.yaml
```

This will process the samples, generate the mapping ambiguity graph, filter it according to the conditions, process the contig labels, run the iterative GRASS algorithm to improve the graph and cluster the resuling graph using [MCL](http://micans.org/mcl/) internally.  Once GRASS is finished, the `human_rapclust` directory should exist.  It will contain the following files which are results from GRASS plus the graph file from RapClust:

`seedLabels.txt, grassGraph.txt, grass.mag.clust, grass.mag.flat.clust, grass.stats.json, mag.filt.net, finalLabels.txt, stats.json`

The most important file for downstream processing is `grass.mag.flat.clust`.  It contains the computed cluster information in a "transcript-to-gene" mapping formation that is compatible with downstream tools like [tximport](https://github.com/mikelove/tximport).  The file `seedLabels.txt` contains the initial contig to gene labeling and `finalLabels.txt` contains the labels after running GRASS and a contig may have multiple labels in this file, each with an associated score. The rest of the files are for internal use in the algorithm.


## Citations:
-------------

Experiments in Graph-based Semi-Supervised Learning Methods for Class-Instance Acquisition. Partha Pratim Talukdar, Fernando Pereira, ACL 2010

Differential analysis of gene regulation at transcript resolution with RNA-seq by Cole Trapnell, David G Henderickson, Martin Savageau, Loyal Goff, John L Rinn and Lior Pachter, Nature Biotechnology 31, 46–53 (2013).

Stijn van Dongen. Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, 2000

Charlotte Soneson, Michael I Love, and Mark D Robinson. Differential analyses for rna-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4, 2015.


