beastlier (BEAST Link Inference for Epidemiological Reconstruction)
========

This repository contains software for use with the method for transmission tree reconstruction in [BEAST](http://beast.bio.ed.ac.uk/) described in the paper:

Hall M, Woolhouse M, Rambaut A (2015) Epidemic Reconstruction in a Phylogenetics Framework: Transmission Trees as Partitions of the Node Set. PLoS Comput Biol 11(12): e1004613. doi:[10.1371/journal.pcbi.1004613](http://dx.doi.org/10.1371/journal.pcbi.1004613)

A full [BEAST 2](http://beast2.org/) package is planned. For the present, this repository contains two Python scripts:

* A script (ModifyXML.py) that takes BEAST 1 XML output from BEAUTi and transforms it for transmission tree reconstruction. This will use [lxml](http://lxml.de/) if it is installed, which is recommended as it produces nicer output, but will work without. Command line help is implemented (run the script with the -h argument). The tree prior specified in the original BEAUTi file will be overwritten with the BEASTLIER transmission tree version, so is irrelevent. Other models specified in BEAST (e.g. the nucleotide substitution model and molecular clock) are kept.
* A second script to calculate the maximum parent credibility (MPC) transmission tree from the .net file that BEAST produces when run using this method.

At present, unfortunately, ModifyXML.py works only with the development version of BEAST, which can be downloaded from [https://github.com/beast-dev/beast-mcmc](https://github.com/beast-dev/beast-mcmc). It is not compatible with version 1.8.3. Future BEAST 1 releases will be compatible.

A full tutorial and more detailed documentation will appear soon.
