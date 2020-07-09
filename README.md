# ECCOreproduce

We provide the detailed analysis scripts for ECCO here. 

Those analysis scripts require the R packages: data.table, peer, MatrixEQTL, doParallel and TwoSampleMR.

The R package TwoSampleMR could be installed through:

```
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
```
And here is the detail information about the installation of R package peer:
https://github.com/PMBio/peer/wiki/Installation-instructions

We also provide the detailed analysis scripts for real data in the real data folder.

Data file is the pipeline about how we process the GTEx original data

Detailed analysis procedure: https://github.com/fanyue322/ECCOreproduce/blob/master/analysis%20notebook.pdf

The GTEx eQTL　mapping results are hosting at: <https://github.com/fanyue322/ECCOreproduce/tree/master/eQTL_mapping_results>

The format of these results are:
```
chromosome   tissue   gene_name(Ensembl Gene ID)  number_of_peer_factors_included   beta      p_value      eGene_or_not(0,1)
    1         Liver         ENSG00000013573                 50                     -1.12     1.949e-163          1
    2         Liver         ENSG00000203875                 50                      0.88     1.215e-154          1
    ..         ..                 ..                        ..                       ..          ..             ..
```
