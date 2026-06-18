# Single-cell-ATAQ-seq analysis of chromatin accessibility in the developing mice cerebellum

The project investigates the chromatin accessibility of the developing mice cerebellum using single-cell-ATAQ-seq analysis. The purpose of this project is to identify cis-regulatory elements (CREs) such as promoters and enhancers in several cerebellar cell types. The project adapts the standard ArchR workflow (Granja et al., 2021) which can be accessed with the link: https://www.archrproject.com/bookdown/index.html

The dataset used for the project is based on the published paper:
Integrated single-cell transcriptomic and epigenetic study of cell state transition and lineage commitment in embryonic mouse cerebellum (Khouri-Farah et al., 2022)

Raw scATAC-seq fragments files from E12 - E14 and the metadata can be downloaded using the following links:
<br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178546 <br>
https://zenodo.org/records/6097549#.Ygw5My-B30o

Analysis was performed using a combination of a precompiled environment using a Docker/Apptainer image which can be accessed with this link:
https://hub.docker.com/r/greenleaflab/archr

and the Nibi High Computing Cluster provdided by the Digital Research Alliance of Canada (DRAC).

#### Software used:
- R Studio v4.4.1 (2024-06-14) (based on the Docker image) <br>
- ArchR v1.0.3 (2024-11-26) (based on the Docker image) <br>
- MACS2 (Feng et al., 2012)

#### HPC Modules used:
- apptainer v1.4.5 <br>
- samtools v1.22.1 <br>
- htslib v1.22.1

#### Citations
Khouri-Farah, N., Guo, Q., Morgan, K., Shin, J. & Li, J. Y. H. Integrated single-cell transcriptomic and epigenetic study of cell state transition and lineage commitment in embryonic mouse cerebellum. Science Advances 8, eabl9156.

Granja, J. M. et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics 53, 403–411 (2021).

Feng, J., Liu, T., Qin, B., Zhang, Y. & Liu, X. S. Identifying ChIP-seq enrichment using MACS. Nature Protocols 7, 1728–1740 (2012).
