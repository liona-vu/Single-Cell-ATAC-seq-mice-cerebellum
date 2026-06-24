# Single-cell-ATAQ-seq analysis of chromatin accessibility in the developing mice cerebellum

The project investigates the chromatin accessibility of the developing mice cerebellum using single-cell-ATAQ-seq analysis. The purpose of this project is to identify and compare cis-regulatory elements (CREs) such as promoters and enhancers in several cerebellar cell types. The project adapts the standard ArchR workflow (Granja et al., 2021) which can be accessed with the link: https://www.archrproject.com/bookdown/index.html

The dataset used for the project is based on the published paper:
Integrated single-cell transcriptomic and epigenetic study of cell state transition and lineage commitment in embryonic mouse cerebellum (Khouri-Farah et al., 2022)

Raw scATAC-seq fragments files from E12 - E14 and the metadata can be downloaded using the following links:
<br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178546 <br>
https://zenodo.org/records/6097549#.Ygw5My-B30o

Analysis was performed using a combination of a precompiled environment using a Docker/Apptainer image which can be accessed with this link:
https://hub.docker.com/r/greenleaflab/archr

a host server that runs RStudio, and the Nibi High Computing Cluster provdided by the Digital Research Alliance of Canada (DRAC).

Before beginning analysis, sign into the HPC cluster with your credentials, and convert docker image into Apptainer SIF image:

```
module load apptainer

apptainer build archr_image.sif docker://immanuelazn/archr:latest
```
<br>

Next, download the following software which will be critical to use for peak calling later:

```
module load scipy-stack/2026a

pip install MACS2
```
<br>

Download the fragment files that are required for this project.

```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393630/suppl/GSM5393630%5FE12%5FsnATAC%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393631/suppl/GSM5393631%5FE13%5FsnATAC%5Ffragments%2Etsv%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5393nnn/GSM5393632/suppl/GSM5393632%5FE14%5FsnATAC%5Ffragments%2Etsv%2Egz
```
<br>
Next, convert gzipped files into bgzipped files for further ArchR analysis. Recommended to send as a batch script since it can take some time. Make sure you are in the same directory as to the fragment file locations. (Can go grab a coffee while you wait)

```
module load samtools/1.22.1
module load htslib/1.22.1

gunzip *.tsv.gz

#Create file for fragment files, move fragment files, then cd into directory
mkdir fragment_files
mv *.tsv fragment_files/
cd fragment_files/

#bgzips the file, the -f flag overwrites the older file
bgzip -f *.tsv

#check for whether it actually got bgzipped by checking the file type
htsfile *.tsv.gz
```

Once you are done with the aforementioned steps, move the files to your local computer because the initial analysis will be performed on the RStudio host server.

To run a host server running RStudio, ensure you have the docker desktop app running the archr docker image. 
```
docker image pull immanuelazn/archr:latest
```
<br>

Then, run the following in terminal:

```
docker run -it --rm -v  your_directory_with_the_fragment_files:/home/rstudio -p 8787:8787 immanuelazn/archr:latest
```

<br>

Type http://localhost:8787/ in the search bar of the browser of your choice. In terminal, a password in red should show up. Type "rstudio" as the username and copy and paste the password in the password bar.

Run the R Script on the host R server.

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
