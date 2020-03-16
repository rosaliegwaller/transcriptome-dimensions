# transcriptome-dimensions
Code to clean, normalize, and characterize transcriptome data in multi-dimensional space

### Download
* GTF: ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz
* Transcript based expression counts: MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz from https://research.themmrf.org/ (login required)

Files in data/transcriptome-dimensions

### Code
```
cd data/transcriptome-dimensions

zgrep 'protein_coding' Homo_sapiens.GRCh37.74.gtf.gz > Homo_sapiens.GRCh37.74.protein_coding.gtf
```
