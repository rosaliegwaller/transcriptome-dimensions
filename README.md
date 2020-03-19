# transcriptome-dimensions
Code to clean, normalize, and characterize transcriptome data in multi-dimensional space

### Download
* GTF:  
* Transcript based expression counts: MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz from https://research.themmrf.org/ (login required)

Files in data/transcriptome-dimensions

### Select first sample for each patient

### Select protein coding transcripts
```
cd data/transcriptome-dimensions

zgrep 'protein_coding' Homo_sapiens.GRCh37.74.gtf.gz > Homo_sapiens.GRCh37.74.protein_coding.gtf

grep -E 'transcript_id '$'"ENST'$'[0-9]{11}' -o Homo_sapiens.GRCh37.74.protein_coding.gtf \
| sed -e 's/transcript_id "//g' \
| uniq \
> Homo_sapiens.GRCh37.74.protein_coding.transcript_id.txt

gunzip -k MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz

grep -F --file=Homo_sapiens.GRCh37.74.protein_coding.transcript_id.txt \
MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt \
> MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts_Protein-Coding.txt
```
### Select transcripts with <=5% samples <=100 reads
```
cd data/transcriptome-dimensions

cat Transcripts_LT5_LT100.txt | sed -e 's/"//g' > transcripts_filtered.txt

grep -F --file=Transcripts_LT5_LT100.txt Homo_sapiens.GRCh37.74.protein_coding.gtf > Homo_sapiens.GRCh37.74.protein_coding.filtered.gtf


```
