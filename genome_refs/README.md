# Translated Protein Coding Genome Reference

The basic Homo sapiens `genome.gtf` and `genome.fasta` that are used in nf-core currently (Dec 2023) do contain a protein coding gene. However, it is not sufficiently translated for ORF calling applications to return called ORFs. 

As a result there is a need for a good test example of a protein coding genome that will have a translated protein coding gene with alternative translation also. 

Below I document some notes as I progress through the process. In theory, replicating what is documented below will return the same reference files
## Notes

### Setup Directory
Create the directory structure. If you havent yet, `cd genome_refs`
```bash
mkdir data
mkdir notebooks
mkdir conda_yamls
``` 

### Candidate Region Selection
In order to find a good example region I went to [GWIPS-viz](https://gwips.ucc.ie) and copied the genomic coordinates to from the browser. 

Chose: `chr12:6,487,526-6,616,396` - [See GWIPS-viz](https://gwips.ucc.ie/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A6487526%2D6616396&hgsid=1426937_1a8pBmazwLZICEC52Lqso5QCht9w)

Set chr as an env variable: (note: I did this afterwards so this has not been tested. I hardcoded "chr12" into below commands)
```bash
CHR="chr12"
```
### File Download 
In order to generate a custom mini `genome` I must first download the full genome plus annotation. 

Annotation:
``` bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz -P data/
gzip -d data/gencode.v44.annotation.gtf.gz 
grep $CHR data/gencode.v44.annotation.gtf > data/${CHR}_gencode_v44.annotation.gtf
```

Genome: 
``` bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz -P data/ 
```

### Setup Env
Create the conda yaml for use in the notebook(s)

See [conda_yamls](conda_yamls/genome_ref.yml)

```bash
conda create -f conda_yamls/genome_ref.yml
```

### File Creation
Using a Jupyter notebook for scripting parse the FASTA and GTFs and update their sequence/coordinates as if the specified genome region is the entire genome.

See [notebook](notebooks/genome_ref.ipynb)


### Testing
To validate whether the new genome fasta will be useable I ran a STAR alignment against it with a HEK293 Ribo-Seq Bam I had locally. 

Build Index:
```bash 
conda create -f conda_yamls/star.yml
conda activate sample-dataset-creation_star

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir data/STAR/genome_index --genomeFastaFiles data/chr12.genome.fasta --sjdbGTFfile data/chr12.gtf
```

Get reads from [ENA](https://www.ebi.ac.uk/ena/browser/view/SRR1630831) they are already trimmed. 

Remove rRNA (Thanks to RiboCode Devs rRNA is [here](https://github.com/xryanglab/RiboCode/blob/master/data/rRNA.fa)):
```bash
conda create -f conda_yamls/bowtie.yml
conda activate sample-dataset-creation_bowtie

bowtie-build data/rRNA.fa data/rRNA
bowtie -p 8 data/rRNA -norc --un data/SRR1630831.lessrRNA.fq -q data/SRR1630831.fastq.gz data/HEK_RiboSeq_rRNA_reads.sam
```

Align to Genome:
```bash
conda activate sample-dataset-creation_star

STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir data/STAR/genome_index --readFilesIn data/SRR1630831.lessrRNA.fq  --outFileNamePrefix data/HEK_RiboSeq_nf_core_genome --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
```

Check alignments:
```bash
less data/HEK_RiboSeq_nf_core_genomeLog.final.out
```

Check RiboCode Periodicity:
```bash
conda create -f conda_yamls/ribocode.yml
conda activate sample-dataset-creation_ribocode

prepare_transcripts -g data/chr12.gtf -f data/chr12.genome.fasta -o data/ribocode_annotation
metaplots -a data/ribocode_annotation -r data/HEK_RiboSeq_nf_core_genomeAligned.toTranscriptome.out.bam

```

This failed for `chr12:6,487,526-6,616,396`. Need to scale up to increase depth on metagene

Rerun with full chr22 as it is smallest
```bash
CHR="chr22"
grep $CHR data/gencode.v44.annotation.gtf > data/${CHR}_gencode_v44.annotation.gtf

conda activate sample-dataset-creation_star

STAR --runThreadN 8 --genomeSAindexNbases 11 --runMode genomeGenerate --genomeDir data/STAR/chr22_genome_index --genomeFastaFiles data/full_chr22.genome.fasta --sjdbGTFfile data/chr22.gtf

STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir data/STAR/chr22_genome_index --readFilesIn data/SRR1630831.lessrRNA.fq  --outFileNamePrefix data/chr22_HEK_RiboSeq_nf_core_genome --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
conda deactivate

conda activate sample-dataset-creation_ribocode
prepare_transcripts -g data/chr22.gtf -f data/full_chr22.genome.fasta -o data/chr22_ribocode_annotation
metaplots -a data/chr22_ribocode_annotation -r data/chr22_HEK_RiboSeq_nf_core_genomeAligned.toTranscriptome.out.bam
conda deactivate

```

Again failed... 

Rerun with full genome 
```bash
gzip -d data/GRCh38.primary_assembly.genome.fa.gz --keep

conda activate sample-dataset-creation_star

STAR --runThreadN 8 --genomeSAindexNbases 11 --runMode genomeGenerate --genomeDir data/STAR/full_genome_index --genomeFastaFiles data/GRCh38.primary_assembly.genome.fa --sjdbGTFfile data/gencode.v44.annotation.gtf

STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir data/STAR/full_genome_index --readFilesIn data/SRR1630831.lessrRNA.fq  --outFileNamePrefix data/full_HEK_RiboSeq_nf_core_genome --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
conda deactivate

conda activate sample-dataset-creation_ribocode
prepare_transcripts -g data/gencode.v44.annotation.gtf -f data/GRCh38.primary_assembly.genome.fa -o data/full_ribocode_annotation
metaplots -a data/full_ribocode_annotation -r data/full_HEK_RiboSeq_nf_core_genomeAligned.toTranscriptome.out.bam
conda deactivate

```