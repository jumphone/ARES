<img src="https://github.com/jumphone/PhenoPro/raw/master/IMG/ARES_logo_text.png" width="360">

ARES is designed for detecting RNA editing sites from RNA-seq data.

### Requirements

    python3  = 3.8.2
    pysam    = 0.16.0.1
    numpy    = 1.19.1
    bedtools = v2.26.0    # https://bedtools.readthedocs.io/en/latest/
    blat     = v.36x7     # http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/



### Usage

### ARES-anno (use ARES with repeat annotation file)

Users can download repeat annotation file from: https://sourceforge.net/projects/sprintpy/files/dbRES/, "dbrep.zip"
    
    python3  ares.py  bam_in_path  ref_in_path  OUT_DIR  bedtools_path  blat_path  anno_in_path
    
    # bam_in_path: path to aligned reads in BAM format
    # ref_in_path: path to reference genome in FASTA format
    # OUT_DIR: ARES will generate a folder to store all results
    # bedtools_path: path to bedtools
    # blat_path: path to blat
    # anno_in_path: path to repeat annotation file. Uers can download it from https://sourceforge.net/projects/sprintpy/files/dbRES/, "dbrep.zip"
 
 
### ARES-free (use ARES without repeat annotation file)    
    
    
    python3  ares.py  bam_in_path  ref_in_path  OUT_DIR  bedtools_path  blat_path 
    
    # bam_in_path: path to aligned reads in BAM format
    # ref_in_path: path to reference genome in FASTA format
    # OUT_DIR: ARES will generate a folder to store all results
    # bedtools_path: path to bedtools
    # blat_path: path to blat

### Output folder

    # fc_res_dsrna.bed: RESs identified by dsRNA-based part
    # ff_res_anno.bed: RESs identified by annotation-based part
    # fg_res_all.bed: all RESs identified by ARES


### Recommended alignment procedure (BWA-MEM)

bwa: http://bio-bwa.sourceforge.net/

samtools: http://www.htslib.org/download/

bamUtil: https://github.com/statgen/bamUtil
    
    # Build index
    bwa  index  -a  bwtsw  reference.fasta 
    
    # Do reads aligment
    bwa  mem  -t  cpu_number  reference.fasta  read_1.fastq  read_2.fastq  >  reads.sam

    # Sort & remove PCR duplicates
    samtools  view  -b  reads.sam  >  reads.bam
    samtools  sort  reads.bam  >  reads.sorted.bam
    bam  dedup_LowMem  --rmDups  --in  reads.sorted.bam  --out  reads.sorted.rmdup.bam  --log  log.txt 
    samtools  index  reads.sorted.rmdup.bam 
    samtools  stats  reads.sorted.rmdup.bam  >  reads.stats.txt   
   
   
### Identify hyper-RESs (SPRINT)

sprint: https://github.com/jumphone/SPRINT

    samtools  view  -f4  -b  reads.sorted.rmdup.bam  >  reads.unaligned.bam
    samtools  bam2fq  reads.unaligned.bam  >  reads.unaligned.fastq
    
    sprint  main  -1  reads.unaligned.fastq  reference.fasta  OUT_DIR  bwa_path  samtools_path
    
    # Please check https://github.com/jumphone/SPRINT for the detailed usage of SPRINT
    # when using sprint, the version of bwa should be 0.7.12, and the version of samtools should be 1.2.
  
  
  
  
