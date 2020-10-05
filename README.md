<img src="https://github.com/jumphone/PhenoPro/raw/master/IMG/ARES_logo_text.png" width="360">

### Requirements

    python3  = 3.8.2
    pysam    = 0.16.0.1
    numpy    = 1.19.1
    bedtools = v2.26.0    # https://bedtools.readthedocs.io/en/latest/
    blat     = v.36x7     # http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/

Users can download binary files of bedtools and blat from: https://github.com/jumphone/ARES/tree/master/bin

Users can download repeat annotation file from: https://sourceforge.net/projects/sprintpy/files/dbRES/, "dbrep.zip"

### Usage

### ARES-anno (use ARES with repeat annotation file)
    
    
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


### Recommended alignmnet procedure (BWA-MEM)

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
    samtools index reads.sorted.rmdup.bam 
    samtools stats reads.sorted.rmdup.bam > reads.stats.txt   
   
   
### Identify hyper-RESs (SPRINT)

sprint: https://github.com/jumphone/SPRINT

    samtools  view  -f4  -b  bam_in_path
    samtools  bam2fq 
