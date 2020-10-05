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
    # 
