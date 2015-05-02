For reviewers: This site is updated with more instructions to build the index needed for Batalign to run on. Below, you can also find the dependencies and OS which we have tested Batalign on. Thank you for trying this tool.<br><br>

"BatAlign: An incremental method for accurate gapped alignment".
<br>
This is a README file for the usage of Batalign.
<br><br>
Please go to http://compbio.ddns.comp.nus.edu.sg/~limjingq/BATALIGN/
<br>
1) Batalign.tar.gz (the tarball of the program)<br>
2) INPUT_READ_all_one_million_reads_datasets.tbz (the tarball of some 1-million-reads datasets that were used for the manuscript)<br>
3) hg19.BatAlign.index.tbz (the tarball of the reference-hg19-index that can be used to align 2) using 1) )<br>
<br><br>
INSTALL BatAlign<br>
=-=-=-=-=-=-=-=-=<br>
a) Download 1)<br>
b) tar -zxvf 1)<br>
c) Change directory into the top directory of b) "batindel/"<br>
d) Type "./configure" then type "make" and finally "make copy"<br>
e) batalign (the binary of BatAlign) will be created in bin/<br><br>

GETTING pre-built INDEX<br>
=-=-=-=-=-=-=-=-=-=-=-=<br>
a) Download 3)<br>
b) tar -zxvf 3)<br>
c) "hg19.fa" will be the input-string to represent this index into program<br><br>

BUILDING INDEX<br>
=-=-=-=-=-=-=<br>
a) Have a fastq-formatted file ready<br>
b) Locate the script "build_indexX" in "batindel/bin"<br>
c) Type "./build_indexX GENOME.fa" to make the neccessary pairing data-structure based on FM-index.<br><br>

USAGE<br>
=-=-=-= <br>
Single-end-reads <br>
./bin/batalign -g INDEX -q INPUT -o OUTPUT <br><br>

Paired-end-reads <br>
./bin/batalign -g INDEX -q INPUT_left -q INPUT_right -o OUTPUT <br>
Example: ./bin/batalign -g /data/index/hg19/hg19.fa -q CML_R1_left.fq -q CML_R2_right.fq -o CML.out.sam --threads 20 <br><br>

INDEX was the mentioned "hg19.fa". Make sure all index files reside in the same directory.<br>
INPUT can be any fastq/fasta file from 2).<br>
OUTPUT is a SAM-format alignment-mapping file.<br>
BatAlign only needs -g -q -o, the above mentioned parameters, to execute in DEFAULT mode.<br>
To parallelize, use "--threads INT". <br><br>

Built with...<br>
=-=-=-=-=-=-=-<br>
GNU automake v1.11.1, GNU autoconf v2.63, gcc v4.4.7.<br>
Tested on Centos Release 6.6 (Final), Debian GNU/Linux 7<br>
Must be SSE2-instructions compatible<br><br>

Thank you for your patience.<br>
