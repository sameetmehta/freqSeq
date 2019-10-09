
shopt -s expand_aliases
set -o errexit -o nounset -o pipefail

alias pbm="source /ycga-ba/home/bioinfo/software/knightlab/bin/pbmscript.sh"

if [ ! -e reads1.fastq.gz ] ; then
  echo Concatenating R1 reads...
  cat Unaligned/*_R1_*.fastq.gz > reads1.fastq.gz

  echo Concatenating R2 reads...
  cat Unaligned/*_R2_*.fastq.gz > reads2.fastq.gz
fi

if [ ! -e trim1.fastq.gz ] ; then
  echo Trimming linkers...

  python ../findLinkers.py reads1.fastq.gz reads2.fastq.gz trim1.fastq.gz trim2.fastq.gz trimlog.txt linkerfail.fna
fi

pbm hs37d5bundle
pbm bwa
pbm samtools

if [ ! -e aligned.bam ] ; then
  echo Aligning reads...

  bwa mem -M -t 8 $REF trim1.fastq.gz trim2.fastq.gz | samtools view -hS - > aligned.bam
fi

echo Counting locations...

samtools view aligned.bam | ../countlocations $REF countlog.txt unmapped.sam multimap.sam wronglength.sam other.sam

echo Making final lists...

python ../msort.py < PyPu.txt > PyPu_sorted.txt
python ../dupcnt.py PyPu_sorted.txt PyPu_single.txt PyPu_recurrent.txt PyPu_duplicate.txt PyPu | tee -a countlog.txt
python ../msort.py < PuPu.txt > PuPu_sorted.txt
python ../dupcnt.py PuPu_sorted.txt PuPu_single.txt PuPu_recurrent.txt PuPu_duplicate.txt PuPu | tee -a countlog.txt
python ../msort.py < PuPy.txt > PuPy_sorted.txt
python ../dupcnt.py PuPy_sorted.txt PuPy_single.txt PuPy_recurrent.txt PuPy_duplicate.txt PuPy | tee -a countlog.txt
python ../msort.py < PyPy.txt > PyPy_sorted.txt
python ../dupcnt.py PyPy_sorted.txt PyPy_single.txt PyPy_recurrent.txt PyPy_duplicate.txt PyPy | tee -a countlog.txt


