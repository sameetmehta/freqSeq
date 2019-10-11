# freqSeq: A statistical method to quantify rare DNA lesions at single-base resolution across the genome.

## Structure of the analysis

The method is effectively divided into a number of steps. One of the primary
requirements for this method to work is availability of a high-performance
computing cluster. The basic steps involved in this method are as follows:

 * Determine which read of the fastq read pair has the linker. Designate the
   read having a linker as R1, and the other read as R2.
 * Trim the linker, and any preceding nucleotides. The first base of the R1 now
   denotes the position immediately downstream (3') of the incision nick. In the 
   case of the cyclobutane pyrimidine dimer (CPD), the first base of R1 is the 3'
   base of this two base PyPy lesion and is termed the "+1 base". If no linker
   sequence is found, the read pair is saved to a different file, and not used
   in the current analysis.
 * R1 after this step is always the proximal end (nearest to the lesion), and R2
   is always distal to the lesion.
 * Align the reads to the reference genome.
 * Read through the bam files, and for each R1 alignment re-create the upstream
   10 bases on the correct strand (bases -1 through -10). In the case of CPD,
   base -1 is the 5' base of this two-base PyPy lesion.
 * Save this information to a file, and then divide the file according to the
   class of lesion. For our purposes there are four classes of legions viz.,
   * PuPu -- The first base of R1, and immediate preceding base both are
     purines.
   * PuPy -- The first base of R1 is a pyrimidine and immediate preceding base
     is a purine.
   * PyPu -- The first base of R1 is a purine, and immediate preceding base is a
     pyrimidine.
   * PyPy -- The first base of the R1 is pyrimidine, and the immediate preceding
     base is also a pyrimdine.
 
 * For each of these files (classified according to the lesion dinucleotide) we
   determine whether the given read pair is:
   * Unique mapping -- Exactly one read-pair matches to the given genomic
     location, and that read-pair does not match to any other location in the
     genome.
   * PCR duplicate -- More than one read pair in which R1 and R2 map to the
     exactly same location, and the R2s have exactly same end base. Because
     barcodes are not used, this criterion overestimates PCR duplicates.
   * Recurrent -- More than one R1 match to the same location, but the
     respective R2s match to different locations.
   * deaminated -- The reference base is a T, yet the first base of the aligned
     R1 read is a C. The result is assumed to reflect the known rapid
     deamination of a cytosine within a CPD.

### The steps for Sequence Analysis - I

Here we will describe the exact steps that we use to analyze the data.

 1. First find the linkers, and designate the read with the linker to be R1, and
    the other read to be R2.
 
    This step is accomplished by running the script `findLinkers.py`. This is a
    `Python` script that takes seven arguments. The linker sequence is
    hard-coded in the script. The script finds the linker (searching within
    first few bases of each read), and trims the read to remove everything upto
    the end of the linker sequence
    
 2. Align the trimmed reads.
 
    The trimmed read pair is aligned to the reference genome using `bwa mem`
    aligner.
    
 3. Count locations.
 
    From the aligned files count how many times the reads mapped correctly, and
    also keep counts of all the read pairs that:
    * Did not map.
    * Mapped multiple times.
    * Were of wrong insert size.
    * Were not correctly mapped, or mapped to a non-canonical chromosome, or a
      part of read maps outside of the chromosome boundary.
    * Reconstruct 10 upstream bases for each R1 alignment, and classify each
      alignement according to the type, as PuPu, PuPy, PyPu, or PyPy.
    
    a. Sort and count duplicates
 
       In this step, each class of lesions is sorted, and further classified as
       
       * Single
       * Recurrent
       * PCR duplicate
       
       And saved to individual files.
 
 4. Put Single and Recurrent Files together.
 
After step 4, we have 16 output files for each sample. For these we use the
files that have all the lesions that are in the single+recurrent file for each
of the sample. These files are `SAM` files with a few extra columns. This is a
tab-separated-file, and is quite amenable to analysis using `R`.

### Steps for Analysis - II

The single+recurrent lesions for each sample can be analyzed using `R`. Most the
analysis is possible in a interactive manner.

 
