# freqseq: A statistical method to quantify rare lesions at single-base resolution across the genome.

## Structure of the analysis

The method is effectively divided into a number of steps. One of the primary
requirements for this method to work is availablity of a high-performance
computing cluster. The basic steps involved in this method are as follows:

 * Deterine which read for the pared fastq reads has the linker.
 * Trim the linker, and any preceding nucleotides. The first base of the R1 now
   denotes the position immediately downstream of the lesion. If no linker
   sequence is found, the read pair is saved to a different file, and not used
   in the current analysis.
 * Designate the read with linker as R1, and the other as R2. R1 after this step
   is always the proximal end (nearest to the lesion), and R2 is always distal
   to the legion.
 * Align the reads to the reference genome.
 * Read through the bam files, and for all the R1 alignments re-create upstream
   10 bases on the correct strand.
 * Save this information to a file, and then divide the file by class of lesion.
   For our purposes there are four classes of legions viz.,
   * PuPu -- The first base of R1, and immediate preceding base both are
     purines.
   * PuPy -- The first base of R1 is a purine but immediate preceding base is a
     pyrimidine.
   * PyPu -- The first base of R1 is a pyrimidine, and immediate preceding base
     is a purine.
   * PyPy -- The first base of the R1 is pyrimidine, and the immediate preceding
     base is also a pyrimidne.
 
 * For each of these files we determine of the given read pair is
   * Unique mapping -- Exactly one read-pair matches to the given genomic
     location.
   * PCR duplicate -- If both R1 and R2 map to the exactly same location.
   * Recurrent -- If more than one R1 match to the same location, but R2 match
     to the different locations.
   * deaminated -- If the reference is a T, and first base of the aligned R1
     read is a C, then we call it de-aminated base.

## The steps

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
    * Mappted multiple times.
    * Were of wrong insert size.
    * Were not correctly mapped but were of other type.
    * Reconstruct 10 upstream bases for each R1 alignment, and classify each
      alignement according to the type, as PuPu, PuPy, PyPu, or PyPy.
    
    a. Sort and count duplicates
 
       In this step, each class of lesions is sorted, and further classified as
       
       * Single 
       * Recurrent
       * PCR duplicatea
       
       And saved to individual files.
 
 4. Put Single and Recurrent Files together.
 
After step 4, we have 16 output files for each sample.

 
