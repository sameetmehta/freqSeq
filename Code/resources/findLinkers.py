import sys
import gzip

#linker = "GATCTGAATTC"
#linker = "GANNNGGGAGATCTGAATTCT"
linker = "GGGAGATCTGAATTCT" # for SP04, 05, and 06, SP12L16, AND SP13L16
# linker = "GGCCGCGATCTGAATTCT" # for SP12L17, and SP13L17
# linker = "GATCTGAATTC" # for SP02 AND 03 # See email from Sanjay sent on May 27

if len(sys.argv) != 7:
    sys.stderr.write("Usage:  findLinkers read1.gz read2.gz trim1.gz trim2.gz trimlog.txt linkerfail.txt\n")
    sys.exit(-1)

def checkLinker(s):
    # We will keep the details in a dictionary.
    d = {3:{'cnt' : 0, 'last' : 0}, \
         4:{'cnt' : 0, 'last' : 0}, \
         0:{'cnt' : 0, 'last' : 0}, \
         1:{'cnt' : 0, 'last' : 0}, \
         2:{'cnt' : 0, 'last' : 0}, \
         5:{'cnt' : 0, 'last' : 0}, \
	 6:{'cnt' : 0, 'last' : 0}, \
	 7:{'cnt' : 0, 'last' : 0}} # We are going to try 8 offsets
    
    ks = d.keys()
    ks.sort() # we will go serially

    for k in ks: # for each key, this is the offset on the sequence
        for i in range(len(linker)):
	    if s[k+i] != linker[i]:
	        d[k]['cnt'] += 1
		d[k]['last'] = i

    counts = [d[k]['cnt'] for k in ks] # gather all the counts
    myind  = counts.index(min(counts)) # get the index of minimum counts,
                                       # basically information on the offset
				       # where minimum errors were seen
    myk = ks[myind] # identify the key with least errors
    return myk, d[myk]['cnt'], d[myk]['last'] # return the offset and other
                                              # numbers

read1cnt = 0
read2cnt = 0
nomatch = 0
readcnt = 0

f1 = gzip.open(sys.argv[1], "rb")
f2 = gzip.open(sys.argv[2], "rb")
out1 = gzip.open(sys.argv[3], "wb")
out2 = gzip.open(sys.argv[4], "wb")
log = open(sys.argv[5], "w")
failfp = open(sys.argv[6], "w")

while True: # each fastq record has 4 lines
    l11 = f1.readline() # read first line of first fastq this is the identifier
    l12 = f1.readline() # read second line of first fastq this is the sequence
    l13 = f1.readline() # read third line of first fastq this the "+"
    l14 = f1.readline() # read fourth line of first fastq this is the quality score

    l21 = f2.readline() # read first line of second fastq this is the identifier
    l22 = f2.readline() # read second line of second fastq this is the identifier
    l23 = f2.readline() # read third line of second fastq this is the identifier
    l24 = f2.readline() # read fourth line of second fastq this is the identifier
    
    if not l11 or not l21:
        if not l11 and not l21:
            break
        sys.stderr.write("Error:  Input files do not match.\n")
        sys.exit(-1)

    if not l12 or not l13 or not l14 or not l22 or not l23 or not l24:
        sys.stderr.write("Error:  Not a fastq file.\n")
        sys.exit(-1)

    if not l11.startswith("@") or not l13.startswith("+") or \
       not l21.startswith("@") or not l23.startswith("+") or \
       len(l12) != len(l14) or len(l22) != len(l24):
        sys.stderr.write("Error:  Not a fastq file.\n")
        sys.exit(-1)
    # sanity checks to make sure that we are dealing with fastq files.

    readcnt += 1
    if readcnt % 20000 == 0:
        sys.stderr.write("   -> %d\n" % readcnt)
        sys.stderr.flush() # print some progress to the stderr

    (pref,suf) = l11[:-1].split(" ", 1)

##     cnt1 = 0 # initialize counts for first fastq file
##     last1 = 0
##     for i in range(len(linker)):
##         if l12[5+i] != linker[i]:
##             cnt1 += 1
##             last1 = i
    
    offset1, cnt1, last1 = checkLinker(l12)


##     cnt2 = 0 # initialize counts for second fastq file
##     last2 = 0
##     for i in range(len(linker)):
##         if l22[5+i] != linker[i]:
##             cnt2 += 1
##             last2 = i
    
    offset2, cnt2, last2 = checkLinker(l22)

    if cnt1 <= 2 and last1 < (len(linker) - 5) + offset1:
#         ### for debugging
# 	print l12 
# 	print offset1 * " " + linker, offset1
# 	print l12[len(linker)+offset1:]
# 	raw_input()
# 	###
        read1cnt += 1
        out1.write(pref + "_A1 " + suf + "\n")
        out1.write(l12[len(linker)+offset1:]) # this change is specific for what we
	                                      # doing on april 14 2016
        out1.write(l13)
        out1.write(l14[len(linker)+offset1:])
        out2.write(pref + "_A1 " + suf + "\n")
        out2.write(l22[10:])
        out2.write(l23)
        out2.write(l24[10:])
    elif cnt2 <= 2 and last2 < (len(linker) - 5) + offset2:
#         ### for debugging
# 	print l22 
# 	print offset2 * " " + linker, offset2
# 	print l22[len(linker)+offset2:]
# 	raw_input()
# 	###
        read2cnt += 1
        out1.write(pref + "_A2 " + suf + "\n")
        out1.write(l12[10:])
        out1.write(l13)
        out1.write(l14[10:])
        out2.write(pref + "_A2 " + suf + "\n")
        out2.write(l22[len(linker)+offset2:])
        out2.write(l23)
        out2.write(l24[len(linker)+offset2:])
    else:
        nomatch += 1
        failfp.write(">" + l11[1:])
        failfp.write(l12)
        failfp.write(">" + l21[1:])
        failfp.write(l22)

f1.close()
f2.close()
out1.close()
out2.close()
failfp.close()

print "Num reads:", readcnt
print "R1 Linker:", read1cnt
print "R2 Linker:", read2cnt
print "No match: ", nomatch

log.write("Num reads:\t%d\n" % readcnt)
log.write("R1 Linker:\t%d\n" % read1cnt)
log.write("R2 Linker:\t%d\n" % read2cnt)
log.write("No match:\t%d\n" % nomatch)
log.close()
