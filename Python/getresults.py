
import sys
import os
import subprocess

if len(sys.argv) < 2:
    sys.stderr.write("Usage:  getresults.py sampleDir...\n")
    sys.exit(-1)

files = [ "linkerfail.fna", "multimap.sam", "other.sam",
          "PuPu_recurrent.txt", "PuPu_single.txt", "PuPu_duplicate.txt",
          "PuPy_recurrent.txt", "PuPy_single.txt", "PuPy_duplicate.txt",
          "PyPu_recurrent.txt", "PyPu_single.txt", "PyPu_duplicate.txt",
          "PyPy_recurrent.txt", "PyPy_single.txt", "PyPy_duplicate.txt",
          "unmapped.sam", "wronglength.sam", "cdeanimation.sam"
       ]

for dir in sys.argv[1:]:
    if dir.endswith("/"):
        dir = dir[:-1]

    for file in files:
        sys.stdout.write("%s:  %s\n" % (dir, file))
        sys.stdout.flush()

        subprocess.call("gzip -c %s/%s > results/%s_%s.gz" % (dir, file, dir, file), shell=True)    

    fp = open("%s/trimlog.txt" % dir)
    lines = fp.readlines()
    fp.close()
    rawreads = int(lines[0].split(":")[1].strip())
    linkermatch = int (lines[1].split(":")[1].strip()) + int(lines[2].split(":")[1].strip())

    fp = open("%s/countlog.txt" % dir)
    lines = fp.readlines()
    fp.close()

    fp = open("%s/%s_stats.txt" % (dir, dir), "w")
    fp.write("Raw reads          %7d\n" % rawreads)
    fp.write("Linkerfail         %7d\n" % (rawreads - linkermatch))
    for l in lines[2:7]:
        fp.write(l)
    for i in range(11, 23, 3):
        fp.write(lines[i])
    for i in range(12, 24, 3):
        fp.write(lines[i])
    for i in range(13, 25, 3):
        fp.write(lines[i])
    fp.close()

    sys.stdout.write("Raw reads          %7d\n" % rawreads)
    sys.stdout.write("Linkerfail         %7d\n" % (rawreads - linkermatch))
    for l in lines[2:7]:
        sys.stdout.write(l)
    for i in range(11, 23, 3):
        sys.stdout.write(lines[i])
    for i in range(12, 24, 3):
        sys.stdout.write(lines[i])
    for i in range(13, 25, 3):
        sys.stdout.write(lines[i])
