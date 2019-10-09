
import sys

if len(sys.argv) != 6:
    sys.stderr.write("Usage:  dupcnt.py sorted.txt single.txt recurrent.txt pcrdup.txt name\n")
    sys.exit(-1)

fp = open(sys.argv[1])
singlefp = open(sys.argv[2], "w")
recurfp = open(sys.argv[3], "w")
dupfp = open(sys.argv[4], "w")
name = sys.argv[5]

singlecnt = 0
recurcnt = 0
dupcnt = 0

nextline = fp.readline()
while nextline:
    lines = []
    fields = nextline.strip().split("\t")
    chr = fields[0]
    pos = fields[1]
    strand = fields[2]
    insertsize = fields[4]

    sizes = {}
    while True:
        lines.append((insertsize, nextline))
        if insertsize not in sizes:
            sizes[insertsize] = 0
        sizes[insertsize] += 1

        nextline = fp.readline()
        if not nextline:
            break

        f = nextline.strip().split("\t")
        c = f[0]
        p = f[1]
        s = f[2]
        i = f[4]
        if not (chr == c and pos == p and strand == s):
            break

        insertsize = i

    for s in sizes:
        if sizes[s] > 1:
	    dupcnt -= 1 # to account for the duplicates being written twice,
	                # once in the single or the recurrent file, and one in
			# the duplicate file
            for t in lines:
                if t[0] == s:
                    dupfp.write(t[1])
                    dupcnt += 1

    if len(sizes) == 1:
        singlefp.write(lines[0][1])
        singlecnt += 1
    else:
        done = set()
        for t in lines:
            if t[0] in done:
                continue

            recurfp.write(t[1])
            done.add(t[0])
            recurcnt += 1

fp.close()
singlefp.close()
recurfp.close()

sys.stdout.write("Duplicate    %s  %7d\n" % (name, dupcnt))
sys.stdout.write("Recurrent    %s  %7d\n" % (name, recurcnt))
sys.stdout.write("Single       %s  %7d\n" % (name, singlecnt))
