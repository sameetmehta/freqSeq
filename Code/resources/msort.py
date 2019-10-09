
import sys

lastlist = []
while True:
    line = sys.stdin.readline()
    if not line:
        break

    line = line.strip()
    fields = line.split("\t")
    
    if "_" in fields[0]: continue

    mychr = fields[0].replace("chr", "")

    if mychr.isdigit():
        chrnum = int(mychr)
    elif mychr[0] == 'M':
        chrnum = 23
    elif mychr[0] == 'X':
        chrnum = 24
    elif mychr[0] == 'Y':
        chrnum = 25
    else:
        sys.stderr.write("Error:  Invalid chromosome:  " + fields[0] + "\n")
        sys.exit(-1)

    lastlist.append((chrnum, int(fields[1]), line))

lastlist.sort()
for t in lastlist:
    sys.stdout.write("%s\n" % (t[2],))
