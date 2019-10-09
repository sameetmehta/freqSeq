
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <vector>
#include <string>
#include <map>

using namespace std;

bool getPairLines(FILE *fp, char *line, char *line2, char *nextline, FILE *otherfp);
void setAlign(char **fields, int chrnum, int& start, int& end, int &alen, char *qalign, char *ralign, char *qscore);
void revcomp(char *seq);
void setChromInfo(char *reffile);
void resetChrSeq(char *chr);
void getSeq(char *chr, int start, int end, char *buf);
int chr2num(char *chr);
void getFlags(char *flagstr, char *buf);
void settabs(char *s, char **fields);

int reflen = 0;
int refnum = -1;
char refseq[300000000];
FILE *reffp;

char *chromseq[1000];
int chromlen[1000];

int main(int argc, char *argv[])
{
	if (argc != 7) {
		fprintf(stderr, "Usage:  countlocations refFile countlog.txt unmapped.sam multimap.sam wronglength.sam other.sam\n");
		exit(-1);
	}

	setChromInfo(argv[1]);

	FILE *log = fopen(argv[2], "w");
	FILE *unmap = fopen(argv[3], "w");
	FILE *multimap = fopen(argv[4], "w");
	FILE *wrong = fopen(argv[5], "w");
	FILE *otherfp = fopen(argv[6], "w");

	int cnt = 0;
	int uniquecnt = 0;
	int unmapcnt = 0;
	int multimapcnt = 0;
	int pupu = 0;
	int pupy = 0;
	int pypu = 0;
	int pypy = 0;
	int deanimated = 0;
	int abasic = 0;
	int other = 0;
	int insertfail = 0;

	FILE *pupufp = fopen("PuPu.txt", "w");
	FILE *pypufp = fopen("PyPu.txt", "w");
	FILE *pupyfp = fopen("PuPy.txt", "w");
	FILE *pypyfp = fopen("PyPy.txt", "w");
	FILE *cdeamfp = fopen("cdeanimation.sam", "w");

	char line[100000], *fields[1000], origline[100000];
	char line2[100000], *fields2[1000], origline2[100000];
	char nextline[100000];
	char qalign[1000], ralign[1000], qscore[1000], origqalign[1000], origralign[1000];
	fgets(nextline, 100000, stdin);
	while (1) { 
		if (!getPairLines(stdin, line, line2, nextline, otherfp)) {
			break;
		}
		if (line[0] == '\0') {
			++other;
			++cnt;
			continue;
		}

		strcpy(origline, line);
		strcpy(origline2, line2);

		settabs(line, fields);
		settabs(line2, fields2);

		int flags = atoi(fields[1]);

		if (flags & 0x100) {
		        ++cnt;
			continue;
		}

		int acclen = strlen(fields[0]);
		if (fields[0][acclen-1] == '1') {
			if (!(flags & 0x40)) {
				continue;
			}
		} else {
			if (!(flags & 0x80)) {
				continue;
			}
		}

		++cnt;
		if (cnt % 50000 == 0) {
			fprintf(stderr, "  -> %d\n", cnt);
			fflush(stderr);
		}

		if (flags & 4) {
			fprintf(unmap, "%s%s", origline, origline2);
			++unmapcnt;
			continue;
		}

		int mapq = atoi(fields[4]);
		if (mapq == 0) {
			fprintf(multimap, "%s%s", origline, origline2);
			++multimapcnt;
			continue;
		}

		int chrnum = chr2num(fields[2]);
		if (chrnum == -1) {
			++other;
			fprintf(otherfp, "%s%s", origline, origline2);
			continue;
		}

		++uniquecnt;

		int insertsize = 0;
		if (sscanf(fields[8], "%d", &insertsize) == 0 ||
			!( (70 <= insertsize && insertsize <= 500) ||
			   (-70 >= insertsize && insertsize >= -500) )) {
			fprintf(wrong, "%s%s", origline, origline2);
			++insertfail;
			continue;
		}
		if (insertsize < 0) {
			insertsize = -insertsize;
		}

		int start = 0;
		int end = 0;
		int alen = 0;
		setAlign(fields, chrnum, start, end, alen, qalign, ralign, qscore);
		strcpy(origqalign, qalign);
		strcpy(origralign, ralign);

		if (flags & 0x10) {
			revcomp(qalign);
			revcomp(ralign);
			int tmp = start;
			start = end;
			end = tmp;
		}

		if (qalign[0] == 'T' && ralign[0] == 'C') {
			++deanimated;
			fprintf(cdeamfp, "%s%s", origline, origline2);
		} else if (qalign[0] == ralign[0]) { 
			char seq[1000];
			char *t = seq;
			char *ref = chromseq[chrnum];
			if (!(flags & 16)) {
				for (int i=-10; i < 10; ++i) {
					if (i == 0) {
						*t++ = '*';
					}
					*t++ = ref[start-1+i];
				}
				*t = '\0';
			} else {
				for (int i=-9; i <= 10; ++i) {
					*t++ = ref[start-1+i];
					if (i == 0) {
						*t++ = '*';
					}
				}
				*t = '\0';
				revcomp(seq);
			}
			t = &seq[11];
			for (int i=0; i < alen; ++i) {
				if (qalign[i] != '-') {
					*t++ = qalign[i];
				}
			}
			*t = '\0';
			
			int pym1flag = (seq[9] == 'C' || seq[9] == 'T');
			int pyp1flag = (seq[11] == 'C' || seq[11] == 'T');

			FILE *fp = NULL;
			if (pym1flag && pyp1flag) {
				fp = pypyfp;
				++pypy;
			} else if (pym1flag && !pyp1flag) {
				fp = pypufp;
				++pypu;
			} else if (!pym1flag && pyp1flag) {
				fp = pupyfp;
				++pupy;
			} else {
				fp = pupufp;
				++pupu;
			}

			fprintf(fp, "%s\t%d\t%s\t%s\t%d", fields[2], start, ((flags & 16) ? "-" : "+"), seq, insertsize);

			char flagstr[1000];
			getFlags(fields[1], flagstr);
			char *mapq = fields[4];
			char *cigar = fields[5];
			char *mdstr = "";
			char *xastr = "-";
			for (int x=11; fields[x] != NULL; ++x) {
				if (strncmp(fields[x], "MD:Z:", 5) == 0) {
					mdstr = fields[x];
				} else if (strncmp(fields[x], "XA:Z:", 5) == 0) {
					xastr = fields[x];
				}
			}

			fprintf(fp, "\t%s\t%s\t%s\t%s", flagstr, mapq, cigar, mdstr);

			char *otherpos = fields2[3];
			char *otherseq = fields2[9];
			char otherflags[1000];
			getFlags(fields2[1], otherflags);
			char *othermapq = fields2[4];
			char *othercigar = fields2[5];
			char *othermdstr = "";
			char *otherxastr = "-";
			for (int x=11; fields2[x] != NULL; ++x) {
				if (strncmp(fields2[x], "MD:Z:", 5) == 0) {
					othermdstr = fields2[x];
				} else if (strncmp(fields2[x], "XA:Z:", 5) == 0) {
					otherxastr = fields2[x];
				}
			}
			
			fprintf(fp, "\t%s\t%s\t%s\t%s\t%s\t%s", otherpos, otherseq, otherflags, othermapq, othercigar, othermdstr);

			fprintf(fp, "\t%s\t%s\n", xastr, otherxastr);

		} else {
			/*
			printf("%s", origline);
			printf("   %s\n   %s\n", origqalign, origralign);
			printf("%s:%d..%d\n   %s\n   %s\n", fields[2], start, end, qalign, ralign);
			*/
			fprintf(otherfp, "%s%s", origline, origline2);
			++other;
		}
	}
	fclose(unmap);
	fclose(multimap);
	fclose(wrong);
	fclose(otherfp);

	fprintf(stderr, "%d reads, %d unmapped, %d multimap, %d wronglength, %d other, %d deanimated, %d pypu, %d pupu, %d pupy, %d pypy.\n", cnt, unmapcnt, multimapcnt, insertfail, other, deanimated, pypu, pupu, pupy, pypy);

	fprintf(log, "Alignment statistics:\n");
	fprintf(log, "Linker Match       %7d\n", cnt);
	fprintf(log, "Unmapped           %7d\n", unmapcnt);
	fprintf(log, "Multimap           %7d\n", multimapcnt);
	fprintf(log, "Wronglength        %7d\n", insertfail);
	fprintf(log, "Other              %7d\n", other);
	fprintf(log, "Deaminated         %7d\n", deanimated);
	fprintf(log, "PyPu               %7d\n", pypu);
	fprintf(log, "PuPu               %7d\n", pupu);
	fprintf(log, "PuPy               %7d\n", pupy);
	fprintf(log, "PyPy               %7d\n", pypy);
	fclose(log);

	return 0;
}

bool getPairLines(FILE *fp, char *line, char *line2, char *nextline, FILE *otherfp)
{
	if (nextline[0] == '\0') {
		return false;
	}

	char accno[1000];
	char *t = accno;
	for (char *s=nextline; *s && !isspace(*s); s++,t++) *t = *s;
	*t = '\0';
	int accnolen = strlen(accno);
	bool switchFlag = (accno[accnolen-1] == '2');

	strcpy((!switchFlag ? line : line2), nextline);

	if (!fgets(nextline, 100000, fp)) {
		fprintf(stderr, "%s", nextline);
		fprintf(stderr, "Error:  Premature end of file.\n");
		exit(-1);
	}
	if (strncmp(nextline, accno, accnolen) != 0) {
		fprintf(stderr, "Error:  Invalid paired-end alignments.\n");
		exit(-1);
	}
	strcpy((!switchFlag ? line2 : line), nextline);

	if (fgets(nextline, 100000, fp)) {
		if (strncmp(nextline, accno, accnolen) == 0) {
			fprintf(otherfp, "%s%s%s", line, line2, nextline);
			line[0] = '\0';
			line2[0] = '\0';
			nextline[0] == '\0';
			while (fgets(nextline, 100000, fp)) {
				if (strncmp(nextline, accno, accnolen) != 0) {
					break;
				}
				fprintf(otherfp, "%s", nextline);
				nextline[0] = '\0';
			}

			return true;
		}
	} else {
		nextline[0] = '\0';
	}

	return true;
}

void setAlign(char **fields, int chrnum, int& start, int& end, int &alen, char *qalign, char *ralign, char *qscore)
{
	start = atoi(fields[3]);
	end = start - 1;
	alen = 0;

	char *readp = fields[9];
	char *scorep = fields[10];

	char *cigar = fields[5];

	char *ref = chromseq[chrnum];
	int refpos = start - 1;
	for (char *s=cigar; *s; ++s) {
		int num = 0;
		for ( ; *s && isdigit(*s); ++s) {
			num = num * 10 + (*s - '0');
		}
		if (*s == 'M') {
			for (int i=0; i < num; ++i) {
				qalign[alen] = *readp++;
				qscore[alen] = *scorep++;
				ralign[alen] = ref[refpos++];
				++end;
				++alen;
			}
		} else if (*s == 'D') {
			for (int i=0; i < num; ++i) {
				qalign[alen] = '-';
				qscore[alen] = '\0';
				ralign[alen] = ref[refpos++];
				++end;
				++alen;
			}
		} else if (*s == 'I') { 
			for (int i=0; i < num; ++i) {
				qalign[alen] = *readp++;
				qscore[alen] = *scorep++;
				ralign[alen] = '-';
				++alen;
			}
		} else if (*s == 'S') {
			for (int i=0; i < num; ++i) {
				readp++;
				scorep++;
			}
		} else if (*s == 'H') {
		} else {
			fprintf(stderr, "Error:  Invalid cigar string:  %s\n", cigar);
			exit(-1);
		}
	}
	qalign[alen] = '\0';
	ralign[alen] = '\0';
	qscore[alen] = '\0';
	
	int last = -1;
	int lastscore = 0;
	for (int i=0; i < alen; ++i) {
		if (qscore[i] == '\0') {
			continue;
		}
		if (last + 1 < i) {
			int ls = (int) lastscore - 33;
			int ns = (int) qscore[i] - 33;
			if (lastscore == 0) {
				ls = ns;
			}

			int newscore = (ls + ns) / 2;
			char score = (char) (newscore + 33);
			while (last + 1 < i) {
				qscore[last+1] = score;
				++last;
			}
		}
		last = i;
	}
	while (last + 1 < alen) {
		qscore[last+1] = qscore[last];
		++last;
	}
}

char rc(char ch)
{
	switch (ch) {
	case 'A':  return 'T';
	case 'C':  return 'G';
	case 'G':  return 'C';
	case 'T':  return 'A';
	case '-':  return '-';
	case '*':  return '*';
	default:   return 'N';
	}
}

void revcomp(char *seq)
{
	int len = strlen(seq);
	for (int i=0,j=len-1; i <= j; ++i,--j) {
		if (i == j) {
			seq[i] = rc(seq[i]);
		} else {
			char tmp = rc(seq[i]);
			seq[i] = rc(seq[j]);
			seq[j] = tmp;
		}
	}
}

typedef struct {
	string chr;
	off_t offset;
	int seqlen;
} CHR_INFO;

vector<CHR_INFO> chromlist;

void getSeq(char *chr, int start, int end, char *buf)
{
	int chrnum = chr2num(chr);
	if (chrnum == -1) {
		buf[0] = '\0';
		return;
	}

	char *ref = chromseq[chrnum];
	char *s = buf;
	for (int i=start; i < end; ++i) {
		*s++ = ref[i-1];
	}
	*s = '\0';
}

void setChromInfo(char *reffile)
{
	reffp = fopen64(reffile, "r");
	if (reffp == NULL) {
		fprintf(stderr, "Error:  Cannot open reference file:  %s\n", reffile);
		exit(-1);
	}

	char faifile[10000];
	sprintf(faifile, "%s.fai", reffile);
	FILE *fp = fopen(faifile, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error:  Cannot open .fai file:  %s\n", faifile);
		exit(-1);
	}
	char line[10000];
	char *fields[1000];
	while (fgets(line, 10000, fp)) {
		settabs(line, fields);
		char *s = fields[0];
		for ( ; *s && !isspace(*s); ++s) ;
		*s = '\0';

		CHR_INFO info;
		info.chr = fields[0];
		info.seqlen = atoi(fields[1]);
		info.offset = (uint64_t) 0;
		for (char *s=fields[2]; *s; ++s) {
			info.offset = (info.offset * (off_t) 10) + (off_t) (*s - '0');
		}
		chromlist.push_back(info);
	}
	fclose(fp);

	for (int i=0; i < (int) chromlist.size(); ++i) {
		char *chr = (char *) chromlist[i].chr.c_str();
		int chrnum = chr2num(chr);
		if (chrnum == -1) {
			continue;
		}

		fprintf(stderr, "Reading %s...\n", chr);
		fflush(stderr);

		resetChrSeq(chr);
		chromlen[chrnum] = reflen;
		chromseq[chrnum] = strdup(refseq);
	}
}

void resetChrSeq(char *chr)
{
	int chrnum = chr2num(chr);
	if (chrnum == refnum) {
		return;
	}

	if (chrnum == -1) {
		return;
	}

	string c = chr;
	refnum = chrnum;

	int idx = 0;
	for ( ; idx < chromlist.size(); ++idx) {
		if (c == chromlist[idx].chr) {
			break;
		}
	}
	if (idx == chromlist.size()) {
		fprintf(stderr, "Error:  Chromosome not found in reference file:  %s\n", chr);
		exit(-1);
	}

	fseeko(reffp, chromlist[idx].offset, SEEK_SET);

	reflen = 0;
	char line[10000];
	while (reflen < chromlist[idx].seqlen) {
		if (!fgets(line, 10000, reffp)) {
			fprintf(stderr, "Error:  Reach EOF of reference file when reading chromosome:  %s\n", chr);
			exit(-1);
		}
		for (char *s=line; *s != '\n'; ++s) {
			refseq[reflen++] = toupper(*s);
		}
	}
	if (reflen != chromlist[idx].seqlen) {
		fprintf(stderr, "Error:  Could not ref chromosome from reference file:  %s\n", chr);
		exit(-1);
	}
}
	
int chr2num(char *chr)
{
	int idx = (chr[0] == 'c' ? 3 : 0);
	if (isdigit(chr[idx])) {
		return atoi(chr+idx);
	} else if (chr[idx] == 'M' && (chr[idx+1] == '\0' || (chr[idx+1] == 'T' && chr[idx+2] == '\0'))) {
		return 23;
	} else if (chr[idx] == 'X' && chr[idx+1] == '\0') {
		return 24;
	} else if (chr[idx] == 'Y' && chr[idx+1] == '\0') {
		return 25;
	} else {
		return -1;
	}
}

void getFlags(char *flagstr, char *buf)
{
	int flags = atoi(flagstr);

	int pos = 0;
	for (int i=1; i < 0x1000; i<<=1) {
		if (flags & i) {
			sprintf(buf+pos, "0x%x,", i);
			pos += strlen(buf+pos);
		}
	}
	if (pos > 0) {
		--pos;
	}
	buf[pos++] = '\0';
}

void settabs(char *s, char **fields)
{
	int i = 0;

	fields[i] = s;
	for ( ; *s && *s != '\n'; ++s) {
		if (*s == '\t') {
			++i;
			fields[i] = s+1;
			*s = '\0';
		}
	}
	++i;
	fields[i] = NULL;
	*s = '\0';
}
