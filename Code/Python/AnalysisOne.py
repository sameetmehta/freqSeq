
import gzip
import urllib2

import re
import time
import pygal
from cruzdb import Genome, sequence
from operator import itemgetter
from bisect import bisect_left, bisect_right

from datetime import datetime

class SortedCollection(object):
    '''Sequence sorted by a key function.

    SortedCollection() is much easier to work with than using bisect() directly.
    It supports key functions like those use in sorted(), min(), and max().
    The result of the key function call is saved so that keys can be searched
    efficiently.

    Instead of returning an insertion-point which can be hard to interpret, the
    five find-methods return a specific item in the sequence. They can scan for
    exact matches, the last item less-than-or-equal to a key, or the first item
    greater-than-or-equal to a key.

    Once found, an item's ordinal position can be located with the index() method.
    New items can be added with the insert() and insert_right() methods.
    Old items can be deleted with the remove() method.

    The usual sequence methods are provided to support indexing, slicing,
    length lookup, clearing, copying, forward and reverse iteration, contains
    checking, item counts, item removal, and a nice looking repr.

    Finding and indexing are O(log n) operations while iteration and insertion
    are O(n).  The initial sort is O(n log n).

    The key function is stored in the 'key' attibute for easy introspection or
    so that you can assign a new key function (triggering an automatic re-sort).

    In short, the class was designed to handle all of the common use cases for
    bisect but with a simpler API and support for key functions.

    >>> from pprint import pprint
    >>> from operator import itemgetter

    >>> s = SortedCollection(key=itemgetter(2))
    >>> for record in [
    ...         ('roger', 'young', 30),
    ...         ('angela', 'jones', 28),
    ...         ('bill', 'smith', 22),
    ...         ('david', 'thomas', 32)]:
    ...     s.insert(record)

    >>> pprint(list(s))         # show records sorted by age
    [('bill', 'smith', 22),
     ('angela', 'jones', 28),
     ('roger', 'young', 30),
     ('david', 'thomas', 32)]

    >>> s.find_le(29)           # find oldest person aged 29 or younger
    ('angela', 'jones', 28)
    >>> s.find_lt(28)           # find oldest person under 28
    ('bill', 'smith', 22)
    >>> s.find_gt(28)           # find youngest person over 28
    ('roger', 'young', 30)

    >>> r = s.find_ge(32)       # find youngest person aged 32 or older
    >>> s.index(r)              # get the index of their record
    3
    >>> s[3]                    # fetch the record at that index
    ('david', 'thomas', 32)

    >>> s.key = itemgetter(0)   # now sort by first name
    >>> pprint(list(s))
    [('angela', 'jones', 28),
     ('bill', 'smith', 22),
     ('david', 'thomas', 32),
     ('roger', 'young', 30)]

    '''

    def __init__(self, iterable=(), key=None):
        self._given_key = key
        key = (lambda x: x) if key is None else key
        decorated = sorted((key(item), item) for item in iterable)
        self._keys = [k for k, item in decorated]
        self._items = [item for k, item in decorated]
        self._key = key

    def _getkey(self):
        return self._key

    def _setkey(self, key):
        if key is not self._key:
            self.__init__(self._items, key=key)

    def _delkey(self):
        self._setkey(None)

    key = property(_getkey, _setkey, _delkey, 'key function')

    def clear(self):
        self.__init__([], self._key)

    def copy(self):
        return self.__class__(self, self._key)

    def __len__(self):
        return len(self._items)

    def __getitem__(self, i):
        return self._items[i]

    def __iter__(self):
        return iter(self._items)

    def __reversed__(self):
        return reversed(self._items)

    def __repr__(self):
        return '%s(%r, key=%s)' % (
            self.__class__.__name__,
            self._items,
            getattr(self._given_key, '__name__', repr(self._given_key))
        )

    def __reduce__(self):
        return self.__class__, (self._items, self._given_key)

    def __contains__(self, item):
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return item in self._items[i:j]

    def index(self, item):
        'Find the position of an item.  Raise ValueError if not found.'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].index(item) + i

    def count(self, item):
        'Return number of occurrences of item'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].count(item)

    def insert(self, item):
        'Insert a new item.  If equal keys are found, add to the left'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def insert_right(self, item):
        'Insert a new item.  If equal keys are found, add to the right'
        k = self._key(item)
        i = bisect_right(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def remove(self, item):
        'Remove first occurence of item.  Raise ValueError if not found'
        i = self.index(item)
        del self._keys[i]
        del self._items[i]

    def find(self, k):
        'Return first item with a key == k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i != len(self) and self._keys[i] == k:
            return self._items[i]
        raise ValueError('No item found with key equal to: %r' % (k,))

    def find_le(self, k):
        'Return last item with a key <= k.  Raise ValueError if not found.'
        i = bisect_right(self._keys, k)
        if i:
            return self._items[i-1]
        raise ValueError('No item found with key at or below: %r' % (k,))

    def find_lt(self, k):
        'Return last item with a key < k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i:
            return self._items[i-1]
        raise ValueError('No item found with key below: %r' % (k,))

    def find_ge(self, k):
        'Return first item with a key >= equal to k.  Raise ValueError if not found'
        i = bisect_left(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key at or above: %r' % (k,))

    def find_gt(self, k):
        'Return first item with a key > k.  Raise ValueError if not found'
        i = bisect_right(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key above: %r' % (k,))

def rc (my_sequence):
    my_dictionary = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([my_dictionary[base] for base in reversed(my_sequence)])
def getSeq(s):
    a=1
def detError(ref,seq):
    #print seq
    err=0
    add=-1
    for i in  range (11,len(seq)):
        add=add+1
        #print seq[i]+" "+ref[add]
        #print ref
        #print seq
        if seq[i] <> ref[add]:
            err=err+1
    return err
def makeChart():
    line_chart = pygal.Line()
    line_chart.title = 'Browser usage evolution (in %)'
    line_chart.x_labels = map(str, range(2002, 2013))
    line_chart.add('Firefox', [None, None,    0, 16.6,   25,   31, 36.4, 45.5, 46.3, 42.8, 37.1])
    line_chart.add('Chrome',  [None, None, None, None, None, None,    0,  3.9, 10.8, 23.8, 35.3])
    line_chart.add('IE',      [85.8, 84.6, 84.7, 74.5,   66, 58.6, 54.7, 44.8, 36.2, 26.6, 20.1])
    line_chart.add('Others',  [14.2, 15.4, 15.3,  8.9,    9, 10.4,  8.9,  5.8,  6.7,  6.8,  7.5])
    line_chart.render_to_file('line_chart.svg')
    bar_chart = pygal.Bar()                                            # Then create a bar graph object
    bar_chart.add('Fibonacci', [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55])  # Add some values
    bar_chart.render_to_file('bar_chart.svg')    


def slidingWindowII(regions,DamagePos,file):
    f = open(file+".count",'w')
    f.write("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\"\n")
    for chrom in regions.keys():
        if chrom in DamagePos: 
            #print chrom
            mapp=load_map(chrom)
            #print "hhh"
            DamagePosTemp=DamagePos[chrom]
                        
            for start in sorted(regions[chrom].keys()):
              for end in regions[chrom][start].keys(): 
                #print chrom +" "+start.__str__()+" "+end.__str__()
                 
                 for x in range (start,end,100):
                    #print x
                    
                    add=0 
                    adds=""
                    try:
                        # quite smart - the condition asks whether the dam pos < (x+100) is GREATER/EQUAL x!
                        #dp=DamagePosTemp.find_lt(x+99)[1]
                        if x <= DamagePosTemp.find_lt(x+100)[1]:
                            #print x.__str__()+" "+dp.__str__()
                            
                            
                            for xx in range (x,x+100):
                                #print xx                                  
                                
                                try:
                                    count=DamagePosTemp.find(xx)[0]
                                    #print "st"
                                    #print x
                                    #print xx
                                    #print DamagePosTemp.find(xx)[0]
                                    #print DamagePosTemp.find(xx)[1]
                                    #print "en"
                                    add=add+count
                                    adds=adds+"("+xx.__str__()+","+count.__str__()+"),"
                                    
                                except:
                                    a=1
                    except:
                        a=1
                 
                    if add>0:
                        #print "i"
                        try:
                            #print "i"
                            #print(mapp.find_lt(x+99))
                            #print x
                            #print "ii"
                            if x <= mapp.find_lt(x+99)[1]:
                               #print x.__str__()+" "+mapp.find_lt(x+99)[0].__str__()
                               ind=mapp.index(mapp.find_lt(x+99))
                               tmappability=(mapp[ind-1][0]+mapp[ind][0]+mapp[ind+1][0])/3 
                               #print tmappability
                            else:
                                tmappability=0
                        except:
                            a=1
                        a=1
                        
                        f.write(chrom+"\t"+x.__str__()+"\t"+(x+99).__str__()+"\t"+add.__str__()+"\t"+tmappability.__str__()+"\t"+adds+"\n" )  
                
        #print chrom                       
 
def load_map(chrom):
    #print chrom
    coll=SortedCollection(key=itemgetter(1))
    with open("/Users/mok6/Dropbox/LiClipse/MelArray/mappabilityII/"+chrom+".txt",'rb') as r:             
        #add=0
        for row in r:
          #add=add+1
          #print add  
          items=row.split()
          #print item[3]
          #print item
          #we are not interested in zero mappability
          if (float(items[3])>0):
           coll.insert((float(items[3]),int(items[1])))
    return coll
 
def rowPass(row,flag):
  
    all=row.split()
    
    chr1=all[0]
    pos1=all[1]
    strand1=all[2]
    seq1=all[3]
    spacer=all[4]
    hex1=all[5]
    map1=all[6]
    cigar1=all[7]
    mdz1=all[8]
    pos2=all[9]
    seq2=all[10]
    hex2=all[11]
    map2=all[12]
    cigar2=all[13]
    mdz2=all[14]
    xa1=all[15]
    xa2=all[16]
    
    
    
    if ((mdz1 == "MD:Z:55") & (mdz2 == "MD:Z:66") & (xa1[0] != "X") & (xa2[0] !="X") ):
      #print mdz1
      ##print mdz2
      #print xa1
      #print xa2
      #print pos2  
    
      return "1"
      print row
    else:
      #return "1" 
      return flag
                
def loadDamagePosII(file,flag):
    a={}
    b={}
    
    file_r=file+"_single.txt.gz"
    with gzip.open(file_r) as r:
        for row in r:
         if row[0] == "M":
             row="M"+row[2:len(row)]    
         row="chr"+row   
         if rowPass(row,flag)=="1":   
          
          items=row.split()
          
          #if items[0] == "chrX":
          #   print items[1]
          
          if items[0] not in b:
              b[items[0]]={}
          if items[1] not in b[items[0]]:
              #print items[0]+"   "+items[1]
              b[items[0]][items[1]]=1
          else:
              b[items[0]][items[1]]=b[items[0]][items[1]]+1
     
    file_r=file+"_recurrent.txt.gz"
    with gzip.open(file_r) as r:
        for row in r:
            #strange labeling of MT chromosome
         if row[0] == "M":
             row="M"+row[2:len(row)]  
         row="chr"+row 
         if rowPass(row,flag)=="1":
          #print "i"
          items=row.split()
          if items[0] not in b:
              b[items[0]]={}
          if items[1] not in b[items[0]]:
              b[items[0]][items[1]]=1
          else:
              b[items[0]][items[1]]=b[items[0]][items[1]]+1
     
          
    file_r=file+"_single.txt.gz"      
    with gzip.open(file_r,'rb') as r:
        for row in r:
          if row[0] == "M":
             row="M"+row[2:len(row)]   
          row="chr"+row 
          items=row.split()      
          if items[0] in b:#print items[1]
           if items[0] not in a.keys():    
              a[items[0]]=SortedCollection(key=itemgetter(1))
              add=0
          #if ((items[0][0] == "c")):
        #  print items[1]
         # print items[0]
           if items[1] in b[items[0]]:
          # if items[0] == "chrX":
           #    print b[items[0]][items[1]]
            a[items[0]].insert((b[items[0]][items[1]],int(items[1])))
            
            
    file_r=file+"_recurrent.txt.gz"      
    with gzip.open(file_r,'rb') as r:
        for row in r:
          if row[0] == "M":
             row="M"+row[2:len(row)]   
          row="chr"+row 
          items=row.split()      
          #print items[1]
          if items[0] not in a.keys():    
              a[items[0]]=SortedCollection(key=itemgetter(1))
              add=0
          #if ((items[0][0] == "c")):
          if items[1] in b[items[0]]:
           a[items[0]].insert((b[items[0]][items[1]],int(items[1])))
            #print items[2]
            #print items[0]+" "+items[1].__str__()
            #print a[items[0]].find(int(items[1]))
           
           #print len(a["chr"+items[0]])
    return a




def get_regions(cytoband):
    regions={}
    add=0
    for item in cytoband:
        if [item.chrom] not in regions.keys():
            #print item.chrom
            regions[item.chrom]={}
            regions[item.chrom][0]={}
        regions[item.chrom][0][item.chromEnd]="+"
        #regions[item.chrom]["end"]=item.chromEnd
        #regions[item.chrom]["chrom"]=item.chrom
        #regions[item.chrom]["start"]=0
    regions["chrM"]={}
    regions["chrM"][0]={}
    regions["chrM"][0][16569]="+"
    return (regions)

###MAIN PROGRAM###

file='/Users/mok6/Desktop/Doug_NEW/dataNewII/SP04HU_PyPu'

print "start"

#get start and end of chromosomes, this is from dbcruz packages, local copy
regions=get_regions(Genome("sqlite:////Users/mok6/Dropbox/LiClipse/MelArray/hg19_c.db").cytoband)
print "end load chr boundaries"

#load dimer data, second parameter is 0=filtered, 1=not filtered
DamagePos=loadDamagePosII(file,"1")
print "dimer data loaded"

#perform sliding window approach
slidingWindowII(regions,DamagePos, file)
print "sliding window done"

