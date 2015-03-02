import fileinput, collections
from collections import defaultdict
import os, sys
mylist = defaultdict(list)

fileid=1
intputfilename=str(sys.argv[1]);
filename, fileExt = os.path.splitext(intputfilename)
outputfilename=str(sys.argv[3]);

for line in open(sys.argv[1]):
    line=line.rstrip("\n")
    mylist[fileid].append(line)
    fileid=fileid+1

totalfileno=len(mylist)
No_split=int(sys.argv[2]);
filenopersliot=totalfileno/No_split
for i in range(1,No_split+1):
    splitfilename="{0}_{1}_{2}{3}".format(outputfilename,"split",i,fileExt)
    print splitfilename
    f = open(splitfilename,'wt')
    print >> f, mylist[1+(i-1)*filenopersliot][0],
    for j in range(2+(i-1)*filenopersliot,i*filenopersliot+1):
        print >> f, "\n"+mylist[j][0],
    if (i==No_split):
        for j in range(i*filenopersliot+1,totalfileno+1):
            print >> f, "\n"+mylist[j][0],
    f.close()
