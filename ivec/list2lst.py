import fileinput,collections
import sys

listname=str(sys.argv[1]);
lstname=str(sys.argv[2]);
fp_lst=open(lstname,"wt")
for theline1 in file (listname):
    theline1=theline1.rstrip()
    keyline=theline1.split("/")
    label=-1
    if ("alv" in keyline[-1]):
        label=1
    elif ("fas" in keyline[-1]):
        label=2
    elif ("prs" in keyline[-1]):
        label=3
    elif ("pus" in keyline[-1]):
        label=4
    elif ("urd" in keyline[-1]):
        label=5
    elif ("xxx" in keyline[-1]):
        label=6
    else:
        keyline2=keyline[-1].split("_")
        if(keyline2[0]=="fa"):
        	label=2
        else:
        	print "error  %s" % (theline1)
                   
    print >> fp_lst, "{0}".format(label)        
fp_lst.close()

