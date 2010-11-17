import sys
import pickle
from rpy import *

def main(args):
    print "Testing..."
    #the order of the openings of the files matters, because the first one
    #sets the axis range
    #theFile = open("cdist.txt", "r")
    theFile = open("Normal log pval -9 NA 1.txt", "r")
    p = pickle.load(theFile)    
    x = []
    y = []
    print len(p)
    #print p
    keys = p.keys()
    for i in range(len(keys)):
        x.append(keys[i])
        y.append(p[keys[i]])   
    r.png("compare.png", width=733, height=900)
    r.par(pty="s")
    r.do_points=1
    r.plot(x, y, col="green", xlab="Log Likelihood", ylab="Count")

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open("random_logs 1.txt.results", "r")
    p = pickle.load(theFile)
    #print p
    x = []
    y = []
    print len(p)
    keys = p.keys()
    for i in range(len(keys)):
        x.append(keys[i])
        y.append(p[keys[i]])        
    r.points(x, y, col="black")

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open("REALGOOD_random_logs 1.txt.results", "r")
    p = pickle.load(theFile)
    #print p
    x = []
    y = []
    print len(p)
    keys = p.keys()
    for i in range(len(keys)):
        x.append(keys[i])
        y.append(p[keys[i]])        
    r.points(x, y, col="red") 

##    theFile.close()
##    theFile = open("cdist.txt", "r")
##    #theFile = open("Shuffled log pval -8.txt", "r")
##    p = pickle.load(theFile)
##    #print p
##    x = []
##    y = []
##    print len(p)
##    keys = p.keys()
##    for i in range(len(keys)):
##        x.append(keys[i])
##        y.append(p[keys[i]])        
##    r.points(x, y, col="green")
##
##    r.dev_off()
##    theFile.close()
##    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF


    
