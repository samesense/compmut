import sys, cm_scanTester5, pickle

def findLength(filename):
    f = open("../../Doc/pwm/" + filename + ".pwm")
    s = f.readline()
    s = f.readline()
    s = f.readline()
    count = 0
    while s != "":
        s = f.readline()
        count = count + 1
    f.close()
    return count
    
#finds the correct number of data points for
#the random dist
#by looking at each PWM and checking the lengths
def getSampleCount(skip):
    f = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwms = pickle.load(f)
    length = 0
    for p in pwms:
        length += (findLength(p) - skip - 1)
    return length

def main(args):
    skip = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    pvals = [-9]        
    for j in pvals:        
        print "Doing pval " + str(j)        
        for s in skip:
            cm_scanTester5.mainRunner(getSampleCount(s), s)
    sys.exit(0)
        
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
