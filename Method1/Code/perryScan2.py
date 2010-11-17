import sys, cm_scan3, pickle, cm_scan3NA, makeGraphs, cm_scanTester2

def findLength(filename):
    f = open("../doc/pwm/" + filename)
    s = f.readline()
    s = f.readline()
    s = f.readline()
    count = 0
    while s != "":
        s = f.readline()
        count = count + 1
    f.close()
    return count

def main(args):
    a_file = open("../doc/pickle/alignment.pickle", "r")

    # Load the 5-way alignment
    alignment = pickle.load(a_file)

    # Close relevant files
    a_file.close()

    skip = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    pvals = [-9, -10]
    len1 = 0    
    f1 = ""    
    for j in pvals:        
        print "Doing pval " + str(j)
        pwms = []
        
        for i in range(111):
            if i+1 < 10:                
                pwms.append("MA000" + str(i+1))                
            elif i+1 < 100:
                pwms.append("MA00" + str(i+1))                
            else:
                pwms.append("MA0" + str(i+1))
        for s in skip:                    
            for i in range(len(pwms)):                
                f1 = open("../doc/fasta/" + pwms[i] + "(p " + str(j) + ").gff", "r")            
                m1 = cm_scan3.buildMotifs(f1, alignment)            
                len1 = findLength(pwms[i] + ".pwm")            
                f1.close()
                print len(m1)
                if len1 > s+1:
                    try:                
                        if len(m1) != 0:
                            print pwms[i]
                            cm_scan3NA.mainRunner(pwms[i] + "(p " + str(j) + ").gff", m1, len(m1), "LL all PWM pval " + str(j) + " NA " + str(s) + "/", s)                    
                        else: print "NO HITS"
                    except KeyError:
                        print "\t\tkey problem", pwms[i]
            makeGraphs.make("LL all PWM pval " + str(j) + " NA " + str(s) + "/", j, s)
            cm_scanTester2.mainRunner(s, j)
    sys.exit(0)
        
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
