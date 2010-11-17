import sys, cm_scan3, pickle, cm_scan3NA

def findLength(filename):
    f = open("../doc/fasta/" + filename)
    s = f.readline()
    count = 1
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

    skip = [1, 2, 3, 4, 5, 6]
    pvals = [-9, -10]
    len1 = 0
    len2 = 0
    f1 = ""
    f2 = ""
    offset = 0
    for j in pvals:
        #if j == -9: offset = 60
        #else: offset = 0
        print "Doing pval " + str(j)
        pwms = []
        pwmsShuffled = []
        for i in range(111):
            if i+1 < 10:                
                pwms.append("MA000" + str(i+1) + "(p " + str(j) + ").gff")
                pwmsShuffled.append("MA000" + str(i+1) + "random(p " + str(j) + ").gff")
            elif i+1 < 100:
                pwms.append("MA00" + str(i+1) + "(p " + str(j) + ").gff")
                pwmsShuffled.append("MA00" + str(i+1) + "random(p " + str(j) + ").gff")
            else:
                pwms.append("MA0" + str(i+1) + "(p " + str(j) + ").gff")
                pwmsShuffled.append("MA0" + str(i+1) + "random(p " + str(j) + ").gff")            
        for i in range(len(pwms)-offset):            
            f1 = open("../doc/fasta/" + pwms[i+offset], "r")
            f2 = open("../doc/fasta/" + pwmsShuffled[i+offset], "r")
            m1 = cm_scan3.buildMotifs(f1, alignment)
            m2 = cm_scan3.buildMotifs(f2, alignment)
            len1 = len(m1)
            len2 = len(m2)
            f1.close()
            f2.close()
            try:
                print len1, len2
                if len1 != 0 and len2 != 0:
                    print pwms[i+offset]
                    cm_scan3NA.mainRunner(pwms[i+offset], m1, min([len1, len2]), "LL all PWM pval " + str(j) + "/")
                    print pwmsShuffled[i+offset]
                    cm_scan3NA.mainRunner(pwmsShuffled[i+offset], m2, min([len1, len2]), "LL all PWMRandom pval " + str(j) + "/")
                else: print "NO HITS"
            except KeyError:
                print "\t\tkey problem", pwms[i+offset]               
    sys.exit(0)
        
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
