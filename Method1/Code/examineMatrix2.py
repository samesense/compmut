import sys, pickle
#called with pwm skip pval
#pwm is just a number
def main(args):
    #args pwm skip pval
    pwm = args[1]
    skip = args[2]
    pval = args[3]
    pre = ""
    if pwm < 10:
        pre = "MA000"
    elif pwm < 100:
        pre = "MA00"
    else: pre = "MA0"
    realMatrixFile = open("LL all PWM pval " + str(pval) +\
                          " NA " + str(skip) + "/" + pre + str(pwm) +\
                          "(p " + str(pval) + ").gff.256sq", "r")
    realMatrix = pickle.load(realMatrixFile)
    realMatrixFile.close()
    
    nucs = ["A", "C", "G", "T"]
    edges = []
    for n1 in nucs:
        for n2 in nucs:
            edges.append(n1+n2)
    print "Real Matrix"
    i = 1
    for matrix in realMatrix:
        print "pos", i, "x pos", i+1+int(skip)
        print "\tAA \tAC \tAG \tAT \tCA \tCC \tCG \tCT \tGA \tGC \tGG \tGT \tTA \tTC \tTG \tTT"        
        for n1 in edges:            
            print n1, "\t", matrix[n1+"AA"], "\t", matrix[n1+"AC"], \
                  "\t", matrix[n1+"AG"], "\t", matrix[n1+"AT"], \
                  "\t", matrix[n1+"CA"], "\t", matrix[n1+"CC"], \
                  "\t", matrix[n1+"CG"], "\t", matrix[n1+"CT"], \
                  "\t", matrix[n1+"GA"], "\t", matrix[n1+"GC"], \
                  "\t", matrix[n1+"GG"], "\t", matrix[n1+"GT"], \
                  "\t", matrix[n1+"TA"], "\t", matrix[n1+"TC"], \
                  "\t", matrix[n1+"TG"], "\t", matrix[n1+"TT"]
        i += 1
        
        sys.exit(0)
    
    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
