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
                          "(p " + str(pval) + ").gff.16sq", "r")
    realMatrix = pickle.load(realMatrixFile)
    realMatrixFile.close()
    randomMatrixFile = open("NA " + str(skip) + " p " + str(pval) +\
                            " " + str(pwm) + ".16sq", "r")
    randomMatrix = pickle.load(randomMatrixFile)
    realMatrixFile.close()
    
    nucs = ["A", "C", "G", "T"]
    print "Real Matrix"
    i = 1
    for matrix in realMatrix["human"]:
        print "pos", i, "x pos", i+1+int(skip)
        print "\tA \tC \tG \tT"        
        for n1 in nucs:            
            print n1, "\t", matrix[n1+"A"], "\t", matrix[n1+"C"], \
                  "\t", matrix[n1+"G"], "\t", matrix[n1+"T"]
        i += 1

    print "Random Matrix"
    i = 1
    for matrix in randomMatrix["human"]:
        print "pos", i, "x pos", i+1+int(skip)
        print "\tA \tC \tG \tT"        
        for n1 in nucs:            
            print n1, "\t", matrix[n1+"A"], "\t", matrix[n1+"C"], \
                  "\t", matrix[n1+"G"], "\t", matrix[n1+"T"]
        i += 1
    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
