import sys, pickle

def findPWMLength(filename):
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

def main(args):
##    for skip in range(9):
##        f = open("../../Results/Real/Sig Vals/Gap " + str(skip) + ".sig")
##        s = f.readline()
##        count = 1
##        while s != "":
##            count += 1
##            s = f.readline()
##        print "Gap", skip, count
    f = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwmList = pickle.load(f)
    f.close()
    pwmInfo = dict()
    for p in pwmList:
        pwmInfo[p] = dict()
        pwmInfo[p]["length"] = findPWMLength(p)
    count = 0
    for skip in range(8):        
        for p in pwmList:
            if pwmInfo[p]["length"] > skip-2:
                try:                    
                    pwmHitFile = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/" + p + "(p -9).gff.ll", "r")
                    pwmHitFile.close()
                    count += pwmInfo[p]["length"]-1-skip
                except: pass
    print count

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
