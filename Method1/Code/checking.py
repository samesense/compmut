import sys, pickle

def main(args):
    count = 0
    pwms = []
    for i in range(111):
        if i+1 < 10:                
            pwms.append("LL all PWM pval -9 NA 1/MA000" + str(i+1) + "(p " + str(-9) + ").gff.ll")
        elif i+1 < 100:
            pwms.append("LL all PWM pval -9 NA 1/MA00" + str(i+1) + "(p " + str(-9) + ").gff.ll")
        else:
            pwms.append("LL all PWM pval -9 NA 1/MA0" + str(i+1) + "(p " + str(-9) + ").gff.ll")
    for f in pwms:
        try:
            q = open(f, "r")
            s = q.readline()
            while s != "":
                count += 1
                s = q.readline()
            count -= 1
        except: pass
    print count
    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
