import sys, pickle
def main(args):
    pwmFile = open("PWMlist.txt", "r")
    pwms = pickle.load(pwmFile)
    pwmFile.close() 
    outPut = open("master_results.txt", "w")
    for pwm in pwms:         
        matFile = open("matrix_list.txt", "r")
        l = matFile.readline()
        while l != "":
            if l.split()[0] == pwm:                    
                if l.split()[4] != ";": outPut.write(pwm + " " +l.split()[3] + " " + l.split()[4] +"\n")
                else: outPut.write(pwm + " " + l.split()[3] + "\n")
            l = matFile.readline()
        matFile.close()
        flag = False        
        for s in range(9):
            if s != 7:
                try:
                    hitFile = open("LL all PWM pval -9 NA " + str(s) + "/" + pwm + "(p -9).gff.ll", "r")
                    l = hitFile.readline()
                    while l.split()[0] != "Hit": l = hitFile.readline()
                    hitFile.close()
                    hits = l.split()[2]
                except: pass
                f = open("Results NA " + str(s) + " p -9.txt", "r")
                l = f.readline()            
                while l.split()[0] != "pval":
                    if l.split()[0] == pwm:
                        togo = l.split()
                        outPut.write(togo[0] + "  " + togo[1] + "  " + togo[2]+togo[3] + "  " + \
                                     str(round(float(togo[4])*100)/100) + "  " + str(hits) +"\n")
                        flag = True
                    l = f.readline()
                f.close()
        if flag: outPut.write("\n")
            
    outPut.close()
    sys.exit(0)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
