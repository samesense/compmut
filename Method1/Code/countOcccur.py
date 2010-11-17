import sys, pickle

def main(args):
    pwms = []
    dpwms = dict()
    for i in range(111):
        if i+1 < 10:                
            pwms.append("MA000" + str(i+1) + "(p " + str(-9) + ").gff")
        elif i+1 < 100:
            pwms.append("MA00" + str(i+1) + "(p " + str(-9) + ").gff")
        else:
            pwms.append("MA0" + str(i+1) + "(p " + str(-9) + ").gff")
    #for i in range(11):
    #    try:
    p = open("random_logs_perryTest.MA0005(p -9).gff.txt", "r")
            #p = open("random_logs_perryTest." + pwms[i+1] + ".txt", "r")
            #print pwms[i]
    f2 = pickle.load(p)  
    
    for s2 in f2:
        s2 = round(float(s2)*100)/100
        if s2 > 1.5: print s2
        if dpwms.has_key(s2):
            dpwms[s2] = dpwms[s2] + 1
        else: dpwms[s2] = 1             
    p.close()
        #except: print "NO HITS"
        
    f = open("todump.txt", "w")
    pickle.dump(dpwms, f)
    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
