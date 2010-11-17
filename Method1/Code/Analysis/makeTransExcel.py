import sys, pickle

def findPWMName(pwm):
    f = open("../../Doc/Vertebrate PWM/full_matrix_list.txt", "r")
    s = f.readline()
    while s != "":
        stringList = s.split(";")[0].split()
        if pwm == stringList[0]:
            f.close()
            return stringList[2]
        s = f.readline()    
    f.close()
    return "noName"

def findPWMFam(pwm):    
    f = open("../../Doc/Vertebrate PWM/full_matrix_list.txt", "r")
    s = f.readline()
    while s != "":
        stringList = s.split(";")[0].split()        
        if pwm == stringList[0]:
            f.close()
            classString = ""
            if len(stringList) > 4:
                classString = stringList[3].replace(",","") + "_" + stringList[4]
            else:  classString = stringList[3]            
            return classString
        s = f.readline()    
    f.close()
    return "noFam"

def average(lineSplits):
    acc = float(0)
    for i in range(len(lineSplits)):
        acc += float(lineSplits[i][2])
        
    return str(acc/float(len(lineSplits)))

def accLLs(lineSplits):
    ret = ""
    for i in range(len(lineSplits)):
        ret = ret + lineSplits[i][2] + ","
    return ret

def main(args):
    gaps = [1, 2, 3, 4, 5, 6, 7]
    pwmFile = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwms = pickle.load(pwmFile)
    pwmFile.close()

    for gap in gaps:    
        excelFile = open("../../Results/Transitive/" + str(gap+1) + " Edges.csv", "w")
        for i in range(len(pwms)):
            flag = False
            try:
                singleFile = open("../../Results/Real/Raw Data/Gap 0/" + pwms[i] + "(p " + str(-9) + ").gff.ll", "r")
                edgeFile = open("../../Results/Real/Raw Data/Gap " + str(gap) + "/" + pwms[i] + "(p "+ str(-9)+ ").gff.ll", "r")
                flag = True
            except: pass #print "not here", "../../Results/Real/Raw Data/Gap " + str(gap) + "/" + pwms[i]+ "(p "+ str(9)+ ").gff.ll"
            if flag == True:
                LLFiles = []
                for q in range(gap+1):
                    LLFiles.append(singleFile.readline())                
                edge = edgeFile.readline()
                while LLFiles[gap].split()[0] != "Hit":                
                    lineSplits = []
                    for q in range(gap+1):
                        lineSplits.append(LLFiles[q].split())
                    
                    excelFile.write(findPWMFam(pwms[i]) + "," + findPWMName(pwms[i]) + "," + pwms[i] + "," +\
                                    lineSplits[0][0].replace("(","").replace(",","") + "," + lineSplits[gap][1].replace(")","") + "," + edge.split()[2] + \
                                    "," + accLLs(lineSplits) + average(lineSplits) + "\n")                       
                    for q in range(gap):
                        LLFiles[q] = LLFiles[q+1]
                    LLFiles[gap] = singleFile.readline()
                    edge = edgeFile.readline()
                singleFile.close()
                edgeFile.close()
        excelFile.close()
    sys.exit(0)
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
