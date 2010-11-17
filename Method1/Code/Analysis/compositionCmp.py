import sys, pickle
from rpy import *

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

def addPositions(aPWMDict, pwmName):
    pwmLength = findLength(pwmName)
    for i in range(pwmLength):
        aPWMDict[str(i+1)] = dict()
        aPWMDict[str(i+1)]["info"] = float(0)
        aPWMDict[str(i+1)]["sig"] = False

def parseSig(s):
    ret = []
    ret.append("")
    ret.append("")
    ret.append("")    
    ls = s.split()
    ret[0] = ls[0]
    ret[1] = ls[2][1]
    ret[2] = ls[3][0]
    return ret

def fillInSigAtGap(pwmDict, sigFile, gap):
    sig = open(sigFile + str(gap) + ".sig", "r")
    s = sig.readline()
    while s != "":
        parse = parseSig(s)
        pwmDict[parse[0]][parse[1]]["sig"] = True
        pwmDict[parse[0]][parse[2]]["sig"] = True
        s = sig.readline()
    sig.close()

def fillInSig(pwmDict, sigFile):
    for gap in range(8):
        fillInSigAtGap(pwmDict, sigFile, gap)

def parseInfoLine(pwmDict, s):
    ls = s.split()
    pwmName = ls[0].split(".")[0]    
    for i in range(len(ls)-1-2): #leave out the last two b/c they don't matter
        pwmDict[pwmName][str(i+1)]["info"] = float(ls[i+1])

def fillInInfo(pwmDict):
    infoFile = open("../../Doc/Vertebrate PWM/matrices.vert.IC", "r")    
    s = infoFile.readline()
    while s != "":
        parseInfoLine(pwmDict, s)
        s = infoFile.readline()
    infoFile.close()    

def check(pwmDict):
    pwmKeys = pwmDict.keys()
    pwmKeys.sort()
    for pwmIndex in range(len(pwmKeys)):
        for pos in pwmDict[pwmKeys[pwmIndex]].keys():
            if pwmDict[pwmKeys[pwmIndex]][pos]["sig"] == True:
                print pwmKeys[pwmIndex], pos, pwmDict[pwmKeys[pwmIndex]][pos]["info"]

# ret [domain, range]
def setupPlot(pwmDict):
    ret = []
    ret.append([])
    ret.append([])
    pwmKeys = pwmDict.keys()
    pwmKeys.sort()
    for pwmIndex in range(len(pwmKeys)):
        for pos in pwmDict[pwmKeys[pwmIndex]].keys():
            ret[0].append(pwmDict[pwmKeys[pwmIndex]][pos]["info"])
            if pwmDict[pwmKeys[pwmIndex]][pos]["sig"] == True:
                ret[1].append(1)
            else: ret[1].append(-1)
    return ret

def plotIt(xy):
    r.png("../../Results/Graphs/Information Content/Dot.png", width = 733, height = 900)
    r.par(pty="s")
    r.plot(xy[0], xy[1], xlab="Infomration Content", ylab="Significance")

def plotDist(sig, notSig):
    out = open("sigDist.txt", "w")
    for s in sig:
        out.write(str(s)+"\n")
    out.close()

    out = open("notSigDist.txt", "w")
    for s in notSig:
        out.write(str(s)+"\n")
    out.close()

def findPvalandPlotDist(pwmDict):
    sig = []
    notSig = []
    pwmKeys = pwmDict.keys()
    pwmKeys.sort()
    for pwmIndex in range(len(pwmKeys)):
        for pos in pwmDict[pwmKeys[pwmIndex]].keys():            
            if pwmDict[pwmKeys[pwmIndex]][pos]["sig"] == True:
                sig.append(pwmDict[pwmKeys[pwmIndex]][pos]["info"])
            notSig.append(pwmDict[pwmKeys[pwmIndex]][pos]["info"])
    print "pval", r.wilcox_test(sig, notSig, alternative="greater", paired=False)["p.value"]
    plotDist(sig, notSig)

def main(args):
    #get all possible positions
    pwmFile = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwms = pickle.load(pwmFile)
    pwmFile.close()
    pwmDict = dict()
    for aPWM in pwms:
        pwmDict[aPWM] = dict()
        addPositions(pwmDict[aPWM], aPWM)
    fillInInfo(pwmDict)
    sigFile = "../../Results/Real/Sig Vals/Threshold/Point 5 Threshold/Gap "
    #sigFile = "../../Results/Real/Sig Vals/One Percent Same Site Random/Gap "
    #sigFile = "../../Results/Real/Sig Vals/One Percent Real/Gap "
    fillInSig(pwmDict, sigFile)
    #check(pwmDict)     
    plotIt(setupPlot(pwmDict))
    findPvalandPlotDist(pwmDict)
    sys.exit(0)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
