#this file goes through all the raw data and bins it
#it also calculates pvalues
import sys, pickle, newPlot
from rpy import *

def average(alist):
    acc = float(0)
    for i in alist: acc += float(i)
    return acc/float(len(alist))

def stdDev(alist, avg):
    acc = float(0)
    for i in alist: acc += (float(i) - avg)*(float(i) - avg)
    return math.sqrt(acc/float(len(alist)))    

def binData(f1prefix, pval, skip):
    dpwmsReal = dict()
    dpwmsRand = dict()
    dpwmsPureRandom = dict()
    dpwmsRandomPWM = dict()
    dpwmsChangeCompPWM = dict()
    distReal = []
    distRandom = []
    distRandomPWM = []
    distChangeCompPWM = []
    distPureRandom = []
    
    #do the random    
    p = open("../../Results/Controls/Same Site/Raw Data/same site Gap " + str(skip) + ".txt", "r")    
    f2 = pickle.load(p)  
    
    for s2 in f2:
        s2 = round(float(s2)*100)/100
        distRandom.append(s2)        
        if dpwmsRand.has_key(str(s2)):            
            dpwmsRand[str(s2)] = dpwmsRand[str(s2)] + 1
        else: dpwmsRand[str(s2)] = 1             
    p.close()
    f = open("../../Results/Controls/Same Site/Binned Data/same site Gap " + str(skip) + ".binned", "w")
    pickle.dump(dpwmsRand, f)
    f.close()
    dump = open("../../Results/Controls/Same Site/Binned Data/same site Gap " + str(skip) + ".txt", "w")
    for k in dpwmsRand.keys():
        dump.write(k + "," + str(dpwmsRand[k])+"\n")
    dump.close()
    
    #now distReal and distRandom have information for the stats
    distRandom.sort()
    for anIndex in range(len(distRandom)-1):
        if distRandom[anIndex] > distRandom[anIndex+1]: print "no"
    length = len(distRandom)
    top = int(math.floor(length*.01))
    onePer = math.floor(float(.99)*float(length))    
    cutOffVal = distRandom[int(onePer)]
    print skip, cutOffVal, distRandom[len(distRandom)-1-top], length, top
    #get the mean and stdDev
    distRandomMean = average(distRandom)
    distRandomStdDev = stdDev(distRandom, distRandomMean)
    
    outFile = open("../../Results/Real/Sig Vals/One Percent Same Site Random/Gap " + str(skip) + ".sig", "w")
    outFile2 = open("../../Results/Real/Sig Vals/Threshold/Point 5 Threshold/Gap " + str(skip) + ".sig", "w")
    outFile3 = open("../../Results/Real/Sig Vals/One Percent Real/Gap " + str(skip) + ".sig", "w")
    outFile4 = open("../../Results/Real/Sig Vals/Five Percent Real/Gap " + str(skip) + ".sig", "w")
    allRealVals = []
    allRealVals2 = []
    myVals = []
    #do the realDistribution
    pwmFile = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwms = pickle.load(pwmFile)
    pwmFile.close()                

    for i in range(len(pwms)):
        flag = False
        try:
            f1 = open(f1prefix+pwms[i]+"(p "+str(pval)+").gff.ll", "r")
            flag = True
        except: pass # print "not here", pwms[i]
        if flag == True:
            s1 = f1.readline()                

            while s1.split()[0] != "Hit":
                v1 = round(float(s1.split()[2])*100)/100
                allRealVals.append(v1)
                allRealVals2.append(v1)
                diff = v1 - distRandomMean
                zscore = diff/distRandomStdDev
                lineSplit = s1.split()
                myVals.append(pwms[i] + "\t" + str(skip) + "\t" + lineSplit[0] + " " + lineSplit[1] + \
                                  " " + lineSplit[2] + "\t" + str(zscore) + "\n")
                if v1 > cutOffVal:
                    #get zscore for v1
                    diff = v1 - distRandomMean
                    zscore = diff/distRandomStdDev
                    lineSplit = s1.split()
                    outFile.write(pwms[i] + "\t" + str(skip) + "\t" + lineSplit[0] + " " + lineSplit[1] + \
                                  " " + lineSplit[2] + "\t" + str(zscore) + "\n")
                if v1 > .5:
                    #get zscore for v1
                    diff = v1 - distRandomMean
                    zscore = diff/distRandomStdDev
                    lineSplit = s1.split()
                    outFile2.write(pwms[i] + "\t" + str(skip) + "\t" + lineSplit[0] + " " + lineSplit[1] + \
                                  " " + lineSplit[2] + "\t" + str(zscore) + "\n")
                distReal.append(v1)                    
                if dpwmsReal.has_key(str(v1)):
                    dpwmsReal[str(v1)] = dpwmsReal[str(v1)] + 1
                else: dpwmsReal[str(v1)] = 1                   
                s1 = f1.readline()                   
            f1.close()
        
    #outFile.write("pval " +  str())
    out = open("../../Results/Real/Binned Data/Gap " + str(skip) + ".binned", "w")
    pickle.dump(dpwmsReal, out)
    out.close()    
    outFile.close()
    outFile2.close()

    allRealVals2.sort()
    length = len(allRealVals2)
    top = int(math.floor(length*.01))
    onePer = math.floor(float(.99)*float(length))    
    cutOffVal = allRealVals2[int(onePer)]
    
    for i in range(len(allRealVals)):
        if allRealVals[i] > cutOffVal:
            #print allRealVals[i], cutOffVal
            outFile3.write(myVals[i])
    outFile3.close()

    top = int(math.floor(length*.05))
    fivePer = math.floor(float(.95)*float(length))    
    cutOffVal = allRealVals2[int(fivePer)]
    
    for i in range(len(allRealVals)):
        if allRealVals[i] > cutOffVal:
            #print allRealVals[i], cutOffVal
            outFile4.write(myVals[i])
    outFile4.close()
    
    #now do the randmoPWMDistribution    
    for i in range(len(pwms)):
        flag = False
        try:
            f1 = open("../../Results/Controls/Constructed PWMs/Raw Data/Gap " + str(skip) + "/sortedhuman1ku-corrected.fasta." + \
                      pwms[i]+"random.gff.ll", "r")
            flag = True
        except: pass # print "not here", pwms[i]
        if flag == True:
            s1 = f1.readline()                

            while s1.split()[0] != "Hit":
                v1 = round(float(s1.split()[2])*100)/100
                distRandomPWM.append(v1)                    
                if dpwmsRandomPWM.has_key(str(v1)):
                    dpwmsRandomPWM[str(v1)] = dpwmsRandomPWM[str(v1)] + 1
                else: dpwmsRandomPWM[str(v1)] = 1                   
                s1 = f1.readline()                   
            f1.close()
        
    #outFile.write("pval " +  str(float(r.wilcox_test(distReal, distRandom, alternative="greater", paired=False)["p.value"])))
    out = open("../../Results/Controls/Constructed PWMs/Binned Data/constructed pwms Gap " + str(skip) + ".binned", "w")
    pickle.dump(dpwmsRandomPWM, out)
    out.close()
    
    for i in range(len(pwms)):
        flag = False
        try:
            f1 = open("../../Results/Controls/Comp Permuted PWMs/Raw Data/Gap " + str(skip) + "/sortedhuman1ku-corrected.fasta." + pwms[i]+"cc.gff.ll", "r")
            flag = True
        except: pass # print "not here", pwms[i]
        if flag == True:
            s1 = f1.readline()                

            while s1.split()[0] != "Hit":
                v1 = round(float(s1.split()[2])*100)/100               
                distChangeCompPWM.append(v1)                    
                if dpwmsChangeCompPWM.has_key(str(v1)):
                    dpwmsChangeCompPWM[str(v1)] = dpwmsChangeCompPWM[str(v1)] + 1
                else: dpwmsChangeCompPWM[str(v1)] = 1                   
                s1 = f1.readline()                   
            f1.close()
        
    #outFile.write("pval " +  str(float(r.wilcox_test(distReal, distRandom, alternative="greater", paired=False)["p.value"])))
    out = open("../../Results/Controls/Comp Permuted PWMs/Binned Data/permuted pwms Gap " + str(skip) + ".binned", "w")
    pickle.dump(dpwmsChangeCompPWM, out)
    out.close()

    #do the completely random case
    p = open("../../Results/Controls/Pure Random/Raw Data/pure random Gap " + str(skip) + ".txt", "r")            
    f2 = pickle.load(p)  
    
    for s2 in f2:
        s2 = round(float(s2)*100)/100
        distPureRandom.append(s2)
        if dpwmsPureRandom.has_key(str(s2)):
            dpwmsPureRandom[str(s2)] = dpwmsPureRandom[str(s2)] + 1
        else: dpwmsPureRandom[str(s2)] = 1             
    p.close()
    f = open("../../Results/Controls/Pure Random/Binned Data/pure random Gap " + str(skip) + ".binned", "w")
    pickle.dump(dpwmsPureRandom, f)
    f.close()

    #now find the pvals for the realDist vs. the other dists
    pvalDict = dict()
    pvalDict["same site"] = float(r.wilcox_test(distReal, distRandom, alternative="greater", paired=False)["p.value"])
    pvalDict["pure random"] = float(r.wilcox_test(distReal, distPureRandom, alternative="greater", paired=False)["p.value"])
    pvalDict["change comp"] = float(r.wilcox_test(distReal, distChangeCompPWM, alternative="greater", paired=False)["p.value"])
    pvalDict["constructed"] = float(r.wilcox_test(distReal, distRandomPWM, alternative="greater", paired=False)["p.value"])

    f = open("../../Results/Controls/Pvals/pvals Gap " + str(skip) + ".txt", "w")
    pickle.dump(pvalDict, f)
    f.close()
    return pvalDict

def main(args):
    skip = range(9)    
    pvalList = [-9]    
    for j in pvalList:       
        for s in skip:
            pvalDict = binData("../../Results/Real/Raw Data/Gap " + str(s) + "/", j, s)            
            pvals = []
            pvals.append(pvalDict["same site"])
            pvals.append(pvalDict["constructed"])
            pvals.append(pvalDict["change comp"])
            pvals.append(pvalDict["pure random"])
                        
            newPlot.binnedPlot("../../Results/Real/Binned Data/Gap " + str(s) + ".binned", \
                               ["../../Results/Controls/Same Site/Binned Data/same site Gap " + str(s) + ".binned",\
                                "../../Results/Controls/Constructed PWMs/Binned Data/constructed pwms Gap " + str(s) + ".binned",\
                                "../../Results/Controls/Comp Permuted PWMs/Binned Data/permuted pwms Gap " + str(s) + ".binned",\
                                "../../Results/Controls/Pure Random/Binned Data/pure random Gap " + str(s) + ".binned"],
                               pvals,\
                               ["green", "red", "black", "orange", "blue"],\
                               "../../Results/Graphs/Line/Gap " + str(s) + ".png")
            newPlot.dottedPlot("../../Results/Real/Binned Data/Gap " + str(s) + ".binned", \
                               ["../../Results/Controls/Same Site/Binned Data/same site Gap " + str(s) + ".binned",\
                                "../../Results/Controls/Constructed PWMs/Binned Data/constructed pwms Gap " + str(s) + ".binned",\
                                "../../Results/Controls/Comp Permuted PWMs/Binned Data/permuted pwms Gap " + str(s) + ".binned",\
                                "../../Results/Controls/Pure Random/Binned Data/pure random Gap " + str(s) + ".binned"],
                               pvals,\
                               ["green", "red", "black", "orange", "blue"],\
                               "../../Results/Graphs/Points/Gap " + str(s) + ".png")
    sys.exit(0)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
