import sys, pickle

def partition(array, begin, end, cmp):
    while begin < end:
         while begin < end:
            if cmp(array[begin][1], array[end][1]):
                (array[begin], array[end]) = (array[end], array[begin])
                break
            end -= 1
         while begin < end:
            if cmp(array[begin][1], array[end][1]):
                (array[begin], array[end]) = (array[end], array[begin])
                break
            begin += 1
    return begin

def sort(array, cmp=lambda x, y: x > y, begin=None, end=None):
    if begin is None: begin = 0
    if end   is None: end   = len(array)
    if begin < end:
        i = partition(array, begin, end-1, cmp)
        sort(array, cmp, begin, i)
        sort(array, cmp, i+1, end)

def PWMSigVals(pwm, skipLimit):    
    acc = 0
    for skip in pwm.keys():
        if skip != "name":
            if int(skip) <= skipLimit:
                acc += len(pwm[skip]["sigVals"])
    return acc

def PWMSigValsNormalized(pwm, skipLimit):    
    acc = 0
    possible = 0
    for skip in pwm.keys():
        if skip != "name":
            if int(skip) <= skipLimit and int(skip) < pwm[str(0)]["length"]:
                totalvals = len(pwm[skip]["sigVals"])
                #if totalvals != 0:
                acc += totalvals
                #try:
                possible += pwm[skip]["length"]
                #except:
                #    print pwm["name"], skip
    return [acc, possible]

#return a list of lists
#each list inside [pwm, hits]
def sortPWMsByHits(families, skipLimit):
    pwmSigHits = []
    for fkey in families.keys():
        for pwmKey in families[fkey]:
            pwmSigHits.append([pwmKey, PWMSigVals(families[fkey][pwmKey], skipLimit)])
    sort(pwmSigHits)
    return pwmSigHits

#return a list of lists
#each list inside [pwm, hits]
def sortPWMsByHitsNormalized(families, skipLimit):
    pwmSigHits = []
    for fkey in families.keys():
        for pwmKey in families[fkey]:
            hitObject = PWMSigValsNormalized(families[fkey][pwmKey], skipLimit)
            pwmSigHits.append([pwmKey, float(hitObject[0])/float(hitObject[1])])
    sort(pwmSigHits)
    return pwmSigHits

def sortPWMsByHitsAtSkip(families, skip):
    pwmSigHits = []
    for fkey in families.keys():
        for pwmKey in families[fkey]:
            pwmSigHits.append([pwmKey, len(families[fkey][pwmKey][str(skip)]["sigVals"])])
    sort(pwmSigHits)
    return pwmSigHits

def sortPWMsByHitsAtSkipNormalized(families, skip):
    pwmSigHits = []
    for fkey in families.keys():
        for pwmKey in families[fkey]:
            thelen = len(families[fkey][pwmKey][str(skip)]["sigVals"])
            if thelen == 0: pwmSigHits.append([pwmKey, 0])
            else: pwmSigHits.append([pwmKey, float(thelen)/float(families[fkey][pwmKey][str(skip)]["length"])])
    sort(pwmSigHits)
    return pwmSigHits 

def sortFamiliesByHits(families, skipLimit):
    familyList = []
    for fkey in families.keys():
        acc = 0
        for pwmKey in families[fkey]:
            acc += PWMSigVals(families[fkey][pwmKey], skipLimit)
        familyList.append([fkey, acc])
    sort(familyList)
    return familyList

def sortFamiliesByHitsNormalized(families, skipLimit):
    familyList = []
    for fkey in families.keys():
        acc = 0
        possible = 0
        for pwmKey in families[fkey]:
            hitObject = PWMSigValsNormalized(families[fkey][pwmKey], skipLimit)
            acc += hitObject[0]
            possible += hitObject[1]
        if acc == 0: familyList.append([fkey, float(0)])
        else:
            familyList.append([fkey, float(acc)/float(possible)])
    sort(familyList)
    return familyList

def sortFamiliesByHitsAtSkip(families, skip):
    familyList = []
    for fkey in families.keys():
        acc = 0
        for pwmKey in families[fkey]:
            acc += len(families[fkey][pwmKey][str(skip)]["sigVals"])
        familyList.append([fkey, acc])
    sort(familyList)
    return familyList

def sortFamiliesByHitsAtSkipNormalized(families, skip):
    familyList = []
    for fkey in families.keys():
        acc = 0
        possible = 0
        for pwmKey in families[fkey]:
            totalvals = len(families[fkey][pwmKey][str(skip)]["sigVals"])
            if totalvals != 0:
                acc += totalvals
                possible += families[fkey][pwmKey][str(skip)]["length"]
        if acc == 0: familyList.append([fkey, float(0)])
        else:
            familyList.append([fkey, float(acc)/float(possible)])
    sort(familyList)
    return familyList

def distanceDistributionsByPWM(families):
    skipList = []
    for skip in range(13):
        skipList.append(sortPWMsByHitsAtSkip(families, skip))
    pwmDict = dict()
    for fkey in families.keys():
        for pwmKey in families[fkey]:
            pwmDict[pwmKey] = dict()
            for skip in range(13):
                flag = False
                for pwm in skipList[skip]:
                    if pwm[0] == pwmKey: pwmDict[pwmKey][str(skip)] = pwm[1]
                    flag = True
                if flag == False: pwmDict[pwmKey][str(skip)] = 0
    return pwmDict

def distanceDistributionsByPWMNormalized(families):
    skipList = []
    for skip in range(13):
        skipList.append(sortPWMsByHitsAtSkipNormalized(families, skip))
    pwmDict = dict()
    for fkey in families.keys():
        for pwmKey in families[fkey]:
            pwmDict[pwmKey] = dict()
            for skip in range(13):
                flag = False
                for pwm in skipList[skip]:
                    if pwm[0] == pwmKey: pwmDict[pwmKey][str(skip)] = pwm[1]
                    flag = True
                if flag == False: pwmDict[pwmKey][str(skip)] = 0
    return pwmDict    

#need to find the distance distribution per family
#so many adjacent, so many one away, and so on
def distanceDistributionsByFamily(families):
    skipList = []
    for skip in range(13):
        skipList.append(sortFamiliesByHitsAtSkip(families, skip))
    familyDict = dict()
    for fkey in families.keys():
        familyDict[fkey] = dict()
        for skip in range(13):
            flag = False
            for fam in skipList[skip]:
                if fam[0] == fkey: familyDict[fkey][str(skip)] = fam[1]
                flag = True
            if flag == False: familyDict[fkey][str(skip)] = 0
    return familyDict

#need to find the distance distribution per family
#so many adjacent, so many one away, and so on
def distanceDistributionsByFamilyNormalized(families):
    skipList = []
    for skip in range(13):
        skipList.append(sortFamiliesByHitsAtSkipNormalized(families, skip))
    familyDict = dict()
    for fkey in families.keys():
        familyDict[fkey] = dict()
        for skip in range(13):
            flag = False
            for fam in skipList[skip]:
                if fam[0] == fkey: familyDict[fkey][str(skip)] = fam[1]
                flag = True
            if flag == False: familyDict[fkey][str(skip)] = 0
    return familyDict

#famRet has familyName ->skip->dist of deps at this skip
def distanceDistByFamilyPWMStyle(families):
    famRet = dict()
    #This finds how many deps are at each site
    famDict = distanceDistributionsByFamily(families)
    pwmDict = distanceDistributionsByPWM(families)
    for famTitle in families.keys():
        famRet[famTitle] = dict()
        for s in range(len(famDict[famTitle])): #skip positions
            pwmTotal = 0
            famRet[famTitle][str(s)] = dict() #create a new dist at this skip
            for pwm in families[famTitle].keys():
                pwmTotal += 1                        
                if famRet[famTitle][str(s)].has_key(str(pwmDict[pwm][str(s)])):
                    famRet[famTitle][str(s)][str(pwmDict[pwm][str(s)])] += 1
                else: famRet[famTitle][str(s)][str(pwmDict[pwm][str(s)])] = 1
            for depKey in famRet[famTitle][str(s)].keys():
                famRet[famTitle][str(s)][depKey] /= float(pwmTotal)
    return famRet

#famRet has familyName ->skip->dist of deps at this skip
def distanceDistByFamilyPWMStyleNormalized(families):
    famRet = dict()
    #This finds how many deps are at each site
    famDict = distanceDistributionsByFamilyNormalized(families)
    pwmDict = distanceDistributionsByPWMNormalized(families)
    for famTitle in families.keys():
        famRet[famTitle] = dict()
        for s in range(len(famDict[famTitle])): #skip positions
            pwmTotal = 0
            famRet[famTitle][str(s)] = dict() #create a new dist at this skip
            for pwm in families[famTitle].keys():
                pwmTotal += 1
                val1 = pwmDict[pwm][str(s)]                
                val = round(val1*float(100))/float(100)
                
                if famRet[famTitle][str(s)].has_key(str(val)):
                    famRet[famTitle][str(s)][str(val)] += 1
                else: famRet[famTitle][str(s)][str(val)] = 1
            for depKey in famRet[famTitle][str(s)].keys():
                famRet[famTitle][str(s)][depKey] /= float(pwmTotal)
    return famRet

def countFamilyCmuts(families, fkey):    
    possible = 0
    for pwmKey in families[fkey]:
        for skip in range(8):
            if families[fkey][pwmKey][str(skip)].has_key("length"):
                possible += families[fkey][pwmKey][str(skip)]["length"]
    return possible

def countPWMCmuts(families, fkey, pwmKey):
    possible = 0    
    for skip in range(8):
        if families[fkey][pwmKey][str(skip)].has_key("length"):
            possible += families[fkey][pwmKey][str(skip)]["length"]
    print possible
    return possible

def makeFamilyExcel(families, sig):
    linesTogo = []
    sorts = []
    sorts2 = []
    startSort = sortFamiliesByHitsNormalized(families, 7)
    #startSort = sortFamiliesByHits(families, 7)
    for aSkip in range(8): sorts.append(sortFamiliesByHitsAtSkipNormalized(families, aSkip))
    for aSkip in range(8): sorts2.append(sortFamiliesByHitsAtSkip(families, aSkip))
    for hit in startSort:
        family = hit[0]
        skipList = []
        skipList2 = []
        for aSkip in range(8):
            flag = False
            for aHit in sorts2[aSkip]:
                #print "name", aHit[0]
                if family == aHit[0]:
                    #print "called"
                    flag = True
                    skipList.append(aHit[1])
            for aHit in sorts[aSkip]:
                if family == aHit[0]:
                    skipList2.append(aHit[1])
            if flag == False:
                skipList.append(0)
                print "not found", family, sorts2[aHit][0]
        commaSplit = family.split(",")
        familynewname = family
        if len(commaSplit) > 1:
            newstr = ""
            for strindex in range(len(commaSplit[1])-1):
                newstr = newstr + commaSplit[1][strindex+1]
            familynewname = commaSplit[0] + ":" + newstr
        #print family
        alineTogo = familynewname + "," + str(len(families[family].keys())) + "," + str(countFamilyCmuts(families, family)) + "," + str(hit[1])
        for aSkipScore in range(len(skipList)):
            alineTogo = alineTogo + "," + str(skipList[aSkipScore]) + "," + str(skipList2[aSkipScore])
        linesTogo.append(alineTogo + "\n")
        
    if sig == str(1):    
        f = open("../../Results/Summaries/Top One Percent Same Site Random/tempFamOut.csv", "w")
    elif sig == str(2):    
        f = open("../../Results/Summaries/One Percent Real/tempFamOut.csv", "w")
    else:
        f = open("../../Results/Summaries/Threshold Point 5/tempFamOut.csv", "w")
    for aline in linesTogo:
        f.write(aline)
    f.close

def makePWMExcel(families, sig):
    linesTogo = []
    sorts = []
    startSort = sortPWMsByHitsNormalized(families, 7)
    for aSkip in range(8): sorts.append(sortPWMsByHitsAtSkipNormalized(families, aSkip))
    for hit in startSort:
        pwmName = hit[0]
        skipList = []
        for aSkip in range(8):
            flag = False
            for aHit in sorts[aSkip]:
                #print "name", aHit[0]
                if pwmName == aHit[0]:
                    #print "called"
                    flag = True
                    skipList.append(aHit[1])
            if flag == False:
                skipList.append(0)
                print "not found", family, aHit[0]
        for fkey in families.keys():
            if pwmName in families[fkey]:
                family = fkey
                pwmRealName = families[fkey][pwmName]["name"]
                commaSplit = family.split(",")
        familynewname = family        
        if len(commaSplit) > 1:
            newstr = ""
            for strindex in range(len(commaSplit[1])-1):
                newstr = newstr + commaSplit[1][strindex+1]
            familynewname = commaSplit[0] + ":" + newstr
        cmutTotal = 0
        for aSkipIndex in range(len(skipList)):
            cmutTotal += len(families[family][pwmName][str(aSkipIndex)]["sigVals"])
        alineTogo = familynewname + "," + pwmName + "," + pwmRealName + "," + str(families[family][pwmName][str(0)]["length"]+1) + "," + \
                    str(cmutTotal) + "," + str(hit[1])
               
        for aSkipIndex in range(len(skipList)):
            #if families[family][pwmName][str(skip)].has_key("length"):
            alineTogo = alineTogo + "," + str(len(families[family][pwmName][str(aSkipIndex)]["sigVals"])) + "," + str(skipList[aSkipIndex])
            #else:
            #    alineTogo = alineTogo + "," + str(0) + "," + str(aSkipScore)
            
        linesTogo.append(alineTogo + "\n")
    if sig == str(1):    
        f = open("../../Results/Summaries/Top One Percent Same Site Random/tempPWMOut.csv", "w")
    elif sig == str(2):    
        f = open("../../Results/Summaries/One Percent Real/tempPWMOut.csv", "w")
    else:
        f = open("../../Results/Summaries/Threshold Point 5/tempPWMOut.csv", "w")
    for aline in linesTogo:
        f.write(aline)
    f.close

def findName(pwm):
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

def topFiveZPWMs(families):
    f = open("../../Results/Real/topFivePWM.txt", "w")
    for fkey in families.keys():
        f.write(fkey+"\n")
        for pwm in families[fkey].keys():
            pwmList = []
            for s in families[fkey][pwm].keys():
                if s != "name":
                    for v in families[fkey][pwm][s]["sigVals"]:
                        pwmList.append([s, v["zscore"]])
            sort(pwmList)
            f.write("-> " + pwm + "\n")
            if len(pwmList) > 5:
                for i in range(5):
                    f.write("---> " + pwmList[len(pwmList)-i-1][0] + " " + pwmList[len(pwmList)-i-1][1]  + "\n")
    f.close()

def topFiveZPWMsforSkip(families, s):
    f = open("../../Results/Real/topFivePWM " + s + ".txt", "w")
    for fkey in families.keys():
        f.write(fkey+"\n")
        for pwm in families[fkey].keys():
            pwmList = []            
            for v in families[fkey][pwm][s]["sigVals"]:
                pwmList.append([s, v["zscore"]])
            sort(pwmList)
            f.write("-> " + pwm + "\n")
            if len(pwmList) > 5:
                for i in range(5):
                    f.write("---> " + pwmList[len(pwmList)-i-1][0] + " " + pwmList[len(pwmList)-i-1][1]  + "\n")
            else:
                for i in range(len(pwmList)):
                    f.write("---> " + pwmList[len(pwmList)-i-1][0] + " " + pwmList[len(pwmList)-i-1][1]  + "\n")
    f.close()

def topFiveZFamilies(families):
    f = open("../../Results/Real/topFiveFam.txt", "w")
    for fkey in families.keys():
        famList = []
        for pwm in families[fkey].keys():
            for s in families[fkey][pwm].keys():
                if s != "name":
                    for v in families[fkey][pwm][s]["sigVals"]:
                        famList.append([s, v["zscore"]])
        sort(famList)
        f.write(fkey+"\n")
        if len(famList) > 5:
            for i in range(5):
                f.write("->" + famList[len(famList)-i-1][0] + " " + famList[len(famList)-i-1][1] + "\n")
    f.close()

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

#this needs an argument
#1 for one percent random
#5 for .5 threshold on real data
def main(args):
    if args[1] == str(1):
        sigFile = "../../Results/Real/Sig Vals/One Percent Same Site Random/Gap "
    elif args[1] == str(2):
        sigFile = "../../Results/Real/Sig Vals/One Percent Real/Gap "
    elif args[1] == str(5):
        sigFile = "../../Results/Real/Sig Vals/Threshold/Point 5 Threshold/Gap "
    f = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwmList = pickle.load(f)
    f.close()
    pwmInfo = dict()
    for p in pwmList:
        pwmInfo[p] = dict()
        pwmInfo[p]["length"] = findPWMLength(p)
    
##    pwmInfo = dict()
##    for i in range(112):
##        if i < 10:
##            pwmInfo["MA000"+str(i)] = dict()
##            #pwmInfo["MA000"+str(i)]["name"] = findName("MA000"+str(i))
##        elif i < 100:
##            pwmInfo["MA00"+str(i)] = dict()
##            #pwmInfo["MA00"+str(i)]["name"] = findName("MA00"+str(i))
##        else:
##            pwmInfo["MA0"+str(i)] = dict()
##            #pwmInfo["MA0"+str(i)]["name"] = findName("MA0"+str(i))

    #go through and grab all the significant values
    #from the results file
    for skip in range(8):
        f = open(sigFile + str(skip) + ".sig", "r")
        s = f.readline()
        while s != "":
            splitLine = s.split()
            pwmID = splitLine[0]
            pos = splitLine[2] + " " + splitLine[3]
            LL = splitLine[4]
            zscore = splitLine[5]
            if pwmInfo[pwmID].has_key(str(skip)):
                listLen = len(pwmInfo[pwmID][str(skip)])
                pwmInfo[pwmID][str(skip)].append(dict())
                pwmInfo[pwmID][str(skip)][listLen]["pos"] = pos
                pwmInfo[pwmID][str(skip)][listLen]["LL"] = LL
                pwmInfo[pwmID][str(skip)][listLen]["zscore"] = zscore
            else:
                pwmInfo[pwmID][str(skip)] = []
                pwmInfo[pwmID][str(skip)].append(dict())
                pwmInfo[pwmID][str(skip)][0]["pos"] = pos
                pwmInfo[pwmID][str(skip)][0]["LL"] = LL
                pwmInfo[pwmID][str(skip)][0]["zscore"] = zscore
            s = f.readline()
        f.close()
    
    families = dict()
    count = 0
    f = open("../../Doc/Vertebrate PWM/full_matrix_list.txt", "r")
    s = f.readline()
    while s != "":
        stringList = s.split(";")[0].split()
        
        classString = ""
        for i in range(len(stringList)):
            if i > 2:
                classString += stringList[i]
                if i < len(stringList)-1:
                    classString += "_"
        if families.has_key(classString): pass
        else: families[classString] = dict()
        #print classString
        s = f.readline()
        count += 1
    f.close()
    print "found all", count, "PWMs"

    f = open("../../Doc/Vertebrate PWM/full_matrix_list.txt", "r")
    s = f.readline()
    while s != "":
        colonSplit = s.split(";")
        stringList = colonSplit[0].split()
        
        classString = ""
        pwmID = stringList[0]
        for i in range(len(stringList)):
            if i > 2:
                classString += stringList[i]
                if i < len(stringList)-1:
                    classString += "_"
        for line in colonSplit:
            lineSplit = line.split()
            if lineSplit[0] == "sysgroup":
                if lineSplit[1] == "\"vertebrate\"":
                    families[classString][pwmID] = dict()
                    families[classString][pwmID]["name"] = findName(pwmID)
                    for skip in range(8):
                        families[classString][pwmID][str(skip)] = dict()
                        families[classString][pwmID][str(skip)]["pwmHits"] = 0
                        families[classString][pwmID][str(skip)]["sigVals"] = []
                        #info about all the branch/branch changes found for this PWM at this skip
                        families[classString][pwmID][str(skip)]["16sq"] = []
                    break
                    
        s = f.readline()
        
    f.close()
    #kill all the families with no vertebrate PWMs
    for fkey in families.keys():
        if len(families[fkey]) == 0: del families[fkey]
        
    #pwmHitFile = open("LL all PWM pval -9 NA 0/MA0002(p -9).gff.ll", "r")
    #print "done"
    #sys.exit(0)
    for fkey in families.keys(): #for every family
        for pwmKey in families[fkey].keys(): #get every pwm for this family            
            for skipKey in families[fkey][pwmKey].keys(): #look at all skip positions
                if skipKey != "name":
##                    if pwmKey == "MA0024":
##                            print "have key", skipKey
                    flag = False
                    try:                    
                        pwmHitFile = open("../../Results/Real/Raw Data/Gap " + skipKey + "/" + pwmKey + "(p -9).gff.ll", "r")
                        flag = True
                    except:
                        if skipKey == "0":
                            #remove this PWM b/c it has no hits
                            del families[fkey][pwmKey]
                            #no need to look at the other skips
                            break                         
                         
                    if flag:                        
                        #fill in the sig hits for this pwm at this skip
                        #why would a PWM not have a skipKey:
                            #the PWM is not long enough
                            #or the PWM had no sig hits at this skip
                        if pwmInfo[pwmKey].has_key(skipKey):
                            for infoDict in pwmInfo[pwmKey][skipKey]:
                                families[fkey][pwmKey][skipKey]["sigVals"].append(infoDict)
    ##                    families[fkey][pwmKey][skipKey]["sigVals"][i]["pos"] = s.split(")")[0]+")"
    ##                    families[fkey][pwmKey][skipKey]["sigVals"][i]["LL"] = s.split()[2]
    ##                    families[fkey][pwmKey][skipKey]["sigVals"][i]["zscore"] = s.split()[3]
                            i = 0
                            s = pwmHitFile.readline()                    
                            while s.split()[0] != "Hit":                        
                                s = pwmHitFile.readline()
                                i += 1
                            pwmHitFile.close()                    
                            families[fkey][pwmKey][skipKey]["pwmHits"] = int(s.split()[2])
                            
                            #this is just the number of possible values for the PWM
                            if families[fkey][pwmKey][skipKey].has_key("length"):
                                if i != families[fkey][pwmKey][skipKey]["length"]:
                                    print "error in lengths:", pwmKey                            
                            else:
                                if pwmKey == "MA0024" and skipKey == str(0):
                                    print "assigned" 
                                families[fkey][pwmKey][skipKey]["length"] = i
                        else:
                            i = 0
                            s = pwmHitFile.readline()                    
                            while s.split()[0] != "Hit":                        
                                s = pwmHitFile.readline()
                                i += 1
                            pwmHitFile.close()                    
                            families[fkey][pwmKey][skipKey]["pwmHits"] = int(s.split()[2])
                            
                            #this is just the number of possible values for the PWM
                            if families[fkey][pwmKey][skipKey].has_key("length"):
                                if i != families[fkey][pwmKey][skipKey]["length"]:
                                    print "error in lengths:", pwmKey                            
                            else:
                                families[fkey][pwmKey][skipKey]["length"] = i 
##                    if flag == True and int(skipKey) < 2:                    
##                        pwm16sqFile = open("LL all PWM pval -9 NA " + skipKey + "/" + pwmKey + "(p -9).gff.256sq", "r")
##                        listForAllPos = pickle.load(pwm16sqFile)
##                        for infoKey in families[fkey][pwmKey][skipKey]["sigVals"]:
##                            #print infoKey
##                            families[fkey][pwmKey][skipKey]["16sq"].append(infoKey["pos"][1])                                            
            if families[fkey].has_key(pwmKey):
                if families[fkey][pwmKey]["name"] == "E2F":
                    print pwmKey
                    for s in families[fkey][pwmKey].keys():
                        if s != "name" and int(s) <= int(families[fkey][pwmKey][str(0)]["length"])-1:
                            print s, families[fkey][pwmKey][s]["length"]
##    print families["RUNT"]["MA0002"]["0"]["pwmHits"]
##    for sigVal in families["RUNT"]["MA0002"]["0"]["sigVals"]:
##        print sigVal["pos"], sigVal["LL"], sigVal["zscore"]
    #sortPWMsByHits(families)
##    myans = []
##    for s in range(8):
##        myans.append(str(s) + " Sites Away Total\n")
##        ans = sortFamiliesByHitsAtSkip(families, s)
##        for a in ans:
##            myans.append(str(a)+"\n")
##        ans = sortFamiliesByHitsAtSkipNormalized(families, s)
##        myans.append("\n" + str(s) + " Sites Away Normalized\n")
##        for a in ans:
##            myans.append(str(a)+"\n")
##        myans.append("\n")
##        
##    outFile = open("CmutsCounts.txt", "w")
##    for a in myans:
##        outFile.write(a)
##    outFile.close()
    
##    print "At 0"
##    ans = sortFamiliesByHitsAtSkipNormalized(families, 8)
##    for a in ans: print a
##    for pwmkey in families["Unknown"]:
##        print pwmkey
##    for v in families["Unknown"]["MA0024"]["0"]["sigVals"]: print v


##    ans = sortPWMsByHits(families, 7)
##    for a in ans:
##        for fkey in families.keys():
##            if a[0] in families[fkey].keys():
##                print a, fkey
##    ans = sortPWMsByHitsNormalized(families, 7)
##    for a in ans:
##        for fkey in families.keys():
##            if a[0] in families[fkey].keys():
##                print a, fkey
    
    #need to find the distance distribution per family
    #so many adjacent, so many one away, and so on
##    famDict = distanceDistributionsByFamily(families)
##    famDict2 = distanceDistributionsByFamilyNormalized(families)
##    for f in famDict.keys():
##        print f
##        for t in range(len(famDict[f])):
##            print t, famDict[f][str(t)], famDict2[f][str(t)]    

##    pwmDict = distanceDistributionsByPWMNormalized(families)
##    pwmDict2 = distanceDistributionsByPWM(families)
##    for pwm in families["bZIP"].keys():
##        print pwm
##        for t in range(len(pwmDict[pwm])):        
##            print t, pwmDict[pwm][str(t)], pwmDict2[pwm][str(t)]

##    skipDist = distanceDistByFamilyPWMStyleNormalized(families)
##    family = "RUNT"
##    for pwm in families[family].keys(): print pwm
##    for s in skipDist[family].keys():
##        print "skip", s
##        for d in skipDist[family][str(s)].keys():
##            print "  ", d, skipDist[family][str(s)][str(d)]

##    for fkey in families.keys():
##        print fkey
##        for s in range(2):
##            print "Skip", s
##            for pwmkey in families[fkey].keys():
####                print families[fkey][pwmkey][str(s)]["16sq"]
####                sys.exit(0)
##                for pos in families[fkey][pwmkey][str(s)]["16sq"]:
##                    print pos
##                    ls = []
##                    print families[fkey][pwmkey][str(s)]["16sq"][int(pos)]
##                    for dkey in families[fkey][pwmkey][str(s)]["16sq"][int(pos)]:
##                        #print dkey
##                        ls.append([dkey, families[fkey][pwmkey][str(s)]["16sq"][int(pos)][dkey]])
##                        sort(ls)
##                print pwmkey, len(ls)
##                for index in range(5): print "\t", ls[index][0], ls[index][1]
    makeFamilyExcel(families, args[1])
    makePWMExcel(families, args[1])
    #topFiveZPWMs(families)
##    for s in range(5):
##        topFiveZPWMsforSkip(families, str(s))
##    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
