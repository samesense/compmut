import pickle
from rpy import *

def myMin(theList):
    m = float(10)
    for item in theList:
        if float(item) < float(m):            
            m = float(item)
    return m

def myMax(theList):
    m = float(-1)
    for item in theList:
        if float(item) > float(m):            
            m = float(item)
    return m

def dotPlot(realDistFile, randomDistFile, outputFile):
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)
    allYVals = []
    xReal = []
    yReal = []
    #print len(p)    
    keys = p.keys()
    #if p.has_key("0.0"): print p["0.0"]
    maxXReal = myMax(keys)
    minXReal = myMin(keys)
            
    for i in range(len(keys)):
        allYVals.append(p[keys[i]])
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1    

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open(randomDistFile, "r")
    p = pickle.load(theFile)
    #print p
    xRand = []
    yRand = []
    #print len(p)
    keys = p.keys()
    map(float, keys)
    keys.sort()
    #print keys
    maxXRand = myMax(keys)
    minXRand = myMin(keys)
    for i in range(len(keys)):
        allYVals.append(p[keys[i]])
        xRand.append(keys[i])
        yRand.append(p[keys[i]])    
    
    #r.plot(xReal, yReal, col="green", xlab="Log Likelihood", ylab="Count", ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), max([maxXReal, maxXRand])))
    r.plot(xReal, yReal, col="green", xlab="Log Likelihood", ylab="Count", \
           ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), max([maxXReal, maxXRand])))
    r.points(xRand, yRand, col="black")   
    r.dev_off()
    theFile.close()

    print sum(yReal), sum(yRand)

def dotPlot2(realDistFile, randomDistFile, outputFile):
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)
    allYVals = []
    xReal = []
    yReal = []
    #print len(p)    
    keys = p.keys()
    #if p.has_key("0.0"): print p["0.0"]
    maxXReal = myMax(keys)
    minXReal = myMin(keys)
            
    for i in range(len(keys)):
        allYVals.append(p[keys[i]])
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1    

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open(randomDistFile, "r")
    p = pickle.load(theFile)
    #print p
    xRand = []
    yRand = []
    #print len(p)
    keys = p.keys()
    map(float, keys)
    keys.sort()
    #print keys
    maxXRand = myMax(keys)
    minXRand = myMin(keys)
    for i in range(len(keys)):
        allYVals.append(p[keys[i]])
        xRand.append(keys[i])
        yRand.append(p[keys[i]])    
    
    #r.plot(xReal, yReal, col="green", xlab="Log Likelihood", ylab="Count", ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), max([maxXReal, maxXRand])))
    r.plot(xReal, yReal, col="green", xlab="Log Likelihood", ylab="Count", \
           ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), max([maxXReal, maxXRand])))
    r.points(xRand, yRand, col="black")   
    r.dev_off()
    theFile.close()

    print sum(yReal), sum(yRand)

def binPlot(realDistFile, randomDistFile, outputFile):
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)    
    xReal = []
    yReal = []
    allYVals = []
    #print len(p)
    #print p
    keys = p.keys()        
    #print keys, myMin(keys)
    maxXReal = myMax(keys)
    minXReal = myMin(keys)
            
    for i in range(len(keys)):        
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1    

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open(randomDistFile, "r")
    p = pickle.load(theFile)
    #print p
    xRand = []
    yRand = []
    #print len(p)
    keys = p.keys()   
    maxXRand = myMax(keys)
    minXRand = myMin(keys)
    for i in range(len(keys)):        
        xRand.append(keys[i])
        yRand.append(p[keys[i]])

    increment = (max([maxXReal, maxXRand]) - min([minXReal,minXRand])) / float(20)
    mult = range(21)
    domain = []
    for m in mult: domain.append(min([minXReal,minXRand]) + float(m)*increment)
    realRange = []
    randRange = []
        
    for d in domain:
        totalRealVal = float(0)
        totalRandVal = float(0)
        for index in range(len(xReal)):
            if float(xReal[index]) >= float(d) and float(xReal[index]) < float(d) + float(increment): totalRealVal += float(yReal[index])
        for index in range(len(xRand)):
            if float(xRand[index]) >= float(d) and float(xRand[index]) < float(d) + float(increment): totalRandVal += float(yRand[index])
        allYVals.append(totalRealVal)
        allYVals.append(totalRandVal)
        realRange.append(totalRealVal)
        randRange.append(totalRandVal)
    r.plot(domain, realRange, col="green", type = "l", xlab="Log Likelihood", \
           ylab="Count", ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), \
                                                       max([maxXReal, maxXRand])))
    r.lines(domain, randRange, col="black")   
    r.dev_off()
    theFile.close()

#give me the set of distributions to plot, with the real distribution given first
#the pvals to display for each control in the distribution
#the colors for each distribution
#the output file
def binnedPlot(realDistFile, controlFiles, pvals, colors, outputFile):
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)
    theFile.close()
    xReal = []
    yReal = []
    allYVals = []
    keys = p.keys()
    maxX = myMax(keys)
    minX = myMin(keys)
            
    for i in range(len(keys)):        
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1
     
    xControl = []
    yControl = []
    for aDistFile in range(len(controlFiles)):
        theFile = open(controlFiles[aDistFile], "r")
        p = pickle.load(theFile)
        theFile.close()
        xControl.append([])
        yControl.append([])
        keys = p.keys()
        aMax = myMax(keys)
        aMin = myMin(keys)
        if aMax > maxX:
            maxX = aMax
        if aMin < minX:
            minX = aMin
        for i in range(len(keys)):
            xControl[len(xControl)-1].append(keys[i])
            yControl[len(xControl)-1].append(p[keys[i]])
        
    increment = (maxX - minX) / float(10)
    mult = range(11)
    domain = []
    for m in mult: domain.append(minX + float(m)*increment)
    realRange = []
    controlRange = []
    for c in controlFiles: controlRange.append([])    
        
    for d in domain:
        totalRealVal = float(0)
        controlTotal = []
        for c in controlFiles: controlTotal.append(float(0))
        
        for index in range(len(xReal)):
            if float(xReal[index]) >= float(d) and float(xReal[index]) < float(d) + float(increment):
                totalRealVal += float(yReal[index])
        for cIndex in range(len(controlFiles)):
            for index in range(len(xControl[cIndex])): 
                if float(xControl[cIndex][index]) >= float(d) and float(xControl[cIndex][index]) < float(d) + float(increment):
                    controlTotal[cIndex] += float(yControl[cIndex][index])
        
        allYVals.append(totalRealVal)
        for cIndex in range(len(controlFiles)):
            allYVals.append(controlTotal[cIndex])
        
        realRange.append(totalRealVal)
        for cIndex in range(len(controlFiles)):
            controlRange[cIndex].append(controlTotal[cIndex])
        
    r.plot(domain, realRange, col=colors[0], type = "l", xlab="Log Likelihood", ylab="Count", \
           ylim=(0,max(allYVals)), xlim=(minX, maxX))
    for cIndex in range(len(controlFiles)):
        r.lines(domain, controlRange[cIndex], col=colors[cIndex+1])

    if len(pvals) != 0:
        r.title("Gap " + outputFile[len(outputFile)-5])            
    else:
        r.title("All Gaps")
    #display the pvals
    if len(pvals) == 0:
        r.legend(maxX-.1, max(allYVals)-5, ("Gap 0", "Gap 1", "Gap 2", "Gap 3", "Gap 4", "Gap 5", "Gap 6", "Gap 7"), \
                     col = (colors[0], colors[1], colors[2], colors[3], colors[4], colors[5], colors[6], colors[7]), pch = (0, 0))
    elif len(pvals) == 4:
        if int(outputFile[len(outputFile)-5]) == 8:
            r.legend(maxX-.2, max(allYVals)-5, ("Same Site\n" + str(pvals[0]), "\nConstructed PWMs\n" + str(pvals[1]), \
                                "\nPermuted PWMs\n" + str(pvals[2]), "\nRandom Sites " + str(pvals[3])), \
                     col = (colors[1], colors[2], colors[3], colors[4]), pch = (0, 0))
        elif int(outputFile[len(outputFile)-5]) < 4:
            r.legend(maxX-.35, max(allYVals)-5, ("Same Site\n" + str(pvals[0]), "\nConstructed PWMs\n" + str(pvals[1]), \
                                "\nPermuted PWMs\n" + str(pvals[2]), "\nRandom Sites " + str(pvals[3])), \
                     col = (colors[1], colors[2], colors[3], colors[4]), pch = (0, 0))
        else:
            r.legend(maxX-.25, max(allYVals)-5, ("Same Site\n" + str(pvals[0]), "\nConstructed PWMs\n" + str(pvals[1]), \
                                "\nPermuted PWMs\n" + str(pvals[2]), "\nRandom Sites " + str(pvals[3])), \
                     col = (colors[1], colors[2], colors[3], colors[4]), pch = (0, 0))
    elif len(pvals) == 3:
        r.legend(0.3, 200, (str(pvals[0]), str(pvals[1]), str(pvals[3])), col = (colors[1], colors[2], colors[3]), pch = (0, 0))
    elif len(pvals) == 2:
        r.legend(0.3, 200, (str(pvals[0]), str(pvals[1])), col = (colors[1], colors[2]), pch = (0, 0))
    else:
        r.legend(0.3, 200, (str(pvals[0])), col = (colors[1]), pch = (0, 0))
    r.dev_off()

#give me the set of distributions to plot, with the real distribution given first
#the pvals to display for each control in the distribution
#the colors for each distribution
#the output file
def dottedPlot(realDistFile, controlFiles, pvals, colors, outputFile):
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)
    theFile.close()
    xReal = []
    yReal = []
    allYVals = []
    keys = p.keys()
    maxX = myMax(keys)
    minX = myMin(keys)
    allYVals = []
            
    for i in range(len(keys)):        
        xReal.append(keys[i])
        yReal.append(p[keys[i]])
        allYVals.append(p[keys[i]])
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1
     
    xControl = []
    yControl = []
    for aDistFile in range(len(controlFiles)):
        theFile = open(controlFiles[aDistFile], "r")
        p = pickle.load(theFile)
        theFile.close()
        xControl.append([])
        yControl.append([])
        keys = p.keys()
        aMax = myMax(keys)
        aMin = myMin(keys)
        if aMax > maxX:
            maxX = aMax
        if aMin < minX:
            minX = aMin
        for i in range(len(keys)):
            xControl[len(xControl)-1].append(keys[i])
            yControl[len(xControl)-1].append(p[keys[i]])
            allYVals.append(p[keys[i]])
        
    r.plot(xReal, yReal, col=colors[0], xlab="Log Likelihood", ylab="Count", \
           ylim=(0,max(allYVals)), xlim=(minX, maxX))
    for cIndex in range(len(controlFiles)):
        r.points(xControl[cIndex], yControl[cIndex], col=colors[cIndex+1])

    r.title("Gap " + outputFile[len(outputFile)-5])            

    #display the pvals
    if len(pvals) == 4:        
        r.legend(0.3, max(allYVals)-1, ("Same Site\n" + str(pvals[0]), "\nConstructed PWMs\n" + str(pvals[1]), \
                            "\nPermuted PWMs\n" + str(pvals[2]), "\nRandom Sites " + str(pvals[3])), \
                 col = (colors[1], colors[2], colors[3], colors[4]), pch = (0, 0))
    elif len(pvals) == 3:
        r.legend(0.3, 200, (str(pvals[0]), str(pvals[1]), str(pvals[3])), col = (colors[1], colors[2], colors[3]), pch = (0, 0))
    elif len(pvals) == 2:
        r.legend(0.3, 200, (str(pvals[0]), str(pvals[1])), col = (colors[1], colors[2]), pch = (0, 0))
    else:
        r.legend(0.3, 200, (str(pvals[0])), col = (colors[1]), pch = (0, 0))
    r.dev_off()
    
def binPlot3(realDistFile, ccDistFile, randomDistFile, veryRandomFile, outputFile):
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)    
    xReal = []
    yReal = []
    allYVals = []
    #print len(p)
    #print p
    keys = p.keys()        
    #print keys, myMin(keys)
    maxXReal = myMax(keys)
    minXReal = myMin(keys)
            
    for i in range(len(keys)):        
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1    

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open(randomDistFile, "r")
    p = pickle.load(theFile)
    #print p
    xRand = []
    yRand = []
    #print len(p)
    keys = p.keys()   
    maxXRand = myMax(keys)
    minXRand = myMin(keys)
    for i in range(len(keys)):        
        xRand.append(keys[i])
        yRand.append(p[keys[i]])

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open(veryRandomFile, "r")
    p = pickle.load(theFile)
    #print p
    xVRand = []
    yVRand = []
    #print len(p)
    keys = p.keys()   
    maxXVRand = myMax(keys)
    minXVRand = myMin(keys)
    for i in range(len(keys)):        
        xVRand.append(keys[i])
        yVRand.append(p[keys[i]])

    theFile.close()
    theFile = open(ccDistFile, "r")
    p = pickle.load(theFile)
    #print p
    xCC = []
    yCC = []
    #print len(p)
    keys = p.keys()   
    maxXCC = myMax(keys)
    minXCC = myMin(keys)
    for i in range(len(keys)):        
        xCC.append(keys[i])
        yCC.append(p[keys[i]])

    increment = (max([maxXReal, maxXRand, maxXCC, maxXVRand]) - min([minXReal,minXRand, minXCC, minXVRand])) / float(10)
    mult = range(11)
    domain = []
    for m in mult: domain.append(min([minXReal,minXRand, minXCC, minXVRand]) + float(m)*increment)
    realRange = []
    randRange = []
    vrandRange = []
    ccRange = []
        
    for d in domain:
        totalRealVal = float(0)
        totalRandVal = float(0)
        totalCCVal = float(0)
        totalVRandVal = float(0)
        for index in range(len(xReal)):
            if float(xReal[index]) >= float(d) and float(xReal[index]) < float(d) + float(increment): totalRealVal += float(yReal[index])
        for index in range(len(xRand)):
            if float(xRand[index]) >= float(d) and float(xRand[index]) < float(d) + float(increment): totalRandVal += float(yRand[index])
        for index in range(len(xVRand)):
            if float(xVRand[index]) >= float(d) and float(xVRand[index]) < float(d) + float(increment): totalVRandVal += float(yVRand[index])
        for index in range(len(xCC)):
            if float(xCC[index]) >= float(d) and float(xCC[index]) < float(d) + float(increment): totalCCVal += float(yCC[index])
        allYVals.append(totalRealVal)
        allYVals.append(totalRandVal)
        allYVals.append(totalVRandVal)
        allYVals.append(totalCCVal) 
        realRange.append(totalRealVal)
        randRange.append(totalRandVal)
        ccRange.append(totalCCVal)
        vrandRange.append(totalVRandVal)
    r.plot(domain, realRange, col="green", type = "l", xlab="Log Likelihood", ylab="Count", ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), max([maxXReal, maxXRand])))
    r.lines(domain, randRange, col="black")
    r.lines(domain, ccRange, col="orange")
    r.lines(domain, vrandRange, col="red")
    r.dev_off()
    theFile.close()

    ysum1 = 0
    for y in yReal:
        ysum1 += y
    ysum2 = 0
    for y in yRand:
        ysum2 += y
    ysum3 = 0
    for y in yCC:
        ysum3 += y
    print ysum1, ysum2, ysum3
    

def binPlot2(distFiles, outputFile):
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1 

    domains = []
    ranges = []
    for d in distFiles:
        theFile = open(d, "r")
        p = pickle.load(theFile)    
        x = []
        y = []
    allYVals = []
    #print len(p)
    #print p
    keys = p.keys()        
    #print keys, myMin(keys)
    maxXPWM = myMax(keys)
    minXPWM = myMin(keys)
            
    for i in range(len(keys)):        
        xPWM.append(keys[i])
        yPWM.append(p[keys[i]])   
       

    theFile.close()

    theFile = open(ccFile, "r")
    p = pickle.load(theFile)    
    xReal = []
    yReal = []
    allYVals = []
    #print len(p)
    #print p
    keys = p.keys()        
    #print keys, myMin(keys)
    maxXReal = myMax(keys)
    minXReal = myMin(keys)
            
    for i in range(len(keys)):        
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1    

    theFile.close()
    
    theFile = open(realDistFile, "r")
    p = pickle.load(theFile)    
    xReal = []
    yReal = []
    allYVals = []
    #print len(p)
    #print p
    keys = p.keys()        
    #print keys, myMin(keys)
    maxXReal = myMax(keys)
    minXReal = myMin(keys)
            
    for i in range(len(keys)):        
        xReal.append(keys[i])
        yReal.append(p[keys[i]])   
    r.png(outputFile, width=733, height=900)
    r.par(pty="s")
    r.do_points=1    

    theFile.close()
    #theFile = open("rdist.txt", "r")
    theFile = open(randomDistFile, "r")
    p = pickle.load(theFile)
    #print p
    xRand = []
    yRand = []
    #print len(p)
    keys = p.keys()   
    maxXRand = myMax(keys)
    minXRand = myMin(keys)
    for i in range(len(keys)):        
        xRand.append(keys[i])
        yRand.append(p[keys[i]])

    increment = (max([maxXReal, maxXRand]) - min([minXReal,minXRand])) / float(20)
    mult = range(21)
    domain = []
    for m in mult: domain.append(min([minXReal,minXRand]) + float(m)*increment)
    realRange = []
    randRange = []
        
    for d in domain:
        totalRealVal = float(0)
        totalRandVal = float(0)
        for index in range(len(xReal)):
            if float(xReal[index]) >= float(d) and float(xReal[index]) < float(d) + float(increment): totalRealVal += float(yReal[index])
        for index in range(len(xRand)):
            if float(xRand[index]) >= float(d) and float(xRand[index]) < float(d) + float(increment): totalRandVal += float(yRand[index])
        allYVals.append(totalRealVal)
        allYVals.append(totalRandVal)
        realRange.append(totalRealVal)
        randRange.append(totalRandVal)
    r.plot(domain, realRange, col="green", type = "l", xlab="Log Likelihood", ylab="Count", \
           ylim=(0,max(allYVals)), xlim=(min([minXReal,minXRand]), max([maxXReal, maxXRand])))
    r.lines(domain, randRange, col="black")   
    r.dev_off()
    theFile.close()       

