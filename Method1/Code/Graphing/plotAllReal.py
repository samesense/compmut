import sys, newPlot

def main(args):
    skip = range(9)    
    pvalList = [-9]    
    for j in pvalList:         
        pvals = []            
        newPlot.binnedPlot("../../Results/Real/Binned Data/Gap " + str(0) + ".binned", \
                           ["../../Results/Real/Binned Data/Gap " + str(1) + ".binned",\
                            "../../Results/Real/Binned Data/Gap " + str(2) + ".binned",\
                            "../../Results/Real/Binned Data/Gap " + str(3) + ".binned",\
                            "../../Results/Real/Binned Data/Gap " + str(4) + ".binned",\
                            "../../Results/Real/Binned Data/Gap " + str(5) + ".binned",\
                            "../../Results/Real/Binned Data/Gap " + str(6) + ".binned",\
                            "../../Results/Real/Binned Data/Gap " + str(7) + ".binned"],
                           pvals,\
                           ["green", "red", "black", "orange", "blue", "yellow", "grey", "brown"],\
                           "../../Results/Graphs/Line/All Gaps Real " + ".png")            
    sys.exit(0)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
