import sys, newPlot

def main(args):
    skip = range(9)    
    pvalList = [-9]    
    for j in pvalList:         
        pvals = []            
        newPlot.binnedPlot("../../Results/Controls/Same Site/Binned Data/same site Gap " + str(0) + ".binned", \
                           ["../../Results/Controls/Same Site/Binned Data/same site Gap " + str(1) + ".binned",\
                            "../../Results/Controls/Same Site/Binned Data/same site Gap " + str(2) + ".binned",\
                            "../../Results/Controls/Same Site/Binned Data/same site Gap " + str(3) + ".binned",\
                            "../../Results/Controls/Same Site/Binned Data/same site Gap " + str(4) + ".binned",\
                            "../../Results/Controls/Same Site/Binned Data/same site Gap " + str(5) + ".binned",\
                            "../../Results/Controls/Same Site/Binned Data/same site Gap " + str(6) + ".binned",\
                            "../../Results/Controls/Same Site/Binned Data/same site Gap " + str(7) + ".binned"],
                           pvals,\
                           ["green", "red", "black", "orange", "blue", "yellow", "grey", "brown"],\
                           "../../Results/Graphs/Line/All Gaps Same Site Random.png")            
    sys.exit(0)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
