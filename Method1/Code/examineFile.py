import sys, pickle

def main(args):
    myNucs = ["A", "C", "G", "T"]
    myEdges = []
    for n1 in myNucs:
        for n2 in myNucs:
            myEdges.append(n1+n2)
    f = open("edgeCounts.txt", "r")
    edgeCountMatrix = pickle.load(f)
    f.close()
    problemCount = 0
    for e1 in myEdges:
        for e2 in myEdges:
            if edgeCountMatrix[e1+e2] == 0:
                print e1+e2
                problemCount += 1
            #else: print e1+e2, edgeCountMatrix[e1+e2]
    print problemCount
    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
