import sys, pickle

def length(f):
    count = 0
    s = f.readline()
    while s != "":
        count += 1
        s = f.readline()
    f.close()
    return count

def length2(d):
    count = 0
    for a in d.keys():
        count += d[a]
    return count

def main(args):
    for s in range(8):
        all = open("../../Results/Real/Binned Data/Gap " + str(s) + ".binned", "r")
        a = pickle.load(all)
        all.close()
        allLength = float(length2(a))
        sig = open("../../Results/Real/Sig Vals/Threshold/Point 5 Threshold/Gap " + str(s) + ".sig", "r")
        c = float(length(sig))
        print c/allLength

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
