import sys, pickle

def main(args):
    f = open(args[1], "r")
    p = pickle.load(f)
    print p
    print len(p)
    
    sys.exit(0)
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
