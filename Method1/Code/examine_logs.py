import sys
import pickle

def main(args):
    f = open("REALGOOD_random_logs.txt", "r")
    p = pickle.load(f)
    for l in p: print l
    f.close()
    sys.exit(0)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
