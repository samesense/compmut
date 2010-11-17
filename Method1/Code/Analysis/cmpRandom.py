import pickle
from rpy import *

def main(args):
    for s in range(5):
        r1f = open("random_logs NA " + str(s) + " p -9.txt", "r")
        r2f = open("random_logs NA " + str(s) + " p -9 cmp.txt", "r")
        r1 = pickle.load(r1f)
        r2 = pickle.load(r2f)
        print s, r.wilcox_test(r1, r2, alternative="two.sided", paired=False)["p.value"]

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
