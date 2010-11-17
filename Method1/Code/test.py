import sys, pickle
from rpy import *

def main(args):
    nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    print r.avg(nums)
    sys.exit(0)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
