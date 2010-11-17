import sys

def main(args):
    f = open(args[1], "r")

    count = 0
    s = f.readline()
    while s != "":
        count += 1
        s = f.readline()
    print count
    sys.exit(0)
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
