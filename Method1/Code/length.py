import sys, pickle

def main(args):
    f = open("Normal log pval -9 NA 7.txt","r")
    a = pickle.load(f)
    f.close()
    sum = 0
    keys = a.keys()
    for k in keys:
        sum += a[k]
    print sum
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
