#!/usr/bin/python
# Greg Donahue, 09/23/2005
# -----------------------------------------------------------------------------
import sys, re, pickle, os, random, chisquarematrix
from Numeric import array, zeros
# -----------------------------------------------------------------------------
# GLOBALS
nucleotides = [ "A", "T", "G", "C" ]
costs = dict()
for i in range(len(nucleotides)):
    for j in range(len(nucleotides)):
        costs[nucleotides[i]+nucleotides[j]] = -1
costs["AA"] = -2
costs["TT"] = -2
costs["GG"] = -2
costs["CC"] = -2

# CLASSES
# Class TreeNode represents a single tree node
class TreeNode:

    # Instance variables:
    # self.index is the taxon/ancestor ID of this TreeNode
    # self.children is the child set for this TreeNode
    # self.color is the coloring for this TreeNode

    # Instantiate new TreeNodes
    def __init__(self, i):
        self.index = i
        self.children = []
        self.color = ""

    # Return a string representation of this TreeNode
    def __repr__(self):
        return str(self.index)

    # Returns whether this TreeNode is equal to another TreeNode
    def __eq__(self, other):
        return (self.index == other.index)

# Class Edge represents an edge in this tree
class Edge:

    # Instance variables:
    # self.one is the first TreeNode in this edge
    # self.two is the second TreeNode in this edge

    # Instantiate new Edges
    def __init__(self, o, t):
        self.one = o
        self.two = t

    # Return a string representation of this Edge
    def __repr__(self):
        return str(self.one)+"-"+str(self.two)

# Class Tree represents phylogeny trees, sans taxa
class Tree:

    # Instance variables:
    # self.nodes is the set of nodes
    # self.edges is a set of edges
    # self.r is the root TreeNode

    # Instantiate new Trees
    def __init__(self, edgelines):

        # God, this is a horrible kludge
        # Occurs solely when we make a new Tree via the copy function
        if len(edgelines) == 0:
            self.nodes = []
            self.r = TreeNode("r")
            self.nodes.append(self.r)
            self.edges = []
            return
        
        # Parse the edgelines and build Edges
        p = re.compile("(\\([\\d]+\\)[\\*]*[\\w]*[\\t ]+"+\
                       "-[\\t ]+\\([\\d]+\\)[\\*]*[\\w]*)\\n")
        edgeset = re.findall(p, edgelines)
        self.edges = []
        for edge in edgeset:
            p = re.compile("\\(([\\d]+)\\)([\\*]*)[\\w]*")
            nodes = re.findall(p, edge)
            self.edges.append(Edge(TreeNode(nodes[0][0]+nodes[0][1]),
                                   TreeNode(nodes[1][0]+nodes[1][1])))

        # Build a tree structure out of the Edges
        self.nodes = []
        for e in self.edges:
            n_one = e.one
            n_two = e.two
            if not n_one in self.nodes: self.nodes.append(n_one)
            else: n_one = self.nodes[self.nodes.index(n_one)]
            if not n_two in self.nodes: self.nodes.append(n_two)
            else: n_two = self.nodes[self.nodes.index(n_two)]
            n_one.children.append(n_two)
            n_two.children.append(n_one)

        # Root the tree structure
        # Create the root, bind it to the first two TreeNodes
        self.r = TreeNode("r")
        one = self.edges[0].one
        two = self.edges[0].two
        self.r.children.append(one)
        self.r.children.append(two)
        # Bind the first two TreeNodes to it
        one.children.append(self.r)
        two.children.append(self.r)
        # Unbind the first two TreeNodes from each other
        del one.children[one.children.index(two)]
        del two.children[two.children.index(one)]
        # Add root to the list of nodes
        self.nodes.append(self.r)

        # Remove edges back to parents to make it a true Tree
        self.root(self.r.children[0], self.r)
        self.root(self.r.children[1], self.r)
    
    # Return a string representation of this Tree
    def __repr__(self):
        return str(self.edges)

    # Root a subtree (recursive)
    def root(self, node, parent):
        #print node.children, parent
        del node.children[node.children.index(parent)]
        if len(node.children) == 0: return
        for c in node.children:
            self.root(c, node)

    # Clone this tree
    def copy(self, root):
        ret = Tree([])
        ret.r = self.copyR(root)
        ret.nodes = self.fillnodes([], ret.r)
        return ret

    # Recursive clone function
    def copyR(self, node):
        n = TreeNode(node.index)
        if len(node.children) == 0: return n
        for c in node.children:
            n.children.append(self.copyR(c))
        return n

    # Recursively add nodes to the Tree's nodelist
    def fillnodes(self, nodelist, node):
        nodelist.append(node)
        if len(node.children) == 0: return nodelist
        for c in node.children:
            nodelist = self.fillnodes(nodelist, c)
        return nodelist

# MODULE FUNCTIONS
# Translate a sequence to its reverse-strand complement
def translate(sequence):
    ret = ""
    for i in range(len(sequence)):
        j = len(sequence)-i-1
        if sequence[j] == "T": ret += "A"
        elif sequence[j] == "A": ret += "T"
        elif sequence[j] == "G": ret += "C"
        elif sequence[j] == "C": ret += "G"
    return ret

# Load GFF hits from PWM_SCAN, throw out hits that span gaps, build Motifs
def buildMotifs(gff, alignment):
    p = re.compile("([\\w_-]+)[\\t][\\w]+[\\t][\\w]+[\\t]([\\d]+)"+\
                   "[\\t]([\\d]+)[\\t]([\\d.-]+)[\\t]([\\+\\-])")
    hits = re.findall(p, gff.read())
    ret = []
    for h in hits:

        # Collect the sequence elements from the his
        id = h[0]
        start = int(h[1])
        stop = int(h[2])
        score = float(h[3])
        strand = h[4]
        entry = alignment[id]
        
        # Get the correct locus positions
        sequence = entry["human"]
        count = 0
        for j in range(len(sequence)):
            if sequence[j] != "-" and sequence[j] != ".": count += 1
            if count == start:
                count = j
                break
        stop = count+(stop-start)
        start = count

        # Check the sequences for breaks
        hsequence = sequence[start:stop+1].upper()
        if not len(hsequence.replace("-","")) == len(hsequence): continue
        if not len(hsequence.replace(".","")) == len(hsequence): continue
        csequence = entry["chimp"][start:stop+1].upper()
        if not len(csequence.replace("-","")) == len(csequence): continue
        if not len(csequence.replace(".","")) == len(csequence): continue
        msequence = entry["mouse"][start:stop+1].upper()
        if not len(msequence.replace("-","")) == len(msequence): continue
        if not len(msequence.replace(".","")) == len(msequence): continue
        rsequence = entry["rat"][start:stop+1].upper()
        if not len(rsequence.replace("-","")) == len(rsequence): continue
        if not len(rsequence.replace(".","")) == len(rsequence): continue
        dsequence = entry["dog"][start:stop+1].upper()
        if not len(dsequence.replace("-","")) == len(dsequence): continue
        if not len(dsequence.replace(".","")) == len(dsequence): continue

        # Translate sequences according to strand identity
        if strand == "-":
            hsequence = translate(hsequence)
            csequence = translate(csequence)
            msequence = translate(msequence)
            rsequence = translate(rsequence)
            dsequence = translate(dsequence)

        # If there are no breaks, create the motif
        ret.append([ hsequence, csequence, msequence, rsequence, dsequence ])
        
    return ret

# Write a motif to the Estimator input file
def writeMotif(motif, f):
    for word in motif:
        f.write(word+"\n")
    f.write("\n-2 -1 -1 -1 -1\n-1 -2 -1 -1 -1\n-1 -1 "+\
            "-2 -1 -1\n-1 -1 -1 -2 -1\n-1 -1 -1 -1 -2\n")

# Unbroken? No gaps?
def unbroken(sequence):
    return (len(sequence) == len(sequence.replace("_","")) and
            len(sequence) == len(sequence.replace(".","")))

# Grab a random motif from the alignment
def mkRandomMotif(alignment, length):
    entry = alignment[random.sample(alignment.keys(), 1)[0]]
    for i in range(len(entry["human"])):
        hs = entry["human"][i:i+length]
        cs = entry["chimp"][i:i+length]
        ms = entry["mouse"][i:i+length]
        rs = entry["rat"][i:i+length]
        ds = entry["dog"][i:i+length]
        if unbroken(hs) and unbroken(cs) and unbroken(ms) \
               and unbroken(rs) and unbroken(ds):
            return [ hs, cs, ms, rs, ds ]
    return []

# Load a tree structure from file
def loadTree(filename):
    
    # Parse the Estimator output
    f = open(filename, "r")
    tree_output = f.readlines()
    treelines = []
    treeline = ""
    for i in range(len(tree_output)):
        if tree_output[i] == "\n":
            treelines.append(treeline)
            treeline = ""
        else: treeline += tree_output[i]
    f.close()

    # Use only the first tree structure, since all are equally well conserved
    return Tree(treelines[0])

# AW YEAH BABY (fitchslap)
def fitchslap(node):

    #print "ITERATION:", node

    # Grab the costs matrix
    global costs
    
    # If we're at a leaf, just return
    if len(node.children) == 0: return

    # Take the union of the child sets
    candidates = []
    for c in node.children:
        fitchslap(c)
        for i in range(len(c.color)):
            if not c.color[i] in candidates: candidates.append(c.color[i])

    # Remove all suboptimal candidates
    subcosts = dict()
    smallest = 3000
    for color in candidates:
        cost = 0
        for c in node.children: cost += costs[color+c.color[0]]
        subcosts[color] = cost
        if cost < smallest: smallest = cost
    for color in subcosts.keys():
        if smallest < subcosts[color]:
            del candidates[candidates.index(color)]

    # Add the viable candidates to the color
    node.color = ""
    for c in candidates: node.color += c
    #print "\t", node.color, node.children[0].index, node.children[1].index
    return

# Feeyotch, beeyotch
def feeyotch(node):
    if len(node.children) == 0: return
    node.color = node.color[0]
    for c in node.children: feeyotch(c)
    return

# Color a new tree given a tree and some submotif
def color(tree, submotif):

    #print submotif
    
    # Create a copy of the tree
    ret = tree.copy(tree.r)
    
    # Color the taxa
    for i in range(len(submotif)):
        for j in range(len(ret.nodes)):
            node = ret.nodes[j]
            if str(i+1) == node.index:
                node.color = submotif[i]
                
    # Color the ancestors (Fitch)
    fitchslap(ret.r)
    feeyotch(ret.r)
    return ret

# Sum the CMs across two adjacent position trees
def countCMs(t1, t2, cmuts):
    if len(t1.children) == 0: return cmuts
    left_key = t1.color+"->"+t1.children[0].color+\
               "/"+t2.color+"->"+t2.children[0].color
    cmuts[left_key] += 1
    right_key = t1.color+"->"+t1.children[1].color+\
                "/"+t2.color+"->"+t2.children[1].color
    cmuts[right_key] += 1
    cmuts = countCMs(t1.children[0], t2.children[0], cmuts)
    cmuts = countCMs(t1.children[1], t2.children[1], cmuts)
    return cmuts

# Main function
def main(args):

    # Parse arguments for PWM filename
    if len(args) != 2:
        print "Usage: ./cm_scan.py PWM_FILE"
        sys.exit(1)
    pwm_file = args[1]

    # Hard-coded values
    p = -12.0
    fasta = "../doc/fasta/human1ku-corrected.fasta"
    a_file = open("../doc/pickle/alignment.pickle", "r")

    # Load the 5-way alignment
    alignment = pickle.load(a_file)

    # Close relevant files
    a_file.close()

    # Run PWM_SCAN on the human FASTA, generate a GFF file
    pwm = pwm_file.replace(".pwm","")
    #os.popen("pwm_scan -f "+fasta+" -M "+pwm+" -p "+str(p)+" -s 1")
    #os.popen("pwm_scan -f "+fasta+" -M "+pwm+" -p "+str(p)+" -s 2")
    gff = open(fasta+"."+pwm.replace("../doc/pwm/","")+".gff", "r")

    # For each hit in the GFF, throw out hits that span gaps, build motifs
    motifs = buildMotifs(gff, alignment)
    if len(motifs) == 0:
        print "NO HITS"
        sys.exit(0)
    #print motifs
    #pickle.dump(hits, open("hits.tmp", "w"))
    #motifs = pickle.load(open("hits.tmp", "r"))

    # Create the CM mega-dictionary
    mega = []
    for i in range(len(motifs[0][0])-1):
        cmuts = dict()
        alphabet = [ "A->A", "A->C", "A->G", "A->T",
                     "C->A", "C->C", "C->G", "C->T",
                     "G->A", "G->C", "G->G", "G->T",
                     "T->A", "T->C", "T->G", "T->T" ]
        for i in alphabet:
            for j in alphabet:
                cmuts[i+"/"+j] = 0
        mega.append(cmuts)

    # Create the tree structure
    supertree = Tree([])
    root = TreeNode("r")
    rodent = TreeNode("")
    primate = TreeNode("")
    highermammal = TreeNode("")
    mouse = TreeNode("3")
    rat = TreeNode("4")
    human = TreeNode("1")
    chimp = TreeNode("2")
    dog = TreeNode("5")
    root.children = [ rodent, highermammal ]
    rodent.children = [ mouse, rat ]
    highermammal.children = [ dog, primate ]
    primate.children = [ human, chimp ]
    supertree.r = root
    supertree.node = [ root, highermammal, rodent, mouse,
                       rat, primate, human, chimp ]

    # For each motif:
    for m in motifs:

        # Check the motif: is it any good?
        flag = True
        for subm in m:
            if len(subm) != len(m[0]): flag = False
            if "N" in subm: flag = False
        if not flag: continue
        
        # Find the aligned sequence from other species, write motif to file
        f = open("cm_scan.tmp", "w")
        writeMotif(m, f)
        f.close()
        
        # Run the Ptree estimator to get a tree structure
        #os.popen("./estimator.sh")

        # Load the tree structure from the Estimator output
        #supertree = loadTree("cm_scan.tree")

        # For each position in the hit, color a new tree
        trees = []
        for position in range(len(m[0])):
            submotif = []
            for i in range(len(m)): submotif.append(m[i][position])
            trees.append(color(supertree, submotif))
            
        # For every pair of adjacent trees, count CMs
        for i in range(len(m[0])-1):
            t1 = trees[i]
            t2 = trees[i+1]
            countCMs(t1.r, t2.r, mega[i])

    # Do a chi-squared analysis on each position matrix
    print "DEGREES OF FREEDOM:", 15*15
    alphabet = [ "A->A", "A->C", "A->G", "A->T",
                 "C->A", "C->C", "C->G", "C->T",
                 "G->A", "G->C", "G->G", "G->T",
                 "T->A", "T->C", "T->G", "T->T" ]
    #print mega[0]["A->A/G->G"]
    for i in range(len(motifs[0][0])-1):
        m = zeros((len(alphabet), len(alphabet)))
        for j in range(len(alphabet)):
            for k in range(len(alphabet)):
                #print i, j, k
                #print "\t", alphabet[j], alphabet[k]
                m[j,k] = mega[i][alphabet[j]+"/"+alphabet[k]]
        print "PVALUE at pos (", i, ", ", i+1, ")", chisquarematrix.chisquare(m)

    #for k in cmuts.keys():
    #    print k, cmuts[k]

    # Sample random aligned sequence 10000 times
    #print "Starting random"
    #cmuts_random = dict()
    #for i in alphabet:
    #    for j in alphabet:
    #        cmuts_random[i+"/"+j] = 0
    #count = 0
    #rmotifs = []
    #while count < 10000:
    #    motif = mkRandomMotif(alignment, len(motifs[0][0]))
    #    if len(motif) == 0: continue
    #    rmotifs.append(motif)
    #    count += 1
    #    print count

    #pickle.dump(rmotifs, open("rmotifs.pickle", "w"))
    #rmotifs = pickle.load(open("rmotifs.pickle", "w"))
    
    # Perform the above technique on that sequence, treating it as a hit
    #for m in rmotifs:
        
        # For each position in the hit, color a new tree
    #    trees = []
    #    for position in range(len(m[0])):
    #        submotif = []
    #        for i in range(len(m)): submotif.append(m[i][position])
    #        trees.append(color(supertree, submotif))

        # For every pair of adjacent trees, count CMs
    #    for i in range(len(m[0])-1):
    #        t1 = trees[i]
    #        t2 = trees[i+1]
    #        countCMs(t1.r, t2.r, cmuts_random)

    # Do matrix statistics
    #for k in cmuts_random.keys():
    #    experimental_mean = cmuts[k]/len(motifs)
    #    sample_mean = cmuts_random[k]/len(rmotifs)

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
