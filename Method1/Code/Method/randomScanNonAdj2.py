#!/usr/bin/python
# Greg Donahue, 09/23/2005
# -----------------------------------------------------------------------------
import sys, re, pickle, os, random, copy, math
from Numeric import array, zeros
from rpy import *
# -----------------------------------------------------------------------------
# GLOBALS
# A list representing our nucleotide alphabet
nucleotides = [ "A", "T", "G", "C" ]
# A cost table for doing the Fitch algorithm (dictionary representation)
# Vastly simplified, it assigns -1 to every X->Y edge (X!=Y) and -2 to every
# Y->Y edge (Y==Y)
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

        # This is a horrible kludge
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

    def printIt(self):
        print "\t   ", self.r.color
        print "     ", self.r.children[0].color, "\t", self.r.children[1].color
        print "\t\t    ", self.r.children[1].children[1].color
        print "  ", self.r.children[0].children[0].color, " ", self.r.children[0].children[1].color, \
        "   ", self.r.children[1].children[0].color, "   ", self.r.children[1].children[1].children[0].color, " ", \
        self.r.children[1].children[1].children[1].color
        print "\n\n\n"

    def fullCopy(self, ret):        
        ret.r.color = self.r.color        
        ret.r.children[0].color = self.r.children[0].color
        ret.r.children[1].color = self.r.children[1].color
        ret.r.children[1].children[1].color = self.r.children[1].children[1].color
        ret.r.children[0].children[0].color = self.r.children[0].children[0].color
        ret.r.children[0].children[1].color = self.r.children[0].children[1].color
        ret.r.children[1].children[0].color = self.r.children[1].children[0].color
        ret.r.children[1].children[1].children[0].color = self.r.children[1].children[1].children[0].color
        ret.r.children[1].children[1].children[1].color = self.r.children[1].children[1].children[1].color

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

        # Collect the sequence elements from the hits
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
        if "N" not in hsequence and "N" not in csequence and "N" not in msequence and "N" not in rsequence and "N" not in dsequence:        
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

def getSamples(i, skip, pval):
    try:
        if i+1 < 10:
            f = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/MA000" + str(i+1) + "(p " + str(pval) + ").gff.ll")
        elif i+1 < 100:
            f = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/MA00" + str(i+1) + "(p " + str(pval) + ").gff.ll")
        else:
            f = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/MA00" + str(i+1) + "(p " + str(pval) + ").gff.ll")
        s = f.readline()
        while s.split()[0] != "Hit": s = f.readline()
        hits = int(s.split()[2])
        f.close
        return hits               
    except IOError:
        print "Not here", "MA000" + str(i+1) + "(p " + str(pval) + ").gff.ll"
        return 0

def getPositions(i, skip, pval):
    try:
        if i+1 < 10:
            f = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/MA000" + str(i+1) + "(p " + str(pval) + ").gff.ll")
        elif i+1 < 100:
            f = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/MA00" + str(i+1) + "(p " + str(pval) + ").gff.ll")
        else:
            f = open("../../Results/Real/Raw Data/Gap " + str(skip) + "/MA00" + str(i+1) + "(p " + str(pval) + ").gff.ll")
        s = f.readline()
        length = 0        
        while s != "":
            s = f.readline()
            length += 1
        f.close()
        #print hits * (length - skip - 1)
        print "len calc", length, skip+1
        return (length - 1)        
    except IOError:
        print "Not here", "MA000" + str(i+1) + "(p " + str(pval) + ").gff.ll"
        return 0

def doMega(mega, logs):
    exp = float((mega[0]["CN"]+mega[0]["CC"])*(mega[0]["NC"]+mega[0]["CC"]))/(mega[0]["CN"]+mega[0]["NC"]+mega[0]["CC"]+mega[0]["NN"])
    if mega[0]["CC"] != 0:
        logs.append(math.log(float(mega[0]["CC"])/exp))
    if mega[0]["CC"] == 0 and exp == float(0):
        logs.append(0)
    else: print "no log"

def getExp(mega):
    return float((mega[0]["CN"]+mega[0]["CC"])*(mega[0]["NC"]+mega[0]["CC"]))/(mega[0]["CN"]+mega[0]["NC"]+mega[0]["CC"]+mega[0]["NN"])

def make16Sq():
    r = dict()
    nucs = ["A", "C", "G", "T"]
    edges = []
    for n1 in nucs:
        for n2 in nucs:
            r[n1+n2] = float(0)    
    return r
    
#########################################
################################################
def mainRunner(skip, pval, treeList):
    #pwm_file = argPWM
    #pwm_file = "MA0005(p -9).gff"
##
##    # Hard-coded values
##    p = -12.0
##    fasta = "../doc/fasta/human1ku-corrected.fasta"
    a_file = open("../../Doc/pickle/alignment.pickle", "r")
##
##    # Load the 5-way alignment
    alignment = pickle.load(a_file)
##
##    # Close relevant files
    a_file.close()
##
##    # Run PWM_SCAN on the human FASTA, generate a GFF file
##    # (Actually, the program assumes that you have already done this separately
##    # and just goes ahead and tries to find the GFF file)
    #pwm = pwm_file.replace(".pwm","")
##    #gff = open(fasta+"."+pwm.replace("../doc/pwm/","")+".gff", "r")
    #gff = open("../doc/fasta/" + pwm); 
##
##    # For each hit in the GFF, throw out hits that span gaps, build motifs
    #motifs = buildMotifs(gff, alignment)
    #if len(motifs) == 0:
    #    print "NO HITS"
    #    return
##    # used in analysis of distributions PE
##    # can be commented out if needed
    #random.shuffle(motifs)
    #if len(motifs) > 1000: motifs = motifs[0:1000]
##
    # Create the tree structure
##    supertree = Tree([])
##    root = TreeNode("r")
##    rodent = TreeNode("")
##    primate = TreeNode("")
##    highermammal = TreeNode("")
##    mouse = TreeNode("3")
##    rat = TreeNode("4")
##    human = TreeNode("1")
##    chimp = TreeNode("2")
##    dog = TreeNode("5")
##    root.children = [ rodent, highermammal ]
##    rodent.children = [ mouse, rat ]
##    highermammal.children = [ dog, primate ]
##    primate.children = [ human, chimp ]
##    supertree.r = root
##    supertree.node = [ root, highermammal, rodent, mouse,
##                       rat, primate, human, chimp ]

    #In true genomic fashion, here is a large region of junk
    #I believe that it was originally used to create the random tree-pair
    #distribution that our distribution was superimposed on (the red curve in
    #the graph) ~Greg 01/04/06
    #Construct random pairs of positions
    count = 0
    #only keep one log for mega
    incrementor, logs = 0, []    
    instances = []
    #for each PWM
    
    while incrementor < 111:        
        #This keeps track of all 16 nt combination counts
        #between two positions, for each taxon
        #The order in the submotif is human, chimp, mouse, rat, dog
        #each 16Sq is indexed by pos, pos+1
        labels = ["human", "chimp", "mouse", "rat", "dog"]
        big16 = dict()
        for label in labels:
            big16[label] = []  
            
        print incrementor
        samples = 0        
        #look at each position
        positionMax = getPositions(incrementor, skip, pval)
        pos = 0
        #look at all the hits for this PWM (incrementor)
        samples = getSamples(incrementor, skip, pval)
        #samples = 1000
        if samples == 0:
            print "no samples"
            incrementor += 1
            continue
        else: print "samples", samples
        print "posmax", positionMax
        while pos < positionMax:
            i = 0
            instances = []
            while i < samples:
                # Pick a gene from the alignment
                gene_key = random.sample(alignment, 1)[0]
                sequence = alignment[gene_key]["human"]
                # Pick a working 2-mer motif from the gene
                check = 0
                flag = False
                while True:
                    position = int(random.uniform(0, len(sequence)-2-skip))
                    try:
                        if verify(alignment, gene_key, position, skip):
                            break
                    except: print "index error in align"
                    if check > 10:
                        flag = True
                        break
                    else: check += 1
                if flag: continue
                # Build the instance
                instance = [ [alignment[gene_key]["human"][position].upper(), alignment[gene_key]["human"][position+1+skip].upper()],
                             [alignment[gene_key]["chimp"][position].upper(), alignment[gene_key]["chimp"][position+1+skip].upper()],
                             [alignment[gene_key]["mouse"][position].upper(), alignment[gene_key]["mouse"][position+1+skip].upper()],
                             [alignment[gene_key]["rat"][position].upper(), alignment[gene_key]["rat"][position+1+skip].upper()],
                             [alignment[gene_key]["dog"][position].upper(), alignment[gene_key]["dog"][position+1+skip].upper()] ]
            
                # Add this instance to instances
                instances.append(instance)
                i += 1

            #matricies go here
            for label in labels:
                big16[label].append(make16Sq())
            #big16.append(make16Sq())
            # Create the change/no-change matrix
            mega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet:
                        cmuts[i+j] = 0
                cmuts["CC"] = float(41);
                cmuts["CN"] = float(387);
                cmuts["NC"] = float(404);
                cmuts["NN"] = float(5568);
                mega.append(cmuts)

            #root rodent
            rootRodentMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(4.5);
                cmuts["CN"] = float(65.5);
                cmuts["NC"] = float(60.5);
                cmuts["NN"] = float(669.5);
                rootRodentMega.append(cmuts)

            rootHMMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(4.5);
                cmuts["CN"] = float(65.5);
                cmuts["NC"] = float(60.5);
                cmuts["NN"] = float(669.5);
                rootHMMega.append(cmuts)

            rodentRatMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(1.5);
                cmuts["CN"] = float(54);
                cmuts["NC"] = float(56);
                cmuts["NN"] = float(688.5);
                rodentRatMega.append(cmuts)

            rodentMouseMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(1.5);
                cmuts["CN"] = float(46);
                cmuts["NC"] = float(57);
                cmuts["NN"] = float(695.5);
                rodentMouseMega.append(cmuts)

            HMDogMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(23.75);
                cmuts["CN"] = float(98.75);
                cmuts["NC"] = float(84.75);
                cmuts["NN"] = float(592.75);
                HMDogMega.append(cmuts)

            HMPrimateMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(5.25);
                cmuts["CN"] = float(50.25);
                cmuts["NC"] = float(75.25);
                cmuts["NN"] = float(669.25);
                HMPrimateMega.append(cmuts)

            primateHumanMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(0);
                cmuts["CN"] = float(2);
                cmuts["NC"] = float(3.5);
                cmuts["NN"] = float(794.5);
                primateHumanMega.append(cmuts)

            primateChimpMega = []
            for i in range(len(instances[0][0])-1):
                cmuts = dict()
                alphabet = [ "N", "C" ]
                for i in alphabet:
                    for j in alphabet: cmuts[i+j] = 0
                cmuts["CC"] = float(0);
                cmuts["CN"] = float(5);
                cmuts["NC"] = float(6.5);
                cmuts["NN"] = float(788.5);
                primateChimpMega.append(cmuts)

            # Count edges in the instance set
            for m in instances:

                #update big 16
                for u in range(len(labels)):
                    big16[labels[u]][pos][m[u][0]+m[u][1]] += 1
                # Build a set of trees out of this motif
                trees = []
                for position in range(len(m[0])): #0 and 1 I think
                    submotif = ""
                    for i in range(len(m)): submotif += m[i][position]
                    trees.append(treeList[submotif])

                # Walk down pairs of trees, fill in C/NC matrix
                for position in range(len(m[0])-1):
                    t1_set = trees[position]
                    t2_set = trees[position+1]

                    weight = float(1)/float(len(t1_set)*len(t2_set))
                    if len(t1_set) == 1 and len(t2_set) == 1: weight = 1

                    for t1 in t1_set:
                        for t2 in t2_set:                    
                            # Decompose the tree at position i 
                            t1_r = t1.r
                            t1_rodent = t1_r.children[0]
                            t1_highermammal = t1_r.children[1]
                            t1_mouse = t1_rodent.children[0]
                            t1_rat = t1_rodent.children[1]
                            t1_dog = t1_highermammal.children[0]
                            t1_primate = t1_highermammal.children[1]
                            t1_human = t1_primate.children[0]
                            t1_chimp = t1_primate.children[1]
                                
                            # Decompose the tree at position i+2 
                            t2_r = t2.r
                            t2_rodent = t2_r.children[0]
                            t2_highermammal = t2_r.children[1]
                            t2_mouse = t2_rodent.children[0]
                            t2_rat = t2_rodent.children[1]
                            t2_dog = t2_highermammal.children[0]
                            t2_primate = t2_highermammal.children[1]
                            t2_human = t2_primate.children[0]
                            t2_chimp = t2_primate.children[1]
                            
                            # Compare edges and update the C/NC matrix
                            resolve(position, mega, rootRodentMega, t1_r, t2_r, t1_rodent, t2_rodent, weight)
                            resolve(position, mega, rodentMouseMega, t1_rodent, t2_rodent, t1_mouse, t2_mouse, weight)
                            three = resolve(position, mega, rodentRatMega, t1_rodent, t2_rodent, t1_rat, t2_rat, weight)
                            resolve(position, mega, rootHMMega, t1_r, t2_r, \
                                    t1_highermammal, t2_highermammal, weight)
                            resolve(position, mega, HMPrimateMega, t1_highermammal, t2_highermammal, \
                                    t1_primate, t2_primate, weight)
                            resolve(position, mega, HMDogMega, t1_highermammal, t2_highermammal, \
                                    t1_dog, t2_dog, weight)
                            resolve(position, mega, primateHumanMega, t1_primate, t2_primate, t1_human, t2_human, weight)
                            resolve(position, mega, primateChimpMega, t1_primate, t2_primate, t1_chimp, t2_chimp, weight)

            #exp = float((mega[0]["CN"]+mega[0]["CC"])*(mega[0]["NC"]+mega[0]["CC"]))/(mega[0]["CN"]+mega[0]["NC"]+mega[0]["CC"]+mega[0]["NN"])
            exp = float(0)
            exp += getExp(rootRodentMega)
            exp += getExp(rootHMMega)
            exp += getExp(rodentRatMega)
            exp += getExp(rodentMouseMega)
            exp += getExp(HMDogMega)
            exp += getExp(HMPrimateMega)
            exp += getExp(primateChimpMega)
            exp += getExp(primateHumanMega)
            #if mega[0]["CC"] != 0:
            logs.append(math.log(float(mega[0]["CC"])/exp))
            count += 1
            print count            
            #go to the next "position" for this PWM
            pos += 1
##        newFile = open("NA " + str(skip) + " p " + str(pval) + " " + str(incrementor+1) + ".16sq", "w")
##        pickle.dump(big16, newFile)
##        newFile.close()
        #go to the next PWM
        incrementor += 1

    f = open("../../Results/Controls/Same Site/Raw Data/same site Gap " + str(skip) + ".txt", "w")
    pickle.dump(logs, f)
    f.close()
    
    #for i in range(len(logs)):
    #    print i, logs[i]
    print "logs length", len(logs)
    #sys.exit(0)
        
###########################################################
######################################################

# Main function
def main(args):
    sys.exit(0)

# Resolve edge dependencies into one of four classes
def resolve(position, matrix, matrix2, anc1, anc2, dec1, dec2, weight):
    flag_i, flag_j = "N", "N"
    if anc1.color != dec1.color: flag_i = "C"
    if anc2.color != dec2.color: flag_j = "C"
    matrix[position][flag_i+flag_j] += weight
    matrix2[position][flag_i+flag_j] += weight    

# Verify that the alignment at this position works
def verify(alignment, key, position, skip):
    h = alignment[key]["human"]
    c = alignment[key]["chimp"]
    m = alignment[key]["mouse"]
    r = alignment[key]["rat"]
    d = alignment[key]["dog"]
    a = [ "A", "T", "G", "C" ]
    if h[position].upper() in a and h[position+1+skip].upper() in a and \
       c[position].upper() in a and c[position+1+skip].upper() in a and \
       m[position].upper() in a and m[position+1+skip].upper() in a and \
       r[position].upper() in a and r[position+1+skip].upper() in a and \
       d[position].upper() in a and d[position+1+skip].upper() in a:
        return True
    else: return False

# Print a tree recursively
def treeprint(treenode):
    print treenode.color
    if len(treenode.children) == 0: return
    treeprint(treenode.children[0])
    treeprint(treenode.children[1])

# Validate the motif
def validate(motif):
    length = len(motif[0])
    for i in motif[1:]:
        if length != len(i): return False
    return True

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
