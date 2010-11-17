#!/usr/bin/python
# Greg Donahue, 09/23/2005
# -----------------------------------------------------------------------------
import sys, re, pickle, os, random, chisquarematrix, copy, math
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
        length =  len(alignment.keys())
        r = random.randint(0, length-1)        
        rkey = alignment.keys()[r]
        entry = alignment[rkey]
        flag = False
             
        #entry = alignment[id]
        
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

def intersection(c1, c2):
    intersect = ""
    union = c1
    for e in c2:
        if e not in c1: union+=e
        if e in c1: intersect+=e
    return [intersect, union]

# fitchslap is the Fitch forward pass
# A simpler, cleaner version of this and the backward pass (feeyotch) can be
# found in scanner.py (uniformly-weighted parsimony)
def fitchslap(node):

    #chimp
    #human
    i = intersection(node.children[1].children[1].children[0].color, node.children[1].children[1].children[1].color)
    if i[0] != "": node.children[1].children[1].color = i[0]
    else: node.children[1].children[1].color = i[1]
    #print "primate", node.children[1].children[1].color

    #dog
    #primate
    i = intersection(node.children[1].children[1].color, node.children[1].children[0].color)
    if i[0] != "": node.children[1].color = i[0]
    else: node.children[1].color = i[1]

    #rat
    #mouse
    i = intersection(node.children[0].children[1].color, node.children[0].children[0].color)
    if i[0] != "": node.children[0].color = i[0]
    else: node.children[0].color = i[1]

    i = intersection(node.children[0].color, node.children[1].color)
    if i[0] != "": node.color = i[0]
    else: node.color = i[1]   

    return

def myHash(index, tree, color):    
    if index == 10:
        tree.r.children[0].color = color
    elif index == 11:
        tree.r.children[0].children[0].color = color
    elif index == 12:
        tree.r.children[0].children[1].color = color
    elif index == 1:
        tree.r.children[1].color = color
    elif index == 2:
        tree.r.children[1].children[0].color = color
    elif index == 3:        
        tree.r.children[1].children[1].color = color
    elif index == 4:
        tree.r.children[1].children[1].children[0].color = color
    else:
        tree.r.children[1].children[1].children[1].color = color

def score(tree):
    s = 0
    if tree.r.color != tree.r.children[0].color: s += 1
    if tree.r.color != tree.r.children[1].color: s += 1
    if tree.r.children[0].color != tree.r.children[0].children[0].color: s += 1
    if tree.r.children[0].color != tree.r.children[0].children[1].color: s += 1
    if tree.r.children[1].color != tree.r.children[1].children[0].color: s += 1
    if tree.r.children[1].color != tree.r.children[1].children[1].color: s += 1
    if tree.r.children[1].children[1].color != tree.r.children[1].children[1].children[0].color: s += 1
    if tree.r.children[1].children[1].color != tree.r.children[1].children[1].children[1].color: s += 1
    return s
    
def scoreIt(trees):
    #print "SCORES HERE++++++++++++++++++++++++++++"
    scores = []
    for t in trees:
        scores.append(score(t))
        #t.printIt()
    #print "scores", scores
    #print "END SCORES___________________________"    
    m = min(scores)
    ret = []
    for i in range(len(trees)):
        if scores[i] == m: ret.append(trees[i])
    return ret

# Color a new tree given a tree and some submotif
def color(tree, submotif):
    trees = []   
    # Create a copy of the tree
    ret = tree.copy(tree.r)
    #add this copy to the tree list
       
    # Color the taxa
    for i in range(len(submotif)):
        for j in range(len(ret.nodes)):
            node = ret.nodes[j]
            if str(i+1) == node.index:
                node.color = submotif[i]
    
    # Color the ancestors (Fitch)    
    fitchslap(ret.r)
    #print "The colored tree:"
    #ret.printIt()
    #root
    node = ret.r.color
    #print "root", ret.r.color
    index = 0
    
    for i in node: #copy in a new tree and set the root        
        trees.append(tree.copy(tree.r))
        ret.fullCopy(trees[index])
        trees[index].r.color = i        
        #trees[index].printIt()
        index += 1

    #look at level 1, rodent and higher mammal
    #grab any root's children b/c at this point they are all the same
    rodentColor = trees[0].r.children[0].color
    hmColor = trees[0].r.children[1].color
    goodSetRodent = [] #the list of intersections with the possible roots
    goodSetHM = []
    roots = []
    interCounter = 0
    for t in trees:
        goodSetRodent.append(intersection(t.r.color, rodentColor))
        if len(goodSetRodent[interCounter][0]) > 0:            
            pass #that's good, keep this set b/c there was an intersection
        else: #no intersection, so just use the colors at this node and pick randomly
            goodSetRodent[interCounter][0] = t.r.children[0].color
        goodSetHM.append(intersection(t.r.color, hmColor))
        if len(goodSetHM[interCounter][0]) > 0:
            pass #that's good, keep this set b/c there was an intersection
        else: #no intersection, so just use the colors at this node and pick randomly
            goodSetHM[interCounter][0] = t.r.children[1].color
        #print "initlal check", goodSetHM[interCounter][0]
        interCounter += 1

    #now you have sets of possible colorings in rodent and hm for each of the given roots
    #need to find out how many more new trees you need to add
    newTreeCount = 0
    startTreeLength = len(trees)
    #print "look at the set", len(goodSetHM[0][0]), goodSetHM[0][0]
    for i in range(startTreeLength):
        if 1 == len(goodSetRodent[i][0])*len(goodSetHM[i][0]): #you only have one tree, is it enough?
            #print "only need 1"
            trees[i].r.children[0].color = goodSetRodent[i][0]
            trees[i].r.children[1].color = goodSetHM[i][0]
            #print trees[i].r.children[1].color
        else: #need more trees
            #print "need more trees"
            currentTreeLength = len(trees)
            #start by assigning the one already in trees
            #here's the problem
            #you need to choose one from rodent and one from hm,
            #but you still need all matchings
##            trees[i].r.children[0].color = goodSetRodent[i][0]
##            trees[i].r.children[1].color = goodSetHM[i][0]
##            goodSetRodent[i][0] = goodSetRodent[i][0][1:len(goodSetRodent[i][0])]
##            goodSetHM[i][0] = goodSetHM[i][0][1:len(goodSetHM[i][0])]
            moreIndex = -1            
            #print "check length again", len(goodSetHM[i][0]), len(goodSetRodent[i][0])
            for r in range(len(goodSetRodent[i][0])):
                for hm in range(len(goodSetHM[i][0])):
                    if moreIndex == -1:
                        #adding to a tree in the list
                        #print "(", r, hm, ")"
                        trees[i].r.children[0].color = goodSetRodent[i][0][r]
                        trees[i].r.children[1].color = goodSetHM[i][0][hm]
                        #print "adding to list tree", goodSetHM[i][hm], goodSetHM[i][0], goodSetHM[i][1]
                    else:
                        trees.append(tree.copy(tree.r))
                        trees[i].fullCopy(trees[currentTreeLength + moreIndex])
                        trees[currentTreeLength + moreIndex].r.children[0].color = goodSetRodent[i][0][r]
                        trees[currentTreeLength + moreIndex].r.children[1].color = goodSetHM[i][0][hm]
                    moreIndex += 1

    #print "check1", trees[0].r.children[1].color
    #look at level 2, just the primate
    dogColor = trees[0].r.children[1].children[1].color
    #print "dogColor", dogColor
    goodSetDog = []
    interCounter = 0
    for t in trees:
        goodSetDog.append(intersection(t.r.children[1].color, dogColor))
        #print goodSetDog
        if len(goodSetDog[interCounter][0]) > 0:
            pass
        else:
            #print "the assigned color", t.r.children[1].children[1].color
            goodSetDog[interCounter][0] = t.r.children[1].children[1].color        
        #print "the color", goodSetDog[interCounter][0]
        interCounter += 1
    #now do the coloring sets
    #print "check2", trees[0].r.children[1].color
    newTreeCount = 0
    newTreeCount = 0
    startTreeLength = len(trees)
    #print "the length", startTreeLength
    for i in range(startTreeLength):
        #print "doing tree", i
        if 1 == len(goodSetDog[i][0]): #you only have one tree, is it enough?
            trees[i].r.children[1].children[1].color = goodSetDog[i][0]
        else: #need more trees
            currentTreeLength = len(trees)
            #start by assigning the one already in trees
            #print "doglength", len(goodSetDog[i][0])
            trees[i].r.children[1].children[1].color = goodSetDog[i][0][0]            
            goodSetDog[i][0] = goodSetDog[i][0][1:len(goodSetDog[i][0])]
            moreIndex = 0
            for r in range(len(goodSetDog[i][0])):                                   
                trees.append(tree.copy(tree.r))
                trees[i].fullCopy(trees[currentTreeLength + moreIndex])                
                trees[currentTreeLength + moreIndex].r.children[1].children[1].color = goodSetDog[i][0][r]                
                moreIndex += 1
     
    #for t in trees: t.printIt()
    trees = scoreIt(trees)
    if len(trees) > 1:
        ri = random.randint(0, len(trees)-1)
        return trees[ri]
    else: return trees[0]       

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

#########################################
################################################
def mainRunner(argPWM, skip):
    pwm_file = argPWM
    #pwm_file = "MA0005(p -9).gff"
##
##    # Hard-coded values
##    p = -12.0
##    fasta = "../doc/fasta/human1ku-corrected.fasta"
    a_file = open("../doc/pickle/alignment.pickle", "r")
##
##    # Load the 5-way alignment
    alignment = pickle.load(a_file)
##
##    # Close relevant files
    a_file.close()
##    keys = alignment.keys
##    for k in keys:
##        print alignment[k]["human"]
##        sys.exit(0)
##    print "size", len(alignment)
##    sys.exit(0)
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
        #print "NO HITS"
        #return
##    # used in analysis of distributions PE
##    # can be commented out if needed
    #random.shuffle(motifs)
    #if len(motifs) > 1000: motifs = motifs[0:1000]
##
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

    #In true genomic fashion, here is a large region of junk
    #I believe that it was originally used to create the random tree-pair
    #distribution that our distribution was superimposed on (the red curve in
    #the graph) ~Greg 01/04/06
    #Construct random pairs of positions

    incrementor, logs = 0, []
    changes = dict()
    changes["root rodent"] = 0
    changes["rodent mouse"] = 0
    changes["rodent rat"] = 0
    changes["root hm"] = 0
    changes["hm dog"] = 0
    changes["hm primate"] = 0
    changes["primate human"] = 0
    changes["chimp primate"] = 0    
    
    while incrementor < 10: # number of PWMs
        instances = []
        print incrementor
        i = 0
        while i < 100: # len(motifs): #this is the number of hits for a pwm
            # Pick a gene from the alignment
            gene_key = random.sample(alignment, 1)[0]
            gene_key2 = random.sample(alignment, 1)[0]
            #print i, gene_key, gene_key2
            sequence = alignment[gene_key]["human"]
            sequence2 = alignment[gene_key2]["human"]
            # Pick a working 2-mer motif from the gene
            check = 0
            flag = False
            while True:
                position = int(random.uniform(0, len(sequence)-2))
                position2 = int(random.uniform(0, len(sequence2)-2))
                try:
                    if verify(alignment, gene_key, position, skip) and verify(alignment, gene_key2, position2, skip):
                        break
                except: print "index error in align"
                if check > 10:
                    flag = True
                    break
                else: check += 1
            if flag: continue
            # Build the instance
            instance = [ [alignment[gene_key]["human"][position].upper(), alignment[gene_key2]["human"][position2].upper()],
                         [alignment[gene_key]["chimp"][position].upper(), alignment[gene_key2]["chimp"][position2].upper()],
                         [alignment[gene_key]["mouse"][position].upper(), alignment[gene_key2]["mouse"][position2].upper()],
                         [alignment[gene_key]["rat"][position].upper(), alignment[gene_key2]["rat"][position2].upper()],
                         [alignment[gene_key]["dog"][position].upper(), alignment[gene_key2]["dog"][position2].upper()] ]
            
            # Add this instance to instances
            #print "instance"
            #print instance
            instances.append(instance)
            i += 1

        # Create the change/no-change matrix
        mega = []
        for i in range(len(instances[0][0])-1):
            cmuts = dict()
            alphabet = [ "N", "C" ]
            for i in alphabet:
                for j in alphabet: cmuts[i+j] = 0
            mega.append(cmuts)

        # Count edges in the instance set
        for m in instances:
            
            # Build a set of trees out of this motif
            trees = []
            for position in range(len(m[0])):                
                submotif = []
                for i in range(len(m)): submotif.append(m[i][position])
                trees.append(color(supertree, submotif))
                #trees.append(color(supertree, ["T", "T", "G", "G", "C"]))
            #print "Tree array length", len(trees)
            #print trees[0]
            #print trees[1]
            #sys.exit(0)
            # Walk down pairs of trees, fill in C/NC matrix
            for position in range(len(m[0])-1):
                
                # Decompose the tree at position i
                t1 = trees[position]
                #t1.printIt()
                t1_r = t1.r
                #print "root", t1_r
                t1_rodent = t1_r.children[0]
                #print "rodent", t1_rodent
                t1_highermammal = t1_r.children[1]
                t1_mouse = t1_rodent.children[0]
                t1_rat = t1_rodent.children[1]
                t1_dog = t1_highermammal.children[0]
                t1_primate = t1_highermammal.children[1]
                t1_human = t1_primate.children[0]
                t1_chimp = t1_primate.children[1]
                
                # Decompose the tree at position i+1
                t2 = trees[position+1]
                #t2.printIt()
                t2_r = t2.r
                t2_rodent = t2_r.children[0]
                t2_highermammal = t2_r.children[1]
                t2_mouse = t2_rodent.children[0]
                t2_rat = t2_rodent.children[1]
                t2_dog = t2_highermammal.children[0]
                t2_primate = t2_highermammal.children[1]
                t2_human = t2_primate.children[0]
                t2_chimp = t2_primate.children[1]

                #print "trees compared_______________________________________________"
                
                # Compare edges and update the C/NC matrix
                #print "resolving", t1_r.color, t2_r.color
                one = resolve(position, mega, t1_r, t2_r, t1_rodent, t2_rodent)
                two = resolve(position, mega, t1_rodent, t2_rodent, t1_mouse, t2_mouse)
                three = resolve(position, mega, t1_rodent, t2_rodent, t1_rat, t2_rat)
                four = resolve(position, mega, t1_r, t2_r, \
                        t1_highermammal, t2_highermammal)
                five = resolve(position, mega, t1_highermammal, t2_highermammal, \
                        t1_primate, t2_primate)
                six = resolve(position, mega, t1_highermammal, t2_highermammal, \
                        t1_dog, t2_dog)
                seven = resolve(position, mega, t1_primate, t2_primate, t1_human, t2_human)
                eight = resolve(position, mega, t1_primate, t2_primate, t1_chimp, t2_chimp)
                
                if one or two or three or  four or  five or  six or  seven or eight:
                    print "comparing the trees++++++++++++++++++++++++++++++++++++++++++++"
                    t1.printIt()
                    t2.printIt()
                    print "trees compared_______________________________________________"
        print mega[0]["CC"], mega[0]["CN"]
        print mega[0]["NC"], mega[0]["NN"]
        exp = float(8)*float((mega[0]["CN"]+mega[0]["CC"])*(mega[0]["NC"]+mega[0]["CC"]))/(mega[0]["CN"]+mega[0]["NC"]+mega[0]["CC"]+mega[0]["NN"])
        #exp = round(exp)

##        if mega[0]["CC"] != 0:
##            logs.append(math.log(float(mega[0]["CC"])/exp))
##        elif mega[0]["CC"] == 0 and mega[0]["CN"] == 0 and mega[0]["NC"] == 0:
##            logs.append(0)
        print exp
        if mega[0]["CC"] != 0:
            logs.append(math.log(float(mega[0]["CC"])/exp))
        print "exp", exp
        if mega[0]["CC"] < exp: print "yes"
        
        incrementor += 1
    keys = changes.keys()
    for j in keys:
        print j, changes[j]
    f = open("REALGOOD_random_logs.txt", "w")
    pickle.dump(logs, f)
    f.close()
    sums = 0
    for i in range(len(logs)):
        sums += logs[i]
        print i, logs[i]
    print sums/len(logs)
        
###########################################################
######################################################

# Main function
def main(args):

##    # Parse arguments for PWM filename
##    if len(args) != 2:
##        print "Usage: ./cm_scan.py PWM_FILE"
##        sys.exit(1)
    pwm_file = args[1]
    #pwm_file = "MA0005(p -9).gff"
##
##    # Hard-coded values
##    p = -12.0
##    fasta = "../doc/fasta/human1ku-corrected.fasta"
    a_file = open("../doc/pickle/alignment.pickle", "r")
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
    pwm = pwm_file.replace(".pwm","")
##    #gff = open(fasta+"."+pwm.replace("../doc/pwm/","")+".gff", "r")
    gff = open("../doc/fasta/" + pwm); 
##
##    # For each hit in the GFF, throw out hits that span gaps, build motifs
    motifs = buildMotifs(gff, alignment)
    if len(motifs) == 0:
        print "NO HITS"
        sys.exit(0)
##    # used in analysis of distributions PE
##    # can be commented out if needed
    random.shuffle(motifs)
    if len(motifs) > 1000: motifs = motifs[0:1000]
##
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

    # In true genomic fashion, here is a large region of junk
    # I believe that it was originally used to create the random tree-pair
    # distribution that our distribution was superimposed on (the red curve in
    # the graph) ~Greg 01/04/06
    # Construct random pairs of positions
    
    incrementor, logs = 10, []
    while incrementor < 1000:
        instances = []
        #print incrementor
        i = 0
        while i < len(motifs):
            #print "the count", i
            # Pick a gene from the alignment
            gene_key = random.sample(alignment, 1)[0]
            sequence = alignment[gene_key]["human"]
            # Pick a working 2-mer motif from the gene
            check = 0
            flag = False
            while True:
                position = int(random.uniform(0, len(sequence)-2))
                if verify(alignment, gene_key, position):
                    #print "good align"
                    break
                if check > 10:
                    flag = True
                    break
                else: check += 1
            if flag: continue
            # Build the instance
            instance = [ alignment[gene_key]["human"][position:position+2].upper(),
                         alignment[gene_key]["chimp"][position:position+2].upper(),
                         alignment[gene_key]["mouse"][position:position+2].upper(),
                         alignment[gene_key]["rat"][position:position+2].upper(),
                         alignment[gene_key]["dog"][position:position+2].upper() ]
            
            # Add this instance to instances
            instances.append(instance)
            
            i += 1
            #print "at increment", i
            # Create the change/no-change matrix
        mega = []
        for k in range(len(instances[0][0])-1):
            cmuts = dict()
            alphabet = [ "N", "C" ]
            for t in alphabet:
                for j in alphabet: cmuts[t+j] = 0
            mega.append(cmuts)

        # Count edges in the instance set
        for m in instances:
            
            # Build a set of trees out of this motif
            trees = []
            for position in range(len(m[0])):
                submotif = []
                for k in range(len(m)): submotif.append(m[k][position])
                trees.append(color(supertree, submotif))

            # Walk down pairs of trees, fill in C/NC matrix
            for position in range(len(m[0])-1):

                # Decompose the tree at position i
                t1 = trees[position]
                t1_r = t1.r
                t1_rodent = t1_r.children[0]
                t1_highermammal = t1_r.children[1]
                t1_mouse = t1_rodent.children[0]
                t1_rat = t1_rodent.children[1]
                t1_dog = t1_highermammal.children[0]
                t1_primate = t1_highermammal.children[1]
                t1_human = t1_primate.children[0]
                t1_chimp = t1_primate.children[1]
                
                # Decompose the tree at position i+1
                t2 = trees[position+1]
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
                resolve(position, mega, t1_r, t2_r, t1_rodent, t2_rodent)
                resolve(position, mega, t1_rodent, t2_rodent, t1_mouse, t2_mouse)
                resolve(position, mega, t1_rodent, t2_rodent, t1_rat, t2_rat)
                resolve(position, mega, t1_r, t2_r,
                        t1_highermammal, t2_highermammal)
                resolve(position, mega, t1_highermammal, t2_highermammal,
                        t1_dog, t2_dog)
                resolve(position, mega, t1_highermammal, t2_highermammal,
                        t1_primate, t2_primate)
                resolve(position, mega, t1_primate, t2_primate, t1_human, t2_human)
                resolve(position, mega, t1_primate, t2_primate, t1_chimp, t2_chimp)
    
            exp = float((mega[0]["CN"]+mega[0]["CC"])*(mega[0]["NC"]+mega[0]["CC"]))/(mega[0]["CN"]+mega[0]["NC"]+mega[0]["CC"]+mega[0]["NN"])
            if mega[0]["CC"] == 0: logs.append(1)
            else: logs.append(math.log(float(mega[0]["CC"])/exp))
            
            incrementor += 1

        f = open("random_logs_perryTest." + pwm_file + ".txt", "w")
        pickle.dump(logs, f)
        f.close()
        for i in range(len(logs)):
            print i, logs[i]
        sys.exit(0)

    # Create the change/no-change matrix
    mega = []
    for i in range(len(motifs[0][0])-1):
        cmuts = dict()
        alphabet = [ "N", "C" ]
        for i in alphabet:
            for j in alphabet:
                cmuts[i+j] = 0
        mega.append(cmuts)
    
        
    # For each instance
    for m in motifs:
        
        # Build a set of trees out of this motif
        if not validate(m): continue
        trees = []
        for position in range(len(m[0])):
            submotif = []
            for i in range(len(m)): submotif.append(m[i][position])
            trees.append(color(supertree, submotif))
            
        # Walk down pairs of trees, fill in C/NC matrix
        for position in range(len(m[0])-2):
            
            # Decompose the tree at position i
            t1 = trees[position]
            t1_r = t1.r
            t1_rodent = t1_r.children[0]
            t1_highermammal = t1_r.children[1]
            t1_mouse = t1_rodent.children[0]
            t1_rat = t1_rodent.children[1]
            t1_dog = t1_highermammal.children[0]
            t1_primate = t1_highermammal.children[1]
            t1_human = t1_primate.children[0]
            t1_chimp = t1_primate.children[1]
                
            # Decompose the tree at position i+1
            t2 = trees[position+2]
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
            resolve(position, mega, t1_r, t2_r, t1_rodent, t2_rodent)
            resolve(position, mega, t1_rodent, t2_rodent, t1_mouse, t2_mouse)
            resolve(position, mega, t1_rodent, t2_rodent, t1_rat, t2_rat)
            resolve(position, mega, t1_r, t2_r,
                    t1_highermammal, t2_highermammal)
            resolve(position, mega, t1_highermammal, t2_highermammal,
                    t1_dog, t2_dog)
            resolve(position, mega, t1_highermammal, t2_highermammal,
                    t1_primate, t2_primate)
            resolve(position, mega, t1_primate, t2_primate, t1_human, t2_human)
            resolve(position, mega, t1_primate, t2_primate, t1_chimp, t2_chimp)

    newFile = open("LL for all PWM/" + pwm + ".ll", "w")
    for i in range(len(mega)):
        obs = mega[i]["CC"]
        exp = (mega[i]["NC"]+mega[i]["CC"])*(mega[i]["CN"]+mega[i]["CC"])
        exp = float(exp)/(mega[i]["NC"]+mega[i]["CC"]+mega[i]["CN"]+mega[i]["NN"])
        ll = math.log(float(obs)/exp)
        if ll > 1.5:
            newFile.write("(" + str(i+1) +  ", " + str(i+2) + ") " +  str(ll) + " *\n")
        else:
            newFile.write("(" + str(i+1) +  ", " + str(i+2) + ") " +  str(ll) + " \n")
        #print "(", i+1, ",", i+2, ")", ll, r.pchisq(float(obs)/exp, 2, lower_tail=False)
        #print "(", i+1, ",", i+2, ")", ll, r.pchisq(float(obs)/exp, 2, lower_tail=False)
    newFile.close();
        
    sys.exit(0)
        
    # Create the change/no-change matrix
    r_mega = []
    for i in range(len(r_pairs[0][0])-1):
        cmuts = dict()
        alphabet = [ "N", "C" ]
        for i in alphabet:
            for j in alphabet: cmuts[i+j] = 0
        r_mega.append(cmuts)
    
    # For each instance
    for m in r_pairs:

        # Build a set of trees out of this motif
        trees = []
        for position in range(len(m[0])):
            submotif = []
            for i in range(len(m)): submotif.append(m[i][position])
            trees.append(color(supertree, submotif))

        # Walk down pairs of trees, fill in C/NC matrix
        for position in range(len(m[0])-1):

            # Decompose the tree at position i
            t1 = trees[position]
            t1_r = t1.r
            t1_rodent = t1_r.children[0]
            t1_highermammal = t1_r.children[1]
            t1_mouse = t1_rodent.children[0]
            t1_rat = t1_rodent.children[1]
            t1_dog = t1_highermammal.children[0]
            t1_primate = t1_highermammal.children[1]
            t1_human = t1_primate.children[0]
            t1_chimp = t1_primate.children[1]

            # Decompose the tree at position i+1
            t2 = trees[position+1]
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
            resolve(position, r_mega, t1_r, t2_r, t1_rodent, t2_rodent)
            resolve(position, r_mega, t1_rodent, t2_rodent, t1_mouse, t2_mouse)
            resolve(position, r_mega, t1_rodent, t2_rodent, t1_rat, t2_rat)
            resolve(position, r_mega, t1_r, t2_r,
                    t1_highermammal, t2_highermammal)
            resolve(position, r_mega, t1_highermammal, t2_highermammal,
                    t1_dog, t2_dog)
            resolve(position, r_mega, t1_highermammal, t2_highermammal,
                    t1_primate, t2_primate)
            resolve(position, r_mega, t1_primate, t2_primate, t1_human, t2_human)
            resolve(position, r_mega, t1_primate, t2_primate, t1_chimp, t2_chimp)

    # Calculate log-likelihood scores for each matrix
    print "Random Distribution of Log-likelihoods"
    for i in range(len(r_mega)):
        r_obs = r_mega[i]["CC"]
        r_exp = (r_mega[i]["NC"]+r_mega[i]["CC"])*(r_mega[i]["CN"]+r_mega[i]["CC"])
        r_exp = float(r_exp)/(r_mega[i]["NC"]+r_mega[i]["CC"]+r_mega[i]["CN"]+r_mega[i]["NN"])
        print r_mega
        print r_obs
        print r_exp
        print math.fabs((r_obs-r_exp)/r_exp)
        r_distro = math.log(math.fabs((r_obs-r_exp)/r_exp))

    # Calculate log-likelihood scores for each matrix
    print "Experimental Distribution of Log-likelihoods"
    for i in range(len(mega)):
        obs = mega[i]["CC"]
        exp = (mega[i]["NC"]+mega[i]["CC"])*(mega[i]["CN"]+mega[i]["CC"])
        exp = float(exp)/(mega[i]["NC"]+mega[i]["CC"]+mega[i]["CN"]+mega[i]["NN"])
        ll = math.log(math.fabs((obs-exp)/exp))
        print i, ll, r_distro
    

# Resolve edge dependencies into one of four classes
def resolve(position, matrix, anc1, anc2, dec1, dec2):
    flag_i, flag_j = "N", "N"
    if anc1.color != dec1.color: flag_i = "C"
    if anc2.color != dec2.color: flag_j = "C"
    matrix[position][flag_i+flag_j] += 1
    if flag_i == "C" and flag_j == "C": return True
    else: return False

# Verify that the alignment at this position works
def verify(alignment, key, position, skip):
    h = alignment[key]["human"]
    c = alignment[key]["chimp"]
    m = alignment[key]["mouse"]
    r = alignment[key]["rat"]
    d = alignment[key]["dog"]
    a = [ "A", "T", "G", "C" ]
    if h[position].upper() in a  and \
       c[position].upper() in a  and \
       m[position].upper() in a  and \
       r[position].upper() in a  and \
       d[position].upper() in a:
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
