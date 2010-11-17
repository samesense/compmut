#!/usr/bin/python
# Greg Donahue, 09/23/2005
# -----------------------------------------------------------------------------
import sys, re, pickle, os, random, copy, math
from Numeric import array, zeros
from rpy import *

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

#this must change for position I think
def getExp(mega, position):
    return float((mega[position]["CN"]+mega[position]["CC"])*(mega[position]["NC"]+mega[position]["CC"]))/ \
           (mega[position]["CN"]+mega[position]["NC"]+mega[position]["CC"]+mega[position]["NN"])

def make16Sq():
    r = dict()
    nucs = ["A", "C", "G", "T"]
    edges = []
    for n1 in nucs:
        for n2 in nucs:
            r[n1+n2] = float(0)    
    return r

def make256Sq():
    r = dict()
    nucs = ["A", "C", "G", "T"]
    edges = []
    for n1 in nucs:
        for n2 in nucs:
            edges.append(n1+n2)
    for e1 in edges:
        for e2 in edges:
            r[e1+e2] = float(0)
    return r

def checkConservation(submotif):
    nucs = []
    Ascore = 0
    Cscore = 0
    Gscore = 0
    Tscore = 0
    
    for i in range(len(["human", "chimp", "dog", "rat", "mouse"])):        
        nucs.append(submotif[i])
    for n in nucs:
        if n == "A": Ascore += 1
        elif n == "C": Cscore += 1
        elif n == "G": Gscore += 1
        else: Tscore += 1    
    m = max([Ascore, Cscore, Gscore, Tscore])
    #print submotif, m
    if m > 4: return 1
    else: return 0    

#########################################
################################################
def mainRunner(argPWM, motifs, hitLimit, writeFile, skip, treeList):
    pwm_file = argPWM
##    treeFile = open("trees.txt", "r")
##    treeList = pickle.load(treeFile)
##    treeFile.close()
    # Hard-coded values
    p = -12.0
    #fasta = "../doc/fasta/human1ku-corrected.fasta"
    #a_file = open("../doc/pickle/alignment.pickle", "r")

    # Load the 5-way alignment
    #alignment = pickle.load(a_file)

    # Close relevant files
    #a_file.close()

    # Run PWM_SCAN on the human FASTA, generate a GFF file
    # (Actually, the program assumes that you have already done this separately
    # and just goes ahead and tries to find the GFF file)
    pwm = pwm_file.replace(".pwm","")
    #gff = open("../doc/fasta/" + pwm); 

    # For each hit in the GFF, throw out hits that span gaps, build motifs
    #motifs = buildMotifs(gff, alignment)
    if len(motifs) == 0:
        print "NO HITS"
        return
    # used in analysis of distributions PE
    # can be commented out if needed
    #takes top 1000
    if len(motifs) > 1000: motifs = motifs[0:1000]
    #if len(motifs) > hitLimit-1: motifs = motifs[0:hitLimit]
    #else : print "Did not reach limit " + str(hitLimit) + " with " + str(len(motifs))

    # Create the change/no-change matrix
    mega = []
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
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
    for i in range(len(motifs[0][0])-1-skip):
        cmuts = dict()
        alphabet = [ "N", "C" ]
        for i in alphabet:
            for j in alphabet: cmuts[i+j] = 0
        cmuts["CC"] = float(0);
        cmuts["CN"] = float(5);
        cmuts["NC"] = float(6.5);
        cmuts["NN"] = float(788.5);
        primateChimpMega.append(cmuts)

    #This keeps track of all 16 nt combination counts
    #between two positions, for each taxon
    #The order in the submotif is human, chimp, mouse, rat, dog
    #each 16Sq is indexed by pos, pos+1
    labels = ["human", "chimp", "mouse", "rat", "dog"]
    big16 = dict()
    for label in labels:
        big16[label] = []    
        for i in range(len(motifs[0][0])-1-skip):
            big16[label].append(make16Sq())

    #This keeps track of all 256 possible changes
    big256 = []
    for i in range(len(motifs[0][0])-1-skip):
            big256.append(make256Sq())

    conserved = []
    for i in range(len(motifs[0][0])):
        conserved.append(0)
    motifCount = 0    
    # For each instance
    for m in motifs:
        
        # Build a set of trees out of this motif        
        if not validate(m): continue
        #print motifCount
        motifCount = motifCount + 1
        trees = []
        for position in range(len(m[0])):
            submotif = ""
            submotif2 = ""
            for i in range(len(m)): submotif += m[i][position]
            conserved[position] += checkConservation(submotif)
            if position < len(m[0])-1-skip:
                for i in range(len(m)): submotif2 += m[i][position+1+skip]
            trees.append(treeList[submotif])
            #submotif goes 
            #Fill in big16
            #The order in the submotif is human, chimp, mouse, rat, dog
            #chances are that I'll only end up using the human counts
            #but I'll have the others just in case.
            if position < len(m[0])-1-skip:
                for u in range(len(labels)):                    
                    big16[labels[u]][position][submotif[u]+submotif2[u]] += 1
            
##            submotif = []
##            for i in range(len(m)): submotif.append(m[i][position])
##            trees.append(color(supertree, submotif))
            
        # Walk down pairs of trees, fill in C/NC matrix
        #change position range to determine other relations
        
        for position in range(len(m[0])-1-skip):
            t1_set = trees[position]
            print position+skip+1
            t2_set = trees[position+skip+1]
            
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
                    resolve(position, mega, rootRodentMega, t1_r, t2_r, t1_rodent, t2_rodent, weight, big256)
                    resolve(position, mega, rodentMouseMega, t1_rodent, t2_rodent, t1_mouse, t2_mouse, weight, big256)
                    three = resolve(position, mega, rodentRatMega, t1_rodent, t2_rodent, t1_rat, t2_rat, weight, big256)
                    resolve(position, mega, rootHMMega, t1_r, t2_r, \
                            t1_highermammal, t2_highermammal, weight, big256)
                    resolve(position, mega, HMPrimateMega, t1_highermammal, t2_highermammal, \
                            t1_primate, t2_primate, weight, big256)
                    resolve(position, mega, HMDogMega, t1_highermammal, t2_highermammal, \
                            t1_dog, t2_dog, weight, big256)
                    resolve(position, mega, primateHumanMega, t1_primate, t2_primate, t1_human, t2_human, weight, big256)
                    resolve(position, mega, primateChimpMega, t1_primate, t2_primate, t1_chimp, t2_chimp, weight, big256)
                    

    newFile = open(writeFile + pwm + ".ll", "w")
    for i in range(len(mega)):        
        obs = mega[i]["CC"]
        exp = float(0)
        exp += getExp(rootRodentMega, i)
        exp += getExp(rootHMMega, i)
        exp += getExp(rodentRatMega, i)
        exp += getExp(rodentMouseMega, i)
        exp += getExp(HMDogMega, i)
        exp += getExp(HMPrimateMega, i)
        exp += getExp(primateChimpMega, i)
        exp += getExp(primateHumanMega, i)        
        ll = math.log(float(obs)/exp)            
        newFile.write("(" + str(i+1) +  ", " + str(i+2+skip) + ") " +  str(ll) + " \n")
                #print "(", i+1, ",", i+2, ")", ll, r.pchisq(float(obs)/exp, 2, lower_tail=False)
                #print "(", i+1, ",", i+2, ")", ll, r.pchisq(float(obs)/exp, 2, lower_tail=False)
        
        #    newFile.write("(" + str(i+1) +  ", " + str(i+2) + ") " +  str(1) + " \n")
    newFile.write("Hit Count: " + str(motifCount))
    newFile.close()

    #ok, let's pickle the big16 so that it is per PWM and per NA
    #it will go in the correct folder
##    newFile = open(writeFile + pwm + ".16sq", "w")
##    pickle.dump(big16, newFile)
##    newFile.close()

##    newFile = open(writeFile + pwm + ".256sq", "w")
##    pickle.dump(big256, newFile)
##    newFile.close()
    i = 1;
    for p in conserved:
        print "pos " + str(i) + ":", float(p)/float(motifCount)
        i += 1
    
    #sys.exit(0)

###########################################################
######################################################

# Main function
def main(args):
    sys.exit(0)    

# Resolve edge dependencies into one of four classes
def resolve(position, matrix, matrix2, anc1, anc2, dec1, dec2, weight, big256):
    flag_i, flag_j = "N", "N"
    if anc1.color != dec1.color: flag_i = "C"
    if anc2.color != dec2.color: flag_j = "C"
    matrix[position][flag_i+flag_j] += weight
    matrix2[position][flag_i+flag_j] += weight
    big256[position][anc1.color+dec1.color+anc2.color+dec2.color] += weight

# Verify that the alignment at this position works
def verify(alignment, key, position):
    h = alignment[key]["human"]
    c = alignment[key]["chimp"]
    m = alignment[key]["mouse"]
    r = alignment[key]["rat"]
    d = alignment[key]["dog"]
    a = [ "A", "T", "G", "C" ]
    if h[position].upper() in a and h[position+1].upper() in a and \
       c[position].upper() in a and c[position+1].upper() in a and \
       m[position].upper() in a and m[position+1].upper() in a and \
       r[position].upper() in a and r[position+1].upper() in a and \
       d[position].upper() in a and d[position+1].upper() in a:
        return True
    else: return False

# Print a tree recursively
def treeprint(treenode):
    print treenode.color
    if len(treenode.children) == 0: return
    treeprint(treenode.children[0])
    treeprint(treenode.children[1])

# Validate the motif
#just checks the lengths of the alignments against
#human alignment
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
