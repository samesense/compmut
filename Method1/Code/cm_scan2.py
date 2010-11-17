#!/usr/bin/python
# Greg Donahue, 09/23/2005
# -----------------------------------------------------------------------------
import sys, re, pickle, os, random, chisquarematrix, copy, math
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

# fitchslap is the Fitch forward pass
def fitchslap(node):

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
    return

# Feeyotch is the backward pass of Fitch
def feeyotch(node):
    if len(node.children) == 0: return
    node.color = node.color[0]
    for c in node.children: feeyotch(c)
    return

# Color a new tree given a tree and some submotif
def color(tree, submotif):
    
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
    
    # Create the positional uberdictionaries and interdirectories
    uberlist = []
    for i in range(len(motifs[0][0])):
        cdict = dict()
        alphabet = [ "A->A", "A->C", "A->G", "A->T",
                     "C->A", "C->C", "C->G", "C->T",
                     "G->A", "G->C", "G->G", "G->T",
                     "T->A", "T->C", "T->G", "T->T" ]
        for edge in alphabet: cdict[edge] = 0
        d = dict()
        d["root->rodent"] = copy.copy(cdict)
        d["root->highermammal"] = copy.copy(cdict)
        d["rodent->rat"] = copy.copy(cdict)
        d["rodent->mouse"] = copy.copy(cdict)
        d["highermammal->dog"] = copy.copy(cdict)
        d["highermammal->primate"] = copy.copy(cdict)
        d["primate->human"] = copy.copy(cdict)
        d["primate->chimp"] = copy.copy(cdict)
        uberlist.append(d)
    interlist = []
    for i in range(len(motifs[0][0])-1):
        cdict = dict()
        alphabet = [ "A->A", "A->C", "A->G", "A->T",
                     "C->A", "C->C", "C->G", "C->T",
                     "G->A", "G->C", "G->G", "G->T",
                     "T->A", "T->C", "T->G", "T->T" ]
        for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                cdict[alphabet[i]+"/"+alphabet[j]] = 0
        d = dict()
        d["root->rodent"] = copy.copy(cdict)
        d["root->highermammal"] = copy.copy(cdict)
        d["rodent->rat"] = copy.copy(cdict)
        d["rodent->mouse"] = copy.copy(cdict)
        d["highermammal->dog"] = copy.copy(cdict)
        d["highermammal->primate"] = copy.copy(cdict)
        d["primate->human"] = copy.copy(cdict)
        d["primate->chimp"] = copy.copy(cdict)
        interlist.append(d)
    # Keeps track of conditional root probabilities
    rlist = []
    for i in range(len(motifs[0][0])-1):
        alphabet = [ "A|A", "A|C", "A|G", "A|T",
                     "C|A", "C|C", "C|G", "C|T",
                     "G|A", "G|C", "G|G", "G|T",
                     "T|A", "T|C", "T|G", "T|T" ]
        d = dict()
        for i in alphabet: d[i] = 0
        rlist.append(d)
    # Keeps track of background probabilities
    blist = dict()
    alphabet = [ "A", "T", "C", "G" ]
    for i in alphabet: blist[i] = 0

    # We want to randomly split the instances into training and test sets
    # And we want to repeat this for random splits until we get it right
    incrementor = 0
    logs = range(len(motifs[0][0])-1)
    while True:
        
        # Split the instance set randomly
        training = motifs[:int(random.uniform(0, 1)*len(motifs))]
        test = motifs[int(random.uniform(0, 1)*len(motifs)):]
        random.shuffle(motifs)

        if len(training) == 0 or len(test) == 0: continue

        # For each instance:
        for m in training:
            
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
            
            # For every pair of adjacent trees, update the ubers/inters
            for i in range(len(m[0])-1):
                t1 = trees[i]
                t2 = trees[i+1]
                # Update the ubers
                roo = t1.r
                rod = roo.children[0]
                hig = roo.children[1]
                mou = rod.children[0]
                rat = rod.children[1]
                dog = hig.children[0]
                pri = hig.children[1]
                hum = pri.children[0]
                chi = pri.children[1]
                #edge counts for a position in the alignment
                #so this data involves one position
                uberlist[i]["root->rodent"][roo.color+"->"+rod.color] += 1
                uberlist[i]["root->highermammal"][roo.color+"->"+hig.color] += 1
                uberlist[i]["rodent->mouse"][rod.color+"->"+mou.color] += 1
                uberlist[i]["rodent->rat"][rod.color+"->"+rat.color] += 1
                uberlist[i]["highermammal->dog"][hig.color+"->"+dog.color] += 1
                uberlist[i]["highermammal->primate"][hig.color+"->"+pri.color] += 1
                uberlist[i]["primate->human"][pri.color+"->"+hum.color] += 1
                uberlist[i]["primate->chimp"][pri.color+"->"+chi.color] += 1                
                # Update the inters
                roo2 = t2.r
                rod2 = roo2.children[0]
                hig2 = roo2.children[1]
                mou2 = rod2.children[0]
                rat2 = rod2.children[1]
                dog2 = hig2.children[0]
                pri2 = hig2.children[1]
                hum2 = pri2.children[0]
                chi2 = pri2.children[1]
                #edge counts for adjacent positions in the alignment
                #this data involves two positions
                interlist[i]["root->rodent"][roo.color+"->"+rod.color+"/"+roo2.color+"->"+rod2.color] += 1
                interlist[i]["root->highermammal"][roo.color+"->"+hig.color+"/"+roo2.color+"->"+hig2.color] += 1
                interlist[i]["rodent->mouse"][rod.color+"->"+mou.color+"/"+rod2.color+"->"+mou2.color] += 1
                interlist[i]["rodent->rat"][rod.color+"->"+rat.color+"/"+rod2.color+"->"+rat2.color] += 1
                interlist[i]["highermammal->dog"][hig.color+"->"+dog.color+"/"+hig2.color+"->"+dog2.color] += 1
                interlist[i]["highermammal->primate"][hig.color+"->"+pri.color+"/"+hig2.color+"->"+pri2.color] += 1
                interlist[i]["primate->human"][pri.color+"->"+hum.color+"/"+pri2.color+"->"+hum2.color] += 1
                interlist[i]["primate->chimp"][pri.color+"->"+chi.color+"/"+pri2.color+"->"+chi2.color] += 1
                # Update the rlist
                rlist[i][roo.color+"|"+roo2.color] += 1
                # Update the blist
                blist[roo.color] += 1

        # Now turn those counts into ratios
        for i in range(len(training[0][0])-1):
            #find the prob of seeing an edge at a certain position
            #does not calculate seeing a node given another node in the same edge
            for j in uberlist[i].keys():
                ubersum = 0
                for k in uberlist[i][j].keys(): ubersum += uberlist[i][j][k]
                for k in uberlist[i][j].keys():
                    uberlist[i][j][k] = float(uberlist[i][j][k])/ubersum
            #also finds the probabilities of seeing say, ACTG, for the two positions
            #does not find conditional probs
            #wait, make sure it doesn't do cond probs
            #Greg is counting the probs of seeing AC and TG
            #he could use this later for the conditionals, i guess
            for j in interlist[i].keys():
                intersum = 0
                for k in interlist[i][j].keys(): intersum += interlist[i][j][k]
                for k in interlist[i][j].keys():
                    interlist[i][j][k] = float(interlist[i][j][k])/intersum
                rsum = 0
            for j in rlist[i].keys(): rsum += rlist[i][j]
            for j in rlist[i].keys(): rlist[i][j] = float(rlist[i][j])/rsum
        bsum = 0
        for i in blist.keys(): bsum += blist[i]
        for i in blist.keys(): blist[i] = float(blist[i])/bsum

        # Now that we have estimated parameters, calculate the log-likelihood
        # ratio of each tree pair in each instance and average
        ratios = []
        for i in range(len(training[0][0])-1): ratios.append(0)
    
        # For each motif:
        for m in test:

            # Check the motif: is it any good?
            flag = True
            for subm in m:
                if len(subm) != len(m[0]): flag = False
                if "N" in subm: flag = False
            if not flag: continue

            # For each position in the hit, color a new tree
            trees = []
            for position in range(len(m[0])):
                submotif = []
                for i in range(len(m)): submotif.append(m[i][position])
                trees.append(color(supertree, submotif))
            
            # For every pair of adjacent trees, update the ubers/inters
            for i in range(len(m[0])-1):
                t1 = trees[i]
                t2 = trees[i+1]
                roo = t1.r
                rod = roo.children[0]
                hig = roo.children[1]
                mou = rod.children[0]
                rat = rod.children[1]
                dog = hig.children[0]
                pri = hig.children[1]
                hum = pri.children[0]
                chi = pri.children[1]
                roo2 = t2.r
                rod2 = roo2.children[0]
                hig2 = roo2.children[1]
                mou2 = rod2.children[0]
                rat2 = rod2.children[1]
                dog2 = hig2.children[0]
                pri2 = hig2.children[1]
                hum2 = pri2.children[0]
                chi2 = pri2.children[1]
                
                # Compute the probability of seeing t1
                #why do it this way?
                #this does not take dependencies between nodes in an edge into account
                #maybe this is the way it should be done?
                #yeah, this does does not look at the node before it
                #the nodes are assumed to be independent for some reason,
                #yet the prob of seeing a node considers the node before it
                roorod = uberlist[i]["root->rodent"]["A->"+rod.color]+\
                         uberlist[i]["root->rodent"]["G->"+rod.color]+\
                         uberlist[i]["root->rodent"]["C->"+rod.color]+\
                         uberlist[i]["root->rodent"]["T->"+rod.color]
                roohig = uberlist[i]["root->highermammal"]["A->"+hig.color]+\
                         uberlist[i]["root->highermammal"]["G->"+hig.color]+\
                         uberlist[i]["root->highermammal"]["C->"+hig.color]+\
                         uberlist[i]["root->highermammal"]["T->"+hig.color]
                rodmou = uberlist[i]["rodent->mouse"]["A->"+mou.color]+\
                         uberlist[i]["rodent->mouse"]["G->"+mou.color]+\
                         uberlist[i]["rodent->mouse"]["C->"+mou.color]+\
                         uberlist[i]["rodent->mouse"]["T->"+mou.color]
                rodrat = uberlist[i]["rodent->rat"]["A->"+rat.color]+\
                         uberlist[i]["rodent->rat"]["G->"+rat.color]+\
                         uberlist[i]["rodent->rat"]["C->"+rat.color]+\
                         uberlist[i]["rodent->rat"]["T->"+rat.color]
                higdog = uberlist[i]["highermammal->dog"]["A->"+dog.color]+\
                         uberlist[i]["highermammal->dog"]["G->"+dog.color]+\
                         uberlist[i]["highermammal->dog"]["C->"+dog.color]+\
                         uberlist[i]["highermammal->dog"]["T->"+dog.color]
                higpri = uberlist[i]["highermammal->primate"]["A->"+pri.color]+\
                         uberlist[i]["highermammal->primate"]["G->"+pri.color]+\
                         uberlist[i]["highermammal->primate"]["C->"+pri.color]+\
                         uberlist[i]["highermammal->primate"]["T->"+pri.color]
                prichi = uberlist[i]["primate->chimp"]["A->"+chi.color]+\
                         uberlist[i]["primate->chimp"]["G->"+chi.color]+\
                         uberlist[i]["primate->chimp"]["C->"+chi.color]+\
                         uberlist[i]["primate->chimp"]["T->"+chi.color]
                prihum = uberlist[i]["primate->human"]["A->"+hum.color]+\
                         uberlist[i]["primate->human"]["G->"+hum.color]+\
                         uberlist[i]["primate->human"]["C->"+hum.color]+\
                         uberlist[i]["primate->human"]["T->"+hum.color]
                p_t1 = blist[roo.color]*roorod*roohig*rodmou*rodrat*higdog*higpri*prichi*prihum
                #p_t1 = blist[roo.color]* \
                #       uberlist[i]["root->rodent"][roo.color+"->"+rod.color]*uberlist[i]["root->highermammal"][roo.color+"->"+hig.color]*uberlist[i]["rodent->mouse"][rod.color+"->"+mou.color]*uberlist[i]["rodent->rat"][rod.color+"->"+rat.color]*uberlist[i]["highermammal->dog"][hig.color+"->"+dog.color]*uberlist[i]["highermammal->primate"][hig.color+"->"+pri.color]*uberlist[i]["primate->chimp"][pri.color+"->"+chi.color]*uberlist[i]["primate->human"][pri.color+"->"+hum.color]
                
                # Compute the probability of seeing t2 | t1
                #this also assumes the independence of nodes
                #the two positions don't involve conditionals, it is
                #doing the prob of the intersection
                roorod = getprobs(interlist[i]["root->rodent"], rod.color, rod2.color)
                roohig = getprobs(interlist[i]["root->highermammal"], hig.color, hig2.color)
                rodmou = getprobs(interlist[i]["rodent->mouse"], mou.color, mou2.color)
                rodrat = getprobs(interlist[i]["rodent->rat"], rat.color, rat2.color)
                higdog = getprobs(interlist[i]["highermammal->dog"], dog.color, dog2.color)
                higpri = getprobs(interlist[i]["highermammal->primate"], pri.color, pri2.color)
                prichi = getprobs(interlist[i]["primate->chimp"], chi.color, chi2.color)
                prihum = getprobs(interlist[i]["primate->human"], hum.color, hum2.color)
                p_t2 = rlist[i][roo.color+"|"+roo2.color]*roorod*roohig*rodmou*rodrat*higdog*higpri*prichi*prihum
                #p_t2 = rlist[i][roo.color+"|"+roo2.color]* \
                #       interlist[i]["root->rodent"][roo.color+"->"+rod.color+"/"+roo2.color+"->"+rod2.color]*interlist[i]["root->highermammal"][roo.color+"->"+hig.color+"/"+roo2.color+"->"+hig2.color]*interlist[i]["rodent->mouse"][rod.color+"->"+mou.color+"/"+rod2.color+"->"+mou2.color]*interlist[i]["rodent->rat"][rod.color+"->"+rat.color+"/"+rod2.color+"->"+rat2.color]*interlist[i]["highermammal->dog"][hig.color+"->"+dog.color+"/"+hig2.color+"->"+dog2.color]*interlist[i]["highermammal->primate"][hig.color+"->"+pri.color+"/"+hig2.color+"->"+pri2.color]*interlist[i]["primate->chimp"][pri.color+"->"+chi.color+"/"+pri2.color+"->"+chi2.color]*interlist[i]["primate->human"][pri.color+"->"+hum.color+"/"+pri2.color+"->"+hum2.color]
                
                # Compute the log-likelihood ratio
                if p_t2 == 0 or p_t1 == 0: ratios[i] += 0
                else: ratios[i] += -1*math.log(p_t2/p_t1)

        #for i in range(len(ratios)): print i, i+1, ratios[i]
        for i in range(len(ratios)): logs[i] += ratios[i]

        incrementor += 1
        if incrementor >= 200: break

    # Print the log-sums
    for i in range(len(logs)): print i+1, i+2, logs[i]

    # Prematurely terminate - let's not worry about X^2 analysis right now
    sys.exit(0)
    
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

def getprobs(inter, c1, c2):
    ret = 0
    alphabet = [ "A", "T", "C", "G" ]
    for i in alphabet:
        for j in alphabet:
            ret += inter[i+"->"+c1+"/"+j+"->"+c2]
    return ret

# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
