import sys, random, cm_scan3, pickle, scanNonAdj2, randomScanNonAdj2

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

    #look at level 2, just the primate
    dogColor = trees[0].r.children[1].children[1].color
   
    goodSetDog = []
    interCounter = 0
    for t in trees:
        goodSetDog.append(intersection(t.r.children[1].color, dogColor))        
        if len(goodSetDog[interCounter][0]) > 0:
            pass
        else:            
            goodSetDog[interCounter][0] = t.r.children[1].children[1].color
        interCounter += 1
    #now do the coloring sets    
    newTreeCount = 0
    newTreeCount = 0
    startTreeLength = len(trees)
    
    for i in range(startTreeLength):        
        if 1 == len(goodSetDog[i][0]): #you only have one tree, is it enough?
            trees[i].r.children[1].children[1].color = goodSetDog[i][0]
        else: #need more trees
            currentTreeLength = len(trees)
            #start by assigning the one already in trees            
            trees[i].r.children[1].children[1].color = goodSetDog[i][0][0]            
            goodSetDog[i][0] = goodSetDog[i][0][1:len(goodSetDog[i][0])]
            moreIndex = 0
            for r in range(len(goodSetDog[i][0])):                                   
                trees.append(tree.copy(tree.r))
                trees[i].fullCopy(trees[currentTreeLength + moreIndex])                
                trees[currentTreeLength + moreIndex].r.children[1].children[1].color = goodSetDog[i][0][r]                
                moreIndex += 1
     
    trees = scoreIt(trees)
    return trees
##    if len(trees) > 1:
##        ri = random.randint(0, len(trees)-1)
##        return trees[ri]
##    else: return trees[0]

# Resolve edge dependencies into one of four classes
def resolve(matrix, anc1, anc2, dec1, dec2):
    flag_i, flag_j = "N", "N"
    if anc1.color != dec1.color: flag_i = "C"
    if anc2.color != dec2.color: flag_j = "C"
    matrix[flag_i+flag_j] += 1

def findLength(filename):
    f = open("../../Doc/pwm/" + filename)
    s = f.readline()
    s = f.readline()
    s = f.readline()
    count = 0
    while s != "":
        s = f.readline()
        count = count + 1
    f.close()
    return count

def main(args):
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
    a = ["A", "C", "T", "G"]
    
    #first I need to make all the trees
    trees = dict()
    for n1 in a:
        for n2 in a:
            for n3 in a:
                for n4 in a:
                    for n5 in a:
                        trees[n1+n2+n3+n4+n5] = color(supertree, [n1, n2, n3, n4, n5])

    a_file = open("../../Doc/pickle/alignment.pickle", "r")

    # Load the 5-way alignment
    alignment = pickle.load(a_file)

    # Close relevant files
    a_file.close()
    pwmFile = open("../../Doc/Vertebrate PWM/PWMlist.txt", "r")
    pwms = pickle.load(pwmFile)
    pwmFile.close()
    skip = range(1)
    pvals = [-9]
    len1 = 0    
    f1 = ""    
    for j in pvals:        
        print "Doing pval " + str(j)        
        for s in skip:            
            for i in range(len(pwms)):                
                #f1 = open("../doc/fasta/" + pwms[i] + "random(p " + str(j) + ").gff", "r")
                f1 = open("../../Doc/fasta/sortedhuman1ku-corrected.fasta." + pwms[i] + "cc.gff", "r")
                m1 = cm_scan3.buildMotifs(f1, alignment)            
                len1 = findLength(pwms[i] + "cc.pwm")            
                f1.close()
                print len(m1)
                if len1 > s+1:
                    if len(m1) != 0:
                        print pwms[i]
                        scanNonAdj2.mainRunner("sortedhuman1ku-corrected.fasta." + pwms[i] + "cc.gff", m1, len(m1), \
                                               "../../Results/Controls/Comp Permuted PWMs/Raw Data/Gap " + \
                                               str(s) + "/", s, trees)
                        #scanNonAdj2.mainRunner(pwms[i] + "(p " + str(j) + ").gff", m1, len(m1), "LL all PWM pval " + \
                        #                       str(j) + " NA " + str(s) + "/", s, trees)                    
                    #else: print "NO HITS"
            #randomScanNonAdj2.mainRunner(s, j, trees)
    sys.exit(0)
        
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
