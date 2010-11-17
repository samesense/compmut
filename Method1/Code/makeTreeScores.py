import sys, pickle, random

#This will make trees with all possible colorings
#and compare them to make a master list of
# change/nochange matricies.
#The list is called treeScores.txt

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
    if len(trees) > 1:
        ri = random.randint(0, len(trees)-1)
        return trees[ri]
    else: return trees[0]

# Resolve edge dependencies into one of four classes
def resolve(matrix, anc1, anc2, dec1, dec2):
    flag_i, flag_j = "N", "N"
    if anc1.color != dec1.color: flag_i = "C"
    if anc2.color != dec2.color: flag_j = "C"
    matrix[flag_i+flag_j] += 1

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

    print "done with the trees"
##    count = 0                                 
##    for t1 in trees:
##        for t2 in trees:
##            print count
##            count += 1
##            #make the scoring matricies for this combination of trees
##            mat = dict()
##            mat["master"] = dict()
##            mat["rootRodent"] = dict()
##            mat["rodentMouse"] = dict()
##            mat["rodentRat"] = dict()
##            mat["primateHuman"] = dict()
##            mat["primateChimp"] = dict()
##            mat["rootHM"] = dict()
##            mat["HMPrimate"] = dict()
##            mat["HMDog"] = dict()
##
##            m_keys = mat.keys()
##            for m in m_keys:
##                mat[m]["CC"] = 0
##                mat[m]["CN"] = 0
##                mat[m]["NC"] = 0
##                mat[m]["NN"] = 0                                            
##
##            # Decompose the tree at position i                                            
##            t1_r = t1.r
##            t1_rodent = t1_r.children[0]
##            t1_highermammal = t1_r.children[1]
##            t1_mouse = t1_rodent.children[0]
##            t1_rat = t1_rodent.children[1]
##            t1_dog = t1_highermammal.children[0]
##            t1_primate = t1_highermammal.children[1]
##            t1_human = t1_primate.children[0]
##            t1_chimp = t1_primate.children[1]
##                
##            # Decompose the tree at position i+2                                            
##            t2_r = t2.r
##            t2_rodent = t2_r.children[0]
##            t2_highermammal = t2_r.children[1]
##            t2_mouse = t2_rodent.children[0]
##            t2_rat = t2_rodent.children[1]
##            t2_dog = t2_highermammal.children[0]
##            t2_primate = t2_highermammal.children[1]
##            t2_human = t2_primate.children[0]
##            t2_chimp = t2_primate.children[1]
##
##            #now score the two trees
##            resolve(mat["master"], t1_r, t2_r, t1_rodent, t2_rodent)
##            resolve(mat["rootRodent"], t1_r, t2_r, t1_rodent, t2_rodent)
##
##            resolve(mat["master"], t1_rodent, t2_rodent, t1_mouse, t2_mouse)
##            resolve(mat["rodentMouse"], t1_rodent, t2_rodent, t1_mouse, t2_mouse)
##
##            resolve(mat["master"], t1_rodent, t2_rodent, t1_rat, t2_rat)
##            resolve(mat["rootRodent"], t1_rodent, t2_rodent, t1_rat, t2_rat)
##
##            resolve(mat["master"], t1_r, t2_r, \
##                    t1_highermammal, t2_highermammal)
##            resolve(mat["rootHM"], t1_r, t2_r, \
##                    t1_highermammal, t2_highermammal)
##
##            resolve(mat["master"], t1_highermammal, t2_highermammal, \
##                    t1_primate, t2_primate)
##            resolve(mat["HMPrimate"], t1_highermammal, t2_highermammal, \
##                    t1_primate, t2_primate)
##
##            resolve(mat["master"], t1_highermammal, t2_highermammal, \
##                    t1_dog, t2_dog)
##            resolve(mat["HMDog"], t1_highermammal, t2_highermammal, \
##                    t1_dog, t2_dog)
##
##            resolve(mat["master"], t1_primate, t2_primate, t1_human, t2_human)
##            resolve(mat["primateHuman"], t1_primate, t2_primate, t1_human, t2_human)
##
##            resolve(mat["master"], t1_primate, t2_primate, t1_chimp, t2_chimp)
##            resolve(mat["primateChimp"], t1_primate, t2_primate, t1_chimp, t2_chimp)
##
##            ret[t1_rat.color+t1_mouse.color+t1_dog.color+t1_chimp.color+t1_human.color+":"+\
##                t2_rat.color+t2_mouse.color+t2_dog.color+t2_chimp.color+t2_human.color] = mat
##
##    f = open("treeScores.txt", "w")
##    pickle.dump(ret, f)                            
    
    f = open("trees.txt", "w")
    pickle.dump(trees, f)
    print "done"
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)
# -----------------------------------------------------------------------------
# EOF
