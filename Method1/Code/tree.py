#!/usr/bin/python
# tree.py contains Tree and TreeNode classes and some pertinent functions
# Greg Donahue, 10/31/2005
# -----------------------------------------------------------------------------
import sys
# -----------------------------------------------------------------------------
# MODULE CLASSES
# Tree class represents a binary, back-linked tree
class Tree:

    # Instance variables
    # self.root is the root node of this tree
    
    # Instantiate new Trees
    def __init__(self):
        self.root = None

    # Return a string representation of this Tree
    def __repr__(self):
        return str(self.getEdges())

    # Return a collection of all the edges in this Tree
    # Edges are represented as coordinate pairs of TreeNodes
    # The getEdges function is a wrapper around a recursive function
    def getEdges(self):
        return self.recursivelyGetEdges(self.root)
    def recursivelyGetEdges(self, node):
        if node.left == None and node.right == None: return []
        ret = []
        if node.left != None:
            ret.append((node, node.left))
            ret.extend(self.recursivelyGetEdges(node.left))
        if node.right != None:
            ret.append((node, node.right))
            ret.extend(self.recursivelyGetEdges(node.right))
        return ret

    # Return a collection of all the TreeNodes in this Tree
    # The getNodes functions is a wrapper around a recursive function
    def getNodes(self):
        return self.recursivelyGetNodes(self.root)
    def recursivelyGetNodes(self, node):
        if node.left == None and node.right == None: return [ node ]
        ret = [ node ]
        if node.left != None: ret.extend(self.recursivelyGetNodes(node.left))
        if node.right != None: ret.extend(self.recursivelyGetNodes(node.right))
        return ret

    # The fitch function colors all the nodes of this Tree given a coloring of
    # the taxa; the coloring is assumed to describe taxon states going from the
    # leftmost taxon to the rightmost taxon
    def fitch(self, coloring):
        self.colorTaxa(coloring)
        self.colorAncestors()

    # Color the taxa of this Tree with the input coloring
    # The colorTaxa function is a wrapper around a recursive function
    def colorTaxa(self, coloring):
        self.recursivelyColorTaxa(self.root, 0, coloring)
    def recursivelyColorTaxa(self, node, index, coloring):
        if node.left == None and node.right == None:
            node.datum = Datum(coloring[index])
            return index+1
        index = self.recursivelyColorTaxa(node.left, index, coloring)
        index = self.recursivelyColorTaxa(node.right, index, coloring)
        return index

    # Color the ancestors of this Tree using state-assignments at the taxa
    # This is the classical weighted parsimony algorithm (Durbin et al, p174)
    # The colorAncestors function is a wrapper around a recursive function
    def colorAncestors(self):
        self.recursivelyColorAncestors(self.root)
    def recursivelyColorAncestors(self, node):
        if node.left == None and node.right == None:
            node.scores[node.datum.color] = 0
            return
        best_state = "A"
        min_score = 100
        for state in [ "A", "T", "G", "C" ]:
            

# TreeNode class represents nodes in a binary, back-linked tree
class TreeNode:

    # Instance variables
    # self.parent is the parent of this node
    # self.left is the left child of this node
    # self.right is the right child of this node
    # self.datum is a data object bound to this node
    
    # Initialize new TreeNodes
    def __init__(self, p, d, side):
        self.parent = p
        self.datum = d
        self.left, self.right = None, None
        if side == "left":
            p.left = self
        elif side == "right":
            p.right = self

    # Return a string representation of this TreeNode
    def __repr__(self):
        return str(self.datum)

# Datum class represents a set of scores and a state for a given node
class Datum:

    # Instance variables
    # self.color is the state at this node
    # self.scores is a dictionary of the costs putting states at this node

    # Instantiate new Datums
    def __init__(self, state):
        self.color = state
        self.scores = dict()
        self.scores["A"] = None
        self.scores["T"] = None
        self.scores["G"] = None
        self.scores["C"] = None

    # Return a string representation of this Datum
    def __repr__(self):
        return self.color

# -----------------------------------------------------------------------------
# MODULE FUNCTIONS
# Main function
def main(args):

    # Build a Tree
    t = Tree()
    rt = TreeNode(None, None, None)
    hm = TreeNode(rt, None, side="right")
    ro = TreeNode(rt, None, side="left")
    pr = TreeNode(hm, None, side="right")
    do = TreeNode(hm, None, side="left")
    hu = TreeNode(pr, None, side="right")
    ch = TreeNode(pr, None, side="left")
    mo = TreeNode(ro, None, side="right")
    ra = TreeNode(ro, None, side="left")
    t.root = rt

    # Print its nodes and edges
    print "Nodes", t.getNodes()
    print "Edges", t.getEdges()

    # Run the Fitch algorithm on the following colorings:
    coloring = "AATAG"
    print "The Fitch coloring of this tree with", coloring
    t.fitch(coloring)
    
# -----------------------------------------------------------------------------
# The following code executes upon command-line invocation
if __name__ == "__main__": main(sys.argv)

# -----------------------------------------------------------------------------
# EOF
