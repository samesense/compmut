#!/usr/bin/python
# Greg Donahue, 12/06/2005
# scanner.py scans a PWM/alignment for conditional mutations
# -----------------------------------------------------------------------------
import sys, pickle, random, math, copy
import cm_scan
# -----------------------------------------------------------------------------
# MODULE FUNCTIONS
# Main function
def main(args):

    # Parse arguments, get the GFF file
    if len(args) != 2: raise Exception("scanner.py GFF_FILE")
    gff_file = open(args[1], "r")

    # Load the alignment
    alignment_file = open("../doc/pickle/alignment.pickle", "r")
    alignment = pickle.load(alignment_file)
    alignment_file.close()

    # Build aligned sequences (motifs) from the GFF hits
    # All motifs are arranged like this: human, chimp, mouse, rat, dog
    motifs = cm_scan.buildMotifs(gff_file, alignment)

    # Build a set of instances
    # Each instance is a sequence of k trees, (k = length of words in motif)
    instances = []
    for m in motifs:
        pass
    #    instances.extend(build_instances(m))

    # Start iteratively looping
    incrementor = 0
    while incrementor < 1000:

        # Randomly split the instances into a training set and a test set
        cutoff = int(random.uniform(0, 1)*len(instances))
        training, test = instances[:cutoff], instances[cutoff:]
        random.shuffle(instances)

        # Increment and continue
        incrementor += 1

    # Print the summed log-likelihoods over all positions

    # Terminate gracefully
    sys.exit(0)

# Build instances from a motif
def build_instances(motif):

    # Create a tree concordance
    concordance = []

    # Gather trees for each position
    for position in range(len(motif[0])):
        scaffold = Tree()
        trees = scaffold.color([ motif[0][position], motif[1][position],
                                 motif[2][position], motif[3][position],
                                 motif[4][position] ])
        concordance.append(trees)

    # Create an instance for every combination of trees
    return resolve(concordance, 0)

# Resolve the concordance into a list of all possible instances (recursive)
#def resolve(concordance, i):
#    if i == len(concordance):
#        return []
#    for tree in concordance[i]:
#        .append(tree)

# -----------------------------------------------------------------------------
# MODULE CLASSES
# Tree class
class Tree:

    # Note: these are very specific trees - they look like this:
    #
    #                       root
    #                       /  \
    #                      /    \
    #                     /      higher mammal
    #                    /              /    \
    #                   /              /      \
    #              rodent             /        primate
    #              /    \            /         /     \
    #             /      \          /         /       \
    #         mouse       rat     dog     chimp        human

    # Instance variables
    # Tree.root = root
    # Tree.rodent = rodent
    # Tree.higher = higher mammal
    # Tree.primate = primate
    # Tree.mouse = mouse
    # Tree.rat = rat
    # Tree.dog = dog
    # Tree.chimp = chimp
    # Tree.human = human

    # Initialize new Trees:
    def __init__(self):
        pass

    # Color this Tree, return a list of Trees
    # We expect to get colors in the order: human, chimp, mouse, rat, dog
    def color(self, colors):

        # Color the taxa
        self.mouse = colors[2]
        self.rat = colors[3]
        self.dog = colors[4]
        self.chimp = colors[1]
        self.human = colors[0]

        # Color the second-level ancestor nodes
        if self.mouse == self.rat: self.rodent = [ self.mouse ]
        else: self.rodent = [ self.mouse, self.rat ]
        if self.human == self.chimp: self.primate = [ self.human ]
        else: self.primate = [ self.human, self.chimp ]

        # Color the third-level ancestor nodes
        if self.dog in self.primate: self.higher = copy.copy(self.primate)
        else: self.higher = [ self.dog ]; self.higher.extend(self.primate)

        # Color the root node
        self.root = copy.copy(self.rodent)
        for i in self.higher:
            if not i in self.rodent: self.root.append(i)

        # Now generate all possible trees from the root down
        ret = []
        for i1 in self.root:

            # Create the tree scaffold, fix taxa
            t = Tree()
            t.human = self.human
            t.chimp = self.chimp
            t.mouse = self.mouse
            t.rat = self.rat
            t.dog = self.dog
            t.root = i1

            # Fix ancestral nodes
            if t.root in self.rodent:
                t.rodent = t.root
                if t.root in self.higher:
                    t.higher = t.root
                    if t.higher in self.primate:
                        t.primate = t.higher
                        ret.append(t)
                    else:
                        for i4 in self.primate:
                            t.primate = i4
                            ret.append(t)
                else:
                    for i3 in self.higher:
                        t.higher = i3
                        if t.higher in self.primate:
                            t.primate = t.higher
                            ret.append(t)
                        else:
                            for i4 in self.primate:
                                t.primate = i4
                                ret.append(t)
            else:
                for i2 in self.rodent:
                    t.rodent = i2
                    if t.root in self.higher:
                        t.higher = t.root
                        if t.higher in self.primate:
                            t.primate = t.higher
                            ret.append(t)
                        else:
                            for i4 in self.primate:
                                t.primate = i4
                                ret.append(t)
                    else:
                        for i3 in self.higher:
                            t.higher = i3
                            if t.higher in self.primate:
                                t.primate = t.higher
                                ret.append(t)
                            else:
                                for i4 in self.primate:
                                    t.primate = i4
                                    ret.append(t)

        # Return the list of Trees created by this coloring
        return ret

# -----------------------------------------------------------------------------
# The following code is executed upon command-line invocation
if __name__ == "__main__": main(sys.argv)

# -----------------------------------------------------------------------------
# EOF
