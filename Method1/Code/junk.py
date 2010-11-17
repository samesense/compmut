#In true genomic fashion, here is a large region of junk
#I believe that it was originally used to create the random tree-pair
#distribution that our distribution was superimposed on (the red curve in
#the graph) ~Greg 01/04/06
#Construct random pairs of positions

incrementor, logs = 0, []
while incrementor < 1000:
    instances = []
    print incrementor
    i = 0
    while i < len(motifs):
        # Pick a gene from the alignment
        gene_key = random.sample(alignment, 1)[0]
        sequence = alignment[gene_key]["human"]
        # Pick a working 2-mer motif from the gene
        check = 0
        flag = False
        while True:
            position = int(random.uniform(0, len(sequence)-2))
            if verify(alignment, gene_key, position): break
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
    logs.append(math.log(float(mega[0]["CC"])/exp))
    
    incrementor += 1

f = open("random_logs.txt", "w")
pickle.dump(logs, f)
f.close()
for i in range(len(logs)):
    print i, logs[i]
sys.exit(0)
