"""Take PWM gff hits, and make aligned PWM hits. This was done b/c Rithun needs hits for the controls."""
import cm_scan3, pickle, os

gff_file = '../Doc/fasta/human1ku-corrected.fasta.MA0001random.gff'
alignment_file = '../Doc/pickle/alignment.pickle'
with open(alignment_file) as f:
    alignment = pickle.load(f)

species = ('human', 'chimp', 'mouse', 'rat', 'dog')
pval = '-9'
print 'ControlType\tPWM\tPWMhit\tSpecies\tSeq'
for gff_file in os.listdir('../Doc/fasta/'):
    pwm_type = ''
    if 'shuffle' in gff_file and pval in gff_file:
        pwm_type = 'shuffle'
    elif 'random' in gff_file and pval in gff_file:
        pwm_type = 'random'
    if pwm_type:
        with open('../Doc/fasta/' + gff_file) as f:
            motif_hits = cm_scan3.buildMotifs(f, alignment)
            for idx, hit in enumerate(motif_hits):
                for org, seq in zip(species, hit):
                    pwm = gff_file.split(pwm_type)[0]
                    print pwm_type + '\t' + pwm + '\t' + str(idx+1) + '\t' + org + '\t' + seq
            
