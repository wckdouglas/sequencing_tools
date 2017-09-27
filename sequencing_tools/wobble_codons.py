#!/usr/bin/env python


from sequencing_tools.fastq_tools import  reverse_complement

class wobble_codon:
    def __init__(self):
        '''
        From: Solving the riddle of codon usage preferences: a test for translational selection
        tRNA_table: key is anticodon,
                    values: list of codon recognized 
        '''
        self.codon_dict = None
        self.anticodon_dict = { # AA: F
                            'AAA': ['TTT'],
                            'GAA': ['TTT','TTC'],

                            #AA: L
                            'TAA': ['TTA','TTG'],
                            'CAA': ['TTG'],
                            'AAG': ['CTT','CTC', 'CTA'],
                            'GAG': ['CTC','CTT'],
                            'TAG': ['CTA','CTG'],
                            'CAG': ['CTG'],


                            # AA: S
                            'AGA': ['TCT','TCC','TCA'],
                            'GGA': ['TCT','TCC'],
                            'TGA': ['TCA','TCG'],
                            'CGA': ['TCG'],

                            # AA: Y
                            'ATA': ['TAT'],
                            'GTA': ['TAC','TAT'],

                            # AA: *
                            'TTA': ['TAA'],
                            'CTA': ['TAG'],
                            'TCA': ['TGA'],

                            #AA: C
                            'ACA': ['TGT'],
                            'GCA': ['TGT','TGC'],

                            #AA: W
                            'CCA': ['TGG'],

                            #AA: P
                            'AGG': ['CCT', 'CCC','CCA'],
                            'GGG': ['CCT','CCC'],
                            'TGG': ['CCA','CCG'],
                            'CGG': ['CCG'],

                            #AA: H
                            'ATG': ['CAT'],
                            'GTG': ['CAC','CAT'],

                            #AA: Q
                            'TTG': ['CAA','CAG'],
                            'CTG': ['CAG'],

                            #AA: R
                            'ACG': ['CGT','CGC','CGA'],
                            'GCG': ['CGC','CGT'],
                            'TCG': ['CGA','CGG'],
                            'CCG': ['CGG'],

                            #AA: I
                            'AAT': ['ATT','ATC','ATA'],
                            'GAT': ['ATT','ATC','ATA'],
                            'TAT': ['ATA'],

                            #AA: M
                            'CAT': ['ATG'],

                            #AA: T 
                            'AGT': ['ACT','ACC','ACA'],
                            'GGT': ['ACT','ACC'],
                            'TGT': ['ACA','ACG'],
                            'CGT': ['ACG'],

                            #AA: N
                            'ATT': ['AAT'],
                            'GTT': ['AAT','AAC'],

                            #AA: K
                            'TTT': ['AAA','AAG'],
                            'CTT': ['AAG'],

                            #AA: S
                            'ACT': ['AGT'],
                            'GCT': ['AGT','AGC'],

                            #AA: R
                            'TCT': ['AGA','AGG'],
                            'CCT': ['AGG'],

                            #AA: V
                            'AAC': ['GTT','GTC','GTA'],
                            'GAC': ['GTC','GTT'],
                            'TAC': ['GTA','GTG'],
                            'CAC':['GTG'],

                            #AA: A
                            'AGC': ['GCT','GCC','GCA'],
                            'GGC': ['GCT','GCC'],
                            'TGC': ['GCA','GCG'],
                            'CGC': ['GCG'],

                            #AA: D
                            'ATC': ['GAT'],
                            'GTC': ['GAC','GAT'],

                            #AA: E
                            'TTC': ['GAA'],
                            'CTC':['GAG'],

                            #AA: G
                            'ACC': ['GGT','GGC','GGA'],
                            'GCC': ['GGT','GGC'],
                            'TCC': ['GGA', 'GGG'],
                            'CCC': ['GGG']
                            }



    def check_codon(self):
        for anticodon, codons in self.anticodon_dict.iteritems():
            assert reverse_complement(anticodon) in codons, '%s missing codon' %(anticodon)
            matching_codon = reverse_complement(anticodon)
            for codon in codons:
                assert codon[:2] == matching_codon[:2], 'Anticodon: %s matching codon: %s?' %(anticodon,codon)
        print 'All clear'

    def make_codon_dict(self):
        self.codon_dict = {}
        for anticodon, codons in self.anticodon_dict.iteritems():
            for codon in codons:
                codon_dict[codon] = anticodon


class tRNA_adaptation_index:
    def __init__(self, tRNA_count):
    '''
    adapted from https://github.com/smsaladi/tAI/blob/master/tAI/tAI.py

    input:
    * tRNA_count: dictionary with anticodon as keys, raw count as values
    '''

        self.tRNA_dict = tRNA_count
        self.wobble = wobble_codon()
        self.wobble.make_codon_dict()
        self.wobble_dict = self.wobble.codon_dict
        self.p = {'T': 0.59, 'C': 0.72, 'A': 0.0001, 'G': 0.32}
        self.isoleucine_p = 1 - 0.89

        tRNA_availability = defaultdict(int)
        for codon in wobble_dict.iteritems():
            codon_prefix = codon[:2]
            wobble_base = codon[2]

            if wobble_base == 'T':                  # INN -> NNT, NNC, NNA
                weights[codon] = self.tRNA_dict[reverse_complement(codon)] + self.p['T']*self.tRNA_dict[reverse_complement(codon_prefix+'C')]
            elif wobble_base == 'C':                # GNN -> NNT, NNC
                weights[codon] = self.tRNA_dict[reverse_complement(codon)] + self.p['C']*self.tRNA_dict[reverse_complement(codon_prefix+'T')]
            elif wobble_base == 'A':                # TNN -> NNA, NNG
                weights[codon] = self.tRNA_dict[reverse_complement(codon)] + self.p['A']*self.tRNA_dict[reverse_complement(codon_prefix+'T')]
            elif wobble_base == 'G':                # CNN -> NNG
                weights[codon] = self.tRNA_dict[reverse_complement(codon)] + self.p['G']*self.tRNA_dict[reverse_complement(codon_prefix+'A')]

            total_tRNA = sum(tRNA_dict[anticodon] for anticodon in anticodons)
            tRNA_availability[codon] = total_tRNA


   for codon in list(weights.index):
        base = codon[:2]

            raise ValueError('Non-standard codon or notation')