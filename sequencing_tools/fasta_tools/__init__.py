import pandas as pd
import numpy as np
from itertools import product
from matplotlib.patches import Rectangle
from ..stats_tools import hamming_distance


class IUPAC:
    """
    Working with IUPAC bases

    Usage::

        iupac = IUPAC()
        iupac.check_degenerate('ACTGN') # {'N'}
        iupac.check_degenerate('ACTGY') # {'Y'}
        list(IUPAC().expand('ACTGN')) #['ACTGA', 'ACTGC', 'ACTGT', 'ACTGG']
        list(IUPAC().expand('ACTGY')) #['ACTGC', 'ACTGT']
    """

    def __init__(self):
        self.IUPAC = {
            "R": ["A", "G"],
            "Y": ["C", "T"],
            "S": ["G", "C"],
            "W": ["A", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "B": ["C", "G", "T"],
            "D": ["A", "G", "T"],
            "H": ["A", "C", "T"],
            "V": ["A", "C", "G"],
            "N": ["A", "C", "T", "G"],
        }

    def check_degenerate(self, seq):
        """
        Args:
            seq (str): sequence with dengerate base

        Returns:
            set: all degenerated bases from this sequence
        """
        bases = set(seq)
        return set(self.IUPAC.keys()).intersection(bases)

    def expand(self, seq):
        """
        output all possible sequence by looping over the dengerative base

        Args:
            seq (str): sequence with dengerate base

        Returns:
            Generator(str): all possible sequences from the DNA/RNA pool
        """
        degenerative_bases = self.check_degenerate(seq)
        expandable_list = [self.IUPAC[b] for b in degenerative_bases]
        for base_combination in product(*expandable_list):
            new_seq = seq
            for i, db in enumerate(degenerative_bases):
                new_seq = new_seq.replace(db, base_combination[i])
            assert set(new_seq) - {"A", "C", "T", "G"} == set()
            yield new_seq


def readfa(file_handle):
    """
    A fasta reader iterator

    Args:
        fp: file handle of a fasta file

    Returns:
        (str, str): sequence id, sequence

    Usage::

        with open('test.fa') as fasta:
            for seq_id, seq in readfq(fasta):
                print(seq_id, seq)
    """
    seqid = ""
    seq = ""
    seq_count = 0
    for line in file_handle:
        if line.startswith(">"):
            seq_count += 1
            if seq_count > 1:
                yield seqid, seq
                seq = ""
                seqid = ""
            seqid = line[1:].strip()
        else:
            seq += line.strip()
    yield seqid, seq


class MultiAlignments:
    def __init__(self, fa_file, RNA=False):
        """
        Plotting multiple-alignment fasta, sequences must be of the same length

        Args:
            fa_file (str): fasta file

        Example::

            # $ cat test.fa
            # >1
            # GGGGAATTAGCTCAAGCGGTAGAGCGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG
            # >2
            # GGGGAATTAGCTCAAGCGGTAGAGCGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG
            # >3
            # GGGGAATTAGCTCAAGCGGTAAAACGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG
            # >4
            # GGGCAACTAGCTCAAGCGGTAAAACGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG
            # >5
            # GCGCAACTAGCTCAAGCGGTAAAACGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGAGG
            # >6
            # GGGGAATTAGTTCAAGCGGTAGAGCGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG

            ma = multi_alignment('test.fa')

            ax = plt.subplot()
            ma.plot(ax = ax)

            ma.concensus()
            #('GGGGAATTAGCTCAAGCGGTAAAACGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG',
            # [array([1.]),
            # array([0.16666667, 0.83333333]),
            # array([1.]),
            # array([0.33333333, 0.66666667]),
            # array([1.]),
            # array([1.]),
            # array([0.33333333, 0.66666667]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([0.83333333, 0.16666667]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([0.5, 0.5]),
            # array([1.]),
            # array([0.5, 0.5]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([1.]),
            # array([0.16666667, 0.83333333]),
            # array([1.])])
        """
        self.records = []
        with open(fa_file) as fa:
            for seqid, seq in readfa(fa):
                if RNA:
                    seq = seq.replace("T", "U").replace("t", "u")
                self.records.append([seqid] + list(seq))

        self.mul_df = pd.DataFrame(self.records).rename(
            columns={0: "seq_id"}
        )  #: sequence matrix, each column is a position, and nucleotide as value
        self.pairwise = None  #: pairwise matrix computed by :py:meth:`sequencing_tools.fasta_tools.MultiAlignment.PairMatrix`

        self.colors = {
            "A": "red",
            "C": "blue",
            "U": "green",
            "G": "orange",
            "-": "black",
            "a": "red",
            "c": "blue",
            "u": "green",
            "g": "orange",
            "t": "green",
        }  #: color dictionary guiding the multiplex alignment plotting

    def plot(
        self,
        ax,
        min_pos=0,
        max_pos=None,
        fontsize=20,
        labelsize=20,
        sample_regex="[A-Za-z0-9_-]+",
    ):
        """
        Args:
            ax (plt.axes): matplotlib axes
            min_pos (float): start position to plot on the multialignments
            max_pos (float): end position to plot on the mutlialignments
            fontsize (int): fontsize for the nucleotides
            labelsize (int): fontsize for sequence id
            sample_regex (str): regex for including sequecne id

        """

        if not max_pos:
            max_pos = self.mul_df.shape[1]
        ax.plot([min_pos, max_pos], [0, self.mul_df.shape[0]], alpha=0)
        for i, (id, seq) in enumerate(
            self.mul_df.pipe(lambda d: d[d.seq_id.str.contains(sample_regex)])
            .set_index("seq_id")
            .iterrows()
        ):
            ax.text(min_pos - 1, i, id, fontsize=labelsize, ha="right", va="center")

            for j, b in seq.items():
                if min_pos <= j <= max_pos:
                    b = "U" if b == "T" else b
                    t = ax.text(
                        j,
                        i,
                        b,
                        fontsize=fontsize,
                        family="monospace",
                        va="center",
                        ha="center",
                        color=self.colors[b],
                    )
        [ax.spines[s].set_visible(False) for s in ["top", "right", "left", "bottom"]]
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.tick_params(axis=u"both", which=u"both", length=0)

    def concensus(self):
        """
        compute a consensus sequence from highest frequency base at each position

        Returns:
            tuple(str, list of numpy array: (consensus sequence, fraction of base making up the composition of the position)
        """
        sequence_matrix = np.array(self.records)[:, 1:]
        consensus_seq = ""
        scores = []
        for pos in range(sequence_matrix.shape[1]):
            bases = sequence_matrix[:, pos]
            # bases = bases[bases!='-']
            b, bcount = np.unique(bases, return_counts=True)
            cb = b[bcount.argmax()]

            # if cb == '-':
            #    cb = str(b[bcount.argmax()][0])
            # elif len(b) > 1:
            #    cb = str(b[bcount == np.sort(bcount)[-2]][0])

            consensus_seq += cb
            score = bcount / bcount.sum()
            scores.append(score)
        return consensus_seq, scores

    def PairMatrix(self):
        """
        Calculate the hamming distances between each sequence pair
        """
        pairwise = []
        for id1, seq1 in self.mul_df.set_index("seq_id").iterrows():
            for id2, seq2 in self.mul_df.set_index("seq_id").iterrows():
                seq1 = "".join(seq1)
                seq2 = "".join(seq2)
                record = (id1, id2, hamming_distance(seq1, seq2))
                pairwise.append(record)

        self.pairwise = pd.DataFrame(pairwise, columns=["id1", "id2", "distance"])
