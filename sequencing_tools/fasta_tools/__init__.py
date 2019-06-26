from matplotlib.patches import Rectangle
from sequencing_tools.stats_tools import hamming_distance
import pandas as pd


def readfa(file_handle):
    seqid = ''
    seq = ''
    seq_count = 0
    for line in file_handle:
        if line.startswith('>'):
            seq_count += 1
            if seq_count > 1:
                yield seqid, seq
                seq = ''
                seqid = ''
            seqid = line[1:].strip()
        else:
            seq += line.strip()
    yield seqid, seq
       

class MutliAlignments():
    '''
    plotting multiple-alignment fasta
    '''
    def __init__(self, fa_file, RNA=False):
        '''
        input:
            fa_file: fasta file
        
        usage:
            ma = multi_alignment("fasta file")
            ma.plot(ax = ax)
        '''
        
        self.records = []
        with open(fa_file) as fa:
            for seqid, seq in readfa(fa):
                if RNA:
                    seq = seq.replace('T','U').replace('t','u')
                self.records.append([seqid] + list(seq))

        
        self.mul_df = pd.DataFrame(self.records)\
            .rename(columns = {0:'seq_id'})
        self.pairwise = None
        self.colors = {'A':'red','C':'blue','U':'green','G': 'orange', '-': 'black',
                      'a':'red','c':'blue','u':'green','g':'orange','t':'green'}

        
    def plot(self, ax, 
             min_pos=0, 
             max_pos=None, 
             fontsize=20,
             labelsize=20,
             sample_regex = '[A-Za-z0-9_-]+'):
        '''
        input:
            ax: matplotlib axes
            min_pos: start position to plot on the multialignments
            max_pos: end position to plot on the mutlialignments
            fontsize: fontsize for the nucleotides
            labelsize: fontsize for sequence id
            sample_regex: regex for including sequecne id
        
        usage:
            ma = multi_alignment("fasta file")
            ma.plot()
        '''

        if not max_pos:
            max_pos = self.mul_df.shape[1]
        ax.plot([min_pos, max_pos], 
                [0, self.mul_df.shape[0]], 
                alpha=0)
        for i, (id, seq) in enumerate(self.mul_df\
                                .pipe(lambda d: d[d.seq_id.str.contains(sample_regex)])\
                                .set_index('seq_id') \
                                .iterrows()):
            ax.text(min_pos - 1, i, id, 
                    fontsize=labelsize, 
                    ha = 'right', 
                    va = 'center')

            for j, b in seq.items():
                if min_pos <= j <= max_pos:
                    b = 'U' if b == 'T' else b
                    t = ax.text(j,i,b, fontsize=fontsize, 
                            family='monospace', 
                            va = 'center', ha='center',
                            color = self.colors[b])
        [ax.spines[s].set_visible(False) for s in ['top','right','left', 'bottom']]
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.tick_params(axis=u'both', which=u'both',length=0)



    def concensus(self):
        '''
        compute a concensus sequence
        '''
        sequence_matrix = np.array(self.records)[:, 1:]
        concensus_seq = ''
        scores = []
        for pos in range(sequence_matrix.shape[1]):
            bases = sequence_matrix[:,pos]
            #bases = bases[bases!='-']
            b, bcount = np.unique(bases, return_counts=True)
            cb = b[bcount.argmax()]

            #if cb == '-':
            #    cb = str(b[bcount.argmax()][0])
            #elif len(b) > 1:
            #    cb = str(b[bcount == np.sort(bcount)[-2]][0])

            concensus_seq += cb
            score = bcount/bcount.sum()
            scores.append(score)
        return concensus_seq, scores


    def PairMatrix(self):
        pairwise = []
        for id1, seq1 in self.mul_df\
                       .set_index('seq_id')\
                       .iterrows():
            for id2, seq2 in self.mul_df\
                           .set_index('seq_id')\
                           .iterrows():
                seq1 = ''.join(seq1)
                seq2 = ''.join(seq2)
                record = (id1, id2, hamming_distance(seq1, seq2))
                pairwise.append(record)

        self.pairwise = pd.DataFrame(pairwise, columns = ['id1','id2','distance'])

