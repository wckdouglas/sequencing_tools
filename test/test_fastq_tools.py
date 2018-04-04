from sequencing_tools.fastq_tools import *
import six


test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'

def test_complement():
    assert(complement('ACTGNactgn') == 'TGACNtgacn')


def test_reverse_complement():
    assert(reverse_complement('ACTGNactgn') == 'ncagtNCAGT')


def test_fastq():
    fq_iter = readfq(open(test_data_path + '/corrected.qual.fq'))
    record = six.next(fq_iter)
    assert(record.id=="AAAAAA_1_member/1")
    assert(record.seq=='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    assert(record.qual=="EAEE6666EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA")

    record.subseq(0,20)
    assert(record.seq == 'AAAAAAAAAAAAAAAAAAAA')
    assert(record.qual=='EAEE6666EEEEEEEEEEEE')