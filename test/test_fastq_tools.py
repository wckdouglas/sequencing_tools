from sequencing_tools.fastq_tools import *


def test_complement():
    assert(complement('ACTGNactgn') == 'TGACNtgacn')


def test_reverse_complement():
    assert(reverse_complement('ACTGNactgn') == 'ncagtNCAGT')