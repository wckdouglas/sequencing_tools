from sequencing_tools.gene_tools.transcriptome import Transcriptome
from pathlib import Path
import pytest
from sequencing_tools.utils import SeqUtilsError
from more_itertools import only
from hypothesis import given, strategies as st

test_data_path = Path(__file__).parent.resolve() / "data"


@pytest.fixture(scope="module")
def transcriptome():
    return Transcriptome(test_data_path / "test_refflat.txt")

@pytest.fixture(scope="module")
def transcript(transcriptome):
    gene = transcriptome.get_gene("WASHC2C")
    return only(gene.values())


def test_Transcriptome_coding_only():
    transcriptome = Transcriptome(test_data_path / "test_refflat.txt", coding_only=True)
    assert len(transcriptome.transcript_dict) == 2


def test_Transcriptome_bad_gene_name(transcriptome):
    assert len(transcriptome.transcript_dict) == 2
    with pytest.raises(SeqUtilsError) as e:
        transcriptome.get_gene("not_a_gene")
    assert "not_a_gene not in database" in str(e.value), "Didn't capture bad gene name"


def test_transcript(transcript):
    assert transcript is not None, "Didn't get gene WASHC2C"
    assert transcript.transcript_length == 4535


@given(total_block_size=st.integers(min_value=1, max_value=3000),
    tstart=st.integers(min_value=1, max_value=1000))
def test_Transcript_blocks(transcript, total_block_size, tstart):
    block_size = 0
    for block_start, block_end in transcript.blocks(tstart, total_block_size+tstart):
        block_size += (block_end - block_start)
    assert block_size == total_block_size, "Didn't get correct block size"


def test_Exon(transcript):
    exon = transcript.exons[1]
    assert exon.length == 55
    assert exon.contain_cds == 0

    assert transcript.exons[2].after_cds == 0
    assert transcript.exons[3].contain_cds == 1
    assert transcript.exons[4].after_cds == 1

    assert transcript.exons[29].after_cde == 0
    assert transcript.exons[30].contain_cde == 1