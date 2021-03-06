import os
from urllib.request import Request, urlopen
import gzip
from sequencing_tools.gene_tools import transcriptome

DataDir = os.path.dirname(os.path.abspath(__file__))


def test_transcriptome():
    refFlat = DataDir + "/data/test.refFlat"
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz"
    if not os.path.isfile(refFlat):
        req = Request(url)
        response = urlopen(req)
        with open(refFlat, "w") as f:
            for line in gzip.decompress(response.read()).decode().split("\n")[:100]:
                print(line, file=f)

    txome = transcriptome.Transcriptome(refFlat, coding_only=False)
    gene = txome.get_gene("WASH7P")
    tx = gene["NR_024540"]
    assert tx.exon_count == 11
    assert tx.blocks(1, 100) == [(29320, 29369), (24841, 24891)]
