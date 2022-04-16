import urllib.error
import pandas as pd
from Bio import SeqIO, Entrez
from openpyxl import workbook, load_workbook
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plasmidviewer as pv

nucleotide_id = "M77789.2"
Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=nucleotide_id) as handle:
    record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

feats = []
facecolors = []
edgecolors = []
labelcolors = []

for feat in record.features:
    print(feat)
    if feat.type == "source" or feat.type == "primer_bind":
        pass
    else:
        feats.append(feat)
        if feat.type == "CDS" and feat.strand >= 0:
            facecolors.append("#FFDFDF")
            edgecolors.append("#EE0000")
            labelcolors.append("#DD0000")

        elif feat.type == "CDS" and feat.strand == -1:
            facecolors.append("#DFDFFF")
            edgecolors.append("#0000EE")
            labelcolors.append("#0000DD")
        else:
            edgecolors.append(None)
            labelcolors.append(None)
            facecolors.append(None)


fig, ax = pv.visualize(record, feature_list=feats, edgecolor=edgecolors, title="lentiCas9-EGFP")
fig.savefig("test2.png")
