from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
import SequenceCompilerLibrary as SeqC
from pydna.amplify import pcr

Entrez.email = 'royerz.2002@gmail.com'


def assembly_analysis(plasmid_id='L37382.1',
                      plasmid_digestion='XbaI',
                      insert_ids=[('X71080.1', 'XbaI')]):

    # Fetch and Digest plasmid
    with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=plasmid_id) as handle:
        seq_record = SeqIO.read(handle, 'gb')

    dseq_record = Dseqrecord(seq_record, linear=False, circular=True)
    digested_plasmid = dseq_record.linearize(getattr(Restriction, plasmid_digestion))

    if len(insert_ids) == 1:
        single_insert = True

    amplicons = []

    # if there is only one insert the logic chooses to insert digestion sites to the insert homogeneously
    if single_insert:
        with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=insert_ids[0][0]) as handle:
            seq_record = SeqIO.read(handle, 'gb')

        fw_primer, rev_primer, fw_assembly, rev_assembly = SeqC.design_primers(seq_record.seq,
                                                                               plasmid_digestion,
                                                                               plasmid_digestion)

        pcr_prod = pcr(fw_assembly, rev_assembly, seq_record)
        pcr_prod_dseqrecord = Dseqrecord(pcr_prod.seq, linear=True, circular=False)
        amplicons.append(Seq(pcr_prod.seq))

    # If there is only one insert the logic chooses to insert digestion sites to the insert heterogeneously with next
    # restriction site as the first of the next one.
    else:
        i = 0
        if i == 0:
            fw_primer, rev_primer, fw_assembly, rev_assembly = SeqC.design_primers(seq_record.seq,
                                                                                   plasmid_digestion,
                                                                                   insert_ids[i + 1][1])

            with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=insert_ids[i][0]) as handle:
                seq_record = SeqIO.read(handle, 'gb')

            pcr_prod = pcr(fw_assembly, rev_assembly, seq_record)
            amplicons.append(Seq(pcr_prod.seq))

        elif i == len(insert_ids)-1:
            fw_primer, rev_primer, fw_assembly, rev_assembly = SeqC.design_primers(seq_record.seq,
                                                                                   insert_ids[i + 1][1],
                                                                                   plasmid_digestion)

            with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=insert_ids[i][0]) as handle:
                seq_record = SeqIO.read(handle, 'gb')

            pcr_prod = pcr(fw_assembly, rev_assembly, seq_record)
            amplicons.append(Seq(pcr_prod.seq))

        else:
            fw_primer, rev_primer, fw_assembly, rev_assembly = SeqC.design_primers(seq_record.seq,
                                                                                   insert_ids[i + 1][1],
                                                                                   plasmid_digestion)

            with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=insert_ids[i][0]) as handle:
                seq_record = SeqIO.read(handle, 'gb')

            pcr_prod = pcr(fw_assembly, rev_assembly, seq_record)
            amplicons.append(Seq(pcr_prod.seq))

        amplicon = pcr_prod.seq
        amplicons.append(amplicon)

    # Pass values to Assembly module
    i = 0
    for sequence in amplicons:
        amplicons[i] = Dseqrecord(sequence)

    from pydna.common_sub_strings import terminal_overlap

    assembly = Assembly((*amplicons, digested_plasmid), limit=9)

    print('Doing circular assembly.')
    print(digested_plasmid)
    for item in pcr_prod_dseqrecord.cut(getattr(Restriction, plasmid_digestion)):
        print(item)
    assemblyResult = (digested_plasmid+pcr_prod_dseqrecord.cut(getattr(Restriction, plasmid_digestion))[1]).looped()

    return assemblyResult


if __name__ == "__main__":
    assembly_result = assembly_analysis()
    print(assembly_result)
