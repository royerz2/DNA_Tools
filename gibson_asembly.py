from Bio import SeqIO, Entrez, Restriction, SeqRecord
from Bio.Seq import Seq
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
import SequenceCompilerLibrary as SeqC
from pydna.amplify import pcr

Entrez.email = 'royerz.2002@gmail.com'


def prepare_insert(insert_id, insert_type, insert_restriction_5, insert_restriction_3):
    if insert_type == 'nucleotide':
        newInsert = SeqC.entrez_fetch_nucleotide(insert_id, insert_restriction_5, insert_restriction_3)

    elif insert_type == 'protein':
        newInsert = SeqC.entrez_fetch_protein(insert_id, insert_restriction_5, insert_restriction_3)

    else:
        raise IndexError('This function can only accept nucleotide or protein type inserts.')

    # Run PCR with the assembly primers to introduce restriction sites we need.
    pcr_prod = pcr(newInsert.fw_assembly, newInsert.rev_assembly, newInsert.dna)

    # Convert pcr product sequence to Dseqrecord object.
    pcr_prod_dseqrecord = Dseqrecord(pcr_prod.seq, linear=True, circular=False)

    if insert_restriction_3 == insert_restriction_5:
        # Cut 5' and 3' ends.
        restricted_pcr_prod_dseqrecords = pcr_prod_dseqrecord.cut(getattr(Restriction, insert_restriction_5))
        restricted_pcr_prod_dseqrecord = restricted_pcr_prod_dseqrecords[1]

    else:
        # Cut 5' end.
        five_restricted_inserts = pcr_prod_dseqrecord.cut(getattr(Restriction, insert_restriction_5))
        five_restricted_insert = five_restricted_inserts[1]

        # Cut 3' end.
        all_restricted_inserts = five_restricted_insert.cut(getattr(Restriction, insert_restriction_3))
        restricted_pcr_prod_dseqrecord = all_restricted_inserts[1]

    return restricted_pcr_prod_dseqrecord


def assembly_analysis(plasmid_id='L37382.1',
                      plasmid_digestion='HindIII',
                      insert_ids=[('X71080.1', 'HindIII', 'nucleotide'),
                                  ('X71080.1', 'XbaI', 'protein')]):

    # Fetch and Digest plasmid
    with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=plasmid_id) as handle:
        seq_record = SeqIO.read(handle, 'gb')

    dseq_record = Dseqrecord(seq_record, linear=False, circular=True)
    restricted_plasmid = dseq_record.linearize(getattr(Restriction, plasmid_digestion))

    if len(insert_ids) == 1:
        single_insert = True
    else:
        single_insert = False

    # if there is only one insert the logic chooses to insert digestion sites to the insert homogeneously.
    if single_insert:
        insert_id = insert_ids[0][0]
        insert_type = insert_ids[0][2]
        insert_restriction_5 = plasmid_digestion
        insert_restriction_3 = plasmid_digestion

        restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                        insert_type,
                                                        insert_restriction_5,
                                                        insert_restriction_3)

        assemblyResult = (restricted_plasmid + restricted_pcr_prod_dseqrecord).looped()

    # If there is only one insert the logic chooses to insert digestion sites to the insert heterogeneously with next
    # restriction site as the first of the next one.
    else:
        ampliconsList = []
        restricted_ampliconsList = []
        linear_assembly = restricted_plasmid
        i = 0
        # If in the first insert position, digest 5' with plasmid restricting enzyme, 3' with tuple spesifcied enzyme.
        if i == 0:
            insert_id = insert_ids[0][0]
            insert_type = insert_ids[0][2]
            insert_restriction_5 = plasmid_digestion
            insert_restriction_3 = insert_ids[i + 1][1]

            restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                            insert_type,
                                                            insert_restriction_5,
                                                            insert_restriction_3)

            # Add results to the respective lists.
            restricted_ampliconsList.append(restricted_pcr_prod_dseqrecord)
            print(restricted_plasmid, restricted_pcr_prod_dseqrecord)
            linear_assembly = restricted_plasmid + restricted_pcr_prod_dseqrecord

        # Elif the last insert position, digest 5' with tuple spesifcied enzyme, 3' with plasmid restricting enzyme.
        elif i == len(insert_ids)-1:
            insert_id = insert_ids[0][0]
            insert_type = insert_ids[0][2]
            insert_restriction_5 = insert_ids[i][1]
            insert_restriction_3 = plasmid_digestion

            restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                            insert_type,
                                                            insert_restriction_5,
                                                            insert_restriction_3)

            restricted_ampliconsList.append(restricted_pcr_prod_dseqrecord)
            linear_assembly = restricted_plasmid + restricted_pcr_prod_dseqrecord

        # Else in a mid insert position, digest 5' and 3' with tuple spesifcied enzyme.
        else:
            insert_id = insert_ids[0][0]
            insert_type = insert_ids[0][2]
            insert_restriction_5 = insert_ids[i][1]
            insert_restriction_3 = insert_ids[i + 1][1]
            restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                            insert_type,
                                                            insert_restriction_5,
                                                            insert_restriction_3)

            restricted_ampliconsList.append(restricted_pcr_prod_dseqrecord)
            linear_assembly = restricted_plasmid + restricted_pcr_prod_dseqrecord

        assemblyResult = linear_assembly.looped()

    print('Doing circular assembly.')
    print(restricted_plasmid)

    return assemblyResult


if __name__ == "__main__":
    assembly_result = assembly_analysis(insert_ids=[('CAA50398.1', 'HindIII', 'protein')])
    print(assembly_result)
    print(str(assembly_result.seq))
