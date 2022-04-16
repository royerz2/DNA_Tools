import SequenceCompilerLibrary as SeqLib
from Bio import SeqIO, Entrez, Restriction
from pydna.assembly import Assembly

def assembly_analysis(assembly_type ,*inserts):

    assembly = Assembly(inserts, limit=14)

    if assembly_type == 'linear':
        assembly_result = assembly.assemble_linear()
    else:
        assembly_result = assembly.assemble_circular()

    return assembly, assembly_result


if __name__ == "__main__":
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id='M77789.2') as handle:
        plasmid = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"
    insert1 =
    insert2 =
    linear_plasmid = plasmid.linearize('XbaI')

    assembly_analysis(linear_plasmid, insert1, insert2)