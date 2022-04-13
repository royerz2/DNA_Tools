from Bio import SeqIO, Entrez

Entrez.email='royerz.2002@gmail.com'


gb_file = "M77789"
for index, record in enumerate(SeqIO.parse("M77789", "gb")):
    print(
        "index %i, ID = %s, length %i, with %i features"
        % (index, record.id, len(record.seq), len(record.features))
    )