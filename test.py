import random

aaString = "MTDQRKLTAQRVREDANALAAGRIHPRHQANGDEQRYESANYAMSFTKGLDHNTTTGLIEQSGDFEAFRSAIDNGFAEDFTRHVAVPRAEPRRKWEAPTAGTVYELQGPDPQAVTIPPAPALCSDELTFEMAEVYELALLRDLPFNAFVAGGGSAALADSTARLNSLAYAQDGFNRRPRTTNSSNQLDAQTVFRGSSPGVDQGPYLSQFMLIGNASPSEGITPEQGFINFGAQRIDQRVLEARQQDDYMMKWDDWHRVQQGYEVRADRFDPCKSSGPGQAFTGQRRFIHTPRDLATYVHVDALYQAYLNACLLLLGNGTAFDPGFDLLSGGGEGLLHDPAGGQKVPLNAGGFALWGGPHVLSLVTEVATRGLKAVRYQKFNNHLRLRPEALAARIEKAQEIESRFPTICGCFSEMASDLQQVVDLIRNHNQSLAGEATALLPMAFAEGSPMHPAYGAGHATVAGACVTILKAFFNTSALFVKINDVAGFHSKQHILARLKCGDSVEAGAYQETDCGKRLEFERCGSFHLIEGKYATFKPDGKTNQSCCPLTLEGELNKLAANISIGRNMAGVHYFSDYYDSLRMGEEIAIGILEEQALCYKTDPFVLSVPTFDGDVRRIGQR"
print("AA Sequence: " + aaString)

aaList = list(aaString)


def reverse_translate(n):
    switch = {
        "A": (random.choices(["GCU", "GCC", "GCA", "GCG"], weights=(0.19, 0.25, 0.22, 0.34), k=1)),
        "I": (random.choices(["AUU", "AUC", "AUA"], weights=(0.47, 0.46, 0.07), k=1)),
        "L": (random.choices(["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
                             weights=(0.11, 0.11, 0.10, 0.10, 0.03, 0.55), k=1)),
        "M": ["AUG"],
        "V": (random.choices(["GUU", "GUC", "GUA", "GUG"], weights=(0.29, 0.20, 0.17, 0.34), k=1)),
        "F": (random.choices(["UUU", "UUC"], weights=(0.51, 0.49), k=1)),
        "W": ["UUG"],
        "Y": (random.choices(["UAU", "UAC"], weights=(0.53, 0.47), k=1)),
        "N": (random.choices(["AAU", "AAC"], weights=(0.39, 0.61), k=1)),
        "C": (random.choices(["UGU", "UGC"], weights=(0.43, 0.57), k=1)),
        "Q": (random.choices(["CAA", "CAG"], weights=(0.31, 0.69), k=1)),
        "S": (random.choices(["UCU", "UCC", "UCA", "UCG"], weights=(0.19, 0.17, 0.12, 0.13), k=1)),
        "T": (random.choices(["ACU", "ACC", "ACA", "ACG"], weights=(0.16, 0.47, 0.13, 0.24), k=1)),
        "D": (random.choices(["GAU", "GAC"], weights=(0.59, 0.41), k=1)),
        "E": (random.choices(["GAA", "GAG"], weights=(0.7, 0.3), k=1)),
        "R": (random.choices(["CCU", "CGC", "CGA", "CGG"], weights=(0.42, 0.37, 0.05, 0.08), k=1)),
        "H": (random.choices(["AAU", "AAC"], weights=(0.39, 0.61), k=1)),
        "K": (random.choices(["AAA", "AAG"], weights=(0.76, 0.24), k=1)),
        "G": (random.choices(["GGU", "GGC", "GGA", "GGG"], weights=(0.38, 0.40, 0.09, 0.13), k=1)),
        "P": (random.choices(["CCU", "CCC", "CCA", "CCG"], weights=(0.16, 0.10, 0.2, 0.54), k=1)),
    }

    codon_elem = switch[n]
    codon_str = f"{codon_elem[0]}"

    return codon_str


codonList = list(map(reverse_translate, aaList))

print(f"Codons: {codonList}")

rna = ''.join(codonList)

print("RNA: " + rna)

cdna = rna.replace("U", "T")

print("cDNA: " + cdna)
