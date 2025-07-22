from fasta_reader import read_fasta
from pathlib import Path

# Caminho do arquivo
arquivo = "Ecoli_Sakai_cds_from_genomic.fna"

# Tabela do código genético RNA -> proteína (1 letra)
genetic_code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

def traduzir(seq_rna):
    proteina = ""
    for i in range(0, len(seq_rna) - 2, 3):
        codon = seq_rna[i:i+3]
        aa = genetic_code.get(codon, "")
        proteina += aa
    return proteina

# Carrega todas as sequências do arquivo FASTA
sequencias = list(read_fasta(arquivo))

# 1. Gerar arquivo com sequência de RNA
with open("Sakai_RNA.fasta", "w") as saida_rna:
    for header, dna_seq in sequencias:
        rna_seq = dna_seq.replace("T", "U").upper()
        saida_rna.write(f">{header} RNA\n{rna_seq}\n")

# 2. Gerar 6 arquivos de proteínas para os frames 1–6
for frame in range(6):
    frame_file = f"Frame{frame+1}.fasta"
    with open(frame_file, "w") as out:
        for idx, (header, dna_seq) in enumerate(sequencias):
            dna_seq = dna_seq.upper().replace("T", "U")
            if frame < 3:
                frame_seq = dna_seq[frame:]
            else:
                rev_seq = dna_seq[::-1].translate(str.maketrans("AUCG", "UAGC"))
                frame_seq = rev_seq[frame - 3:]

            proteina = traduzir(frame_seq)
            out.write(f">{header} Frame{frame+1} Proteína\n{proteina}\n")
