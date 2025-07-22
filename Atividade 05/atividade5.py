from Bio import SeqIO
from Bio.Seq import Seq
import os

# Número do RA
RA = "164931" 

# Dicionário do código genético padrão
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

# Função para traduzir uma sequência de DNA para proteína (1 letra)
def translate_dna(seq):
    protein = ""
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            protein += codon_table.get(codon, 'X')  # X para codon desconhecido
    return protein

# Função para encontrar regiões codificantes de proteínas
def find_proteins(dna_seq, frame_offset, reverse=False):
    proteins = []
    start_codon = 'M'
    stop_codon = '*'
    seq_len = len(dna_seq)
    frame = dna_seq[frame_offset:]
    protein_seq = translate_dna(frame)
    pos = 0
    protein_index = 1

    while pos < len(protein_seq):
        if protein_seq[pos] == start_codon:
            start = pos
            while pos < len(protein_seq) and protein_seq[pos] != stop_codon:
                pos += 1
            if pos < len(protein_seq):
                prot = protein_seq[start:pos]
                if len(prot) >= 30:  # mínimo de 30 aminoácidos
                    nt_start = frame_offset + start*3
                    nt_end = frame_offset + pos*3 + 3
                    if reverse:
                        nt_start = seq_len - nt_start
                        nt_end = seq_len - (frame_offset + pos*3 + 3)
                        nt_start, nt_end = nt_end, nt_start  # ajustar ordem
                    proteins.append((f"proteína {protein_index}", prot, nt_start+1, nt_end))
                    protein_index += 1
        pos += 1
    return proteins

# Processar cada sequência
def process_file(filename, prefix, genbank_code):
    record = SeqIO.read(filename, "fasta")
    dna = str(record.seq).upper()
    rev_dna = str(record.seq.reverse_complement()).upper()

    for i in range(3):  # frames 1, 2, 3
        proteins = find_proteins(dna, i)
        with open(f"{prefix}_frame{i+1}_ativ5_{RA}.fasta", "w") as f:
            for idx, (name, seq, start, end) in enumerate(proteins, 1):
                f.write(f">{genbank_code}, Frame {i+1}, {name}, [location={start}..{end}], {RA}\n")
                f.write(f"{seq}\n")

    for i in range(3):  # frames 4, 5, 6
        proteins = find_proteins(rev_dna, i, reverse=True)
        with open(f"{prefix}_frame{i+4}_ativ5_{RA}.fasta", "w") as f:
            for idx, (name, seq, start, end) in enumerate(proteins, 1):
                f.write(f">{genbank_code}, Frame {i+4}, {name}, [location={start}..{end}], {RA}\n")
                f.write(f"{seq}\n")

# Processar os três genomas
process_file("Ecoli Sakai sequence.fasta", "gc", "BA000007.3")
process_file("Ecoli Sakai Plasmid 1 sequence.fasta", "p1", "AB011549.2")
process_file("Ecoli Sakai Plasmid 2 sequence.fasta", "p2", "AB011548.2")
