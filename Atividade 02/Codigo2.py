from Bio import SeqIO
import csv
import os

# Nomes corretos dos seus arquivos
arquivos_fasta = {
    "Genoma_completo": "Ecoli Sakai sequence.fasta",
    "Plasmidio_1": "Ecoli Sakai Plasmid 1 sequence.fasta",
    "Plasmidio_2": "Ecoli Sakai plasmid 2 sequence.fasta"
}

def contar_nucleotideos(sequencia):
    sequencia = sequencia.upper()
    a = sequencia.count("A")
    t = sequencia.count("T")
    c = sequencia.count("C")
    g = sequencia.count("G")
    total = a + t + c + g
    return a, c, g, t, total

for nome, caminho in arquivos_fasta.items():
    if not os.path.exists(caminho):
        print(f"Arquivo {caminho} não encontrado.")
        continue

    for record in SeqIO.parse(caminho, "fasta"):
        seq = str(record.seq)
        a, c, g, t, total = contar_nucleotideos(seq)

        nome_saida = f"{nome}_contagem.csv"
        with open(nome_saida, mode="w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow([f"Dados do {nome.replace('_', ' ')}"])
            writer.writerow(["A", "C", "G", "T", "Total"])
            writer.writerow([a, c, g, t, total])

        print(f"[✔] Arquivo criado: {nome_saida}")
