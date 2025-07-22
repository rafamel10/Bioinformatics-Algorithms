from Bio import SeqIO
import matplotlib.pyplot as plt
import csv
import math

# Caminho do arquivo FASTA
arquivo_fasta = "Ecoli_Sakai_cds_from_genomic.fna"

# Constantes para cálculo do Tm
Na_conc = 100  # mM
R = 1.987  # cal/(K*mol), constante universal dos gases
Tm_results = []
GC_results = []
dados_nt = []

def calcular_tm(seq, gc_content):
    """Calcula a temperatura de melting (Tm) baseada no conteúdo GC"""
    seq_len = len(seq)
    if seq_len < 1:
        return 0
    # Fórmula com base na referência BioMath
    Tm = 81.5 + 16.6 * math.log10(Na_conc) + 0.41 * (gc_content * 100) - (675 / seq_len)
    return Tm

# Processa todas as sequências do arquivo FASTA
for record in SeqIO.parse(arquivo_fasta, "fasta"):
    seq = str(record.seq).upper()
    A = seq.count('A')
    T = seq.count('T')
    C = seq.count('C')
    G = seq.count('G')
    total = A + T + C + G

    # GC content
    gc_content = (G + C) / total if total > 0 else 0

    # Temperatura de melting
    tm = calcular_tm(seq, gc_content)

    # Armazena os resultados
    dados_nt.append([record.id, A, T, C, G, total])
    GC_results.append([record.id, round(gc_content, 3)])
    Tm_results.append([round(tm, 3), round(gc_content * 100, 2)])  # em porcentagem para gráfico

# D.1: Salvar contagem de nucleotídeos
with open("Dados das sequencias.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["ID", "A", "T", "C", "G", "Total"])
    writer.writerows(dados_nt)

# D.2: Salvar conteúdo GC
with open("Conteudo_GC.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["ID", "GC_Content"])
    writer.writerows(GC_results)

# D.3: Temperatura x GC (%)
with open("Temperatura _x_GC.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Tm (°C)", "GC (%)"])
    writer.writerows(Tm_results)

# D.4: Gráfico GC x Tm
x_tm = [item[0] for item in Tm_results]
y_gc = [item[1] for item in Tm_results]

plt.figure(figsize=(10, 6))
plt.scatter(x_tm, y_gc, alpha=0.5, s=10)
plt.xlabel("Temperatura de Melting (°C)")
plt.ylabel("Conteúdo GC (%)")
plt.title("GC vs. Temperatura de Melting")
plt.suptitle("Gráfico por Rafael Barros", fontsize=10)
plt.grid(True)
plt.savefig("Grafico_GC_vs_Tm.png", dpi=300)
plt.show()
