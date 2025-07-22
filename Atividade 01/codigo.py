# Abrir o arquivo FASTA para leitura
with open('teste1.fasta', 'r') as arquivo:
    linhas = arquivo.readlines()

# Pegar o nome do gene (primeira linha do FASTA)
nome_gene = linhas[0].strip()

# Juntar as linhas seguintes e tirar espaços em branco/linhas novas
sequencia = ''.join(linhas[1:]).replace('\n', '').upper()

# Contar cada nucleotídeo
num_A = sequencia.count('A')
num_T = sequencia.count('T')
num_C = sequencia.count('C')
num_G = sequencia.count('G')

# Calcular total
total = num_A + num_T + num_C + num_G

# Gerar arquivo de saída
with open('resultado.txt', 'w') as saida:
    saida.write(f'Nome do gene: {nome_gene}\n')
    saida.write(f'A: {num_A}\n')
    saida.write(f'C: {num_C}\n')
    saida.write(f'G: {num_G}\n')
    saida.write(f'T: {num_T}\n')
    saida.write(f'Total: {total} nucleotídeos\n')

print("Arquivo resultado.txt gerado com sucesso!")
