import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


file_path = "testset.csv" 
df = pd.read_csv(file_path, header=None)

# Manually assign column names
df.columns = ['Amino Acids', 'Mapped Char']

amino_acid_sequences = df['Amino Acids']
mapped_alphabet_sequences = df['Mapped Char']

# Create a DataFrame from the sequences
data = []
for aa_seq, map_seq in zip(amino_acid_sequences, mapped_alphabet_sequences):
    for aa, map_char in zip(aa_seq, map_seq):
        data.append((aa, map_char))

df = pd.DataFrame(data, columns=["Amino Acid", "Mapped Char"])

# table to count the occurrences
pivot_table = df.pivot_table(index="Amino Acid", columns="Mapped Char", aggfunc=len, fill_value=0)

# plot
plt.figure(figsize=(20, 10))
sns.heatmap(pivot_table, annot=False, fmt="d", cmap="viridis")
plt.title("Amino Acid to Mapped Character Count Matrix")
plt.xlabel("Mapped Character")
plt.ylabel("Amino Acid")
plt.show()
