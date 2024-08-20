import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_path = "matched_sequences.csv"  
df = pd.read_csv(file_path, header=None)


df.columns = ['Amino Acids', 'Mapped Char']

amino_acid_sequences = df['Mapped Char']
amino_acid_counts = {}
total_count = 0

for sequence in amino_acid_sequences:
    for amino_acid in sequence:
        if amino_acid in amino_acid_counts:
            amino_acid_counts[amino_acid] += 1
        else:
            amino_acid_counts[amino_acid] = 1
        total_count += 1

amino_acid_percentages = {aa: (count / total_count) * 100 for aa, count in amino_acid_counts.items()}


percentages_df = pd.DataFrame(list(amino_acid_percentages.items()), columns=['Amino Acid', 'Percentage'])

plt.figure(figsize=(12, 6))
sns.barplot(x='Amino Acid', y='Percentage', data=percentages_df.sort_values('Percentage', ascending=False))
plt.title('Percentage Abundance ProFlex')
plt.xlabel('ProFlex')
plt.ylabel('Percentage')
plt.show()
