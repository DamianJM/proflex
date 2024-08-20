import matplotlib.pyplot as plt
import seaborn as sns

# File containing the values
values_file = "outputsteps.txt"

with open(values_file, 'r') as file:
    values = [float(line.strip()) for line in file]

# Plot the density plot
sns.kdeplot(values, shade=True)
plt.title('Density Plot of Values')
plt.xlabel('Sequence Length')
plt.ylabel('Density')
plt.grid(False)
plt.show()
