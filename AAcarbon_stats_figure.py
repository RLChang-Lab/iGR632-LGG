########## Statistics, Plotting, and Comparision of in silico and in vitro Amino Acid Carbon Source Results ##########
#   After running both the in silico hypothesis generation and in vitro experiments, 
#   Run statistics and plot the results for the paper 
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_1samp


IV_data = {
    "Ser": [0.155, 0.065, 0.041, 0.149, 0.042, 0.138],
    "Phe": [0.203, 0.165, 0.160],
    "Trp": [0.350, 0.045, 0.041, 0.374],
    "Gly": [0.044, 0.042, 0.170, 0.087, 0.010],
    "Lys": [1.301, 0.933, 0.970, 1.120]
}

IS_data={"Ser": 1.128,
         "Phe": 1.03,
         "Trp": 1.01,
         "Gly": 1.11,
         "Lys": 1.00}



# Compute averages and stats
components = list(IV_data.keys())
IV_means = [np.mean(IV_data[c]) for c in components]
IV_std = [np.std(IV_data[c]) for c in components]

# Statistical test (one-sample t-test against 1)
p_values = [ttest_1samp(IV_data[c], 1).pvalue for c in components]

fig, ax = plt.subplots(figsize=(8,6))

bars = ax.bar(components, IV_means, yerr=IV_std, capsize=5, alpha=0.7, color='skyblue', label='IV average')

for i, c in enumerate(components):
    y = IV_data[c]
    x = np.random.normal(i, 0.05, size=len(y))  # jitter for visibility
    ax.scatter(x, y, color='black', zorder=10)

for i, c in enumerate(components):
    ax.scatter(i, IS_data[c], color='red', marker='*', s=200, zorder=15, label='IS FC' if i==0 else "")

for i, p in enumerate(p_values):
    if p < 0.05:
        ax.text(i, max(IV_data[components[i]])+0.1, "*", ha='center', color='black', fontsize=16)

ax.axhline(1, color='gray', linestyle='--', linewidth=1)
ax.set_ylabel("Fold Change (FC)")
ax.legend()
plt.savefig('AA_carbon_ALL_FC.svg', format='svg') 


plt.tight_layout()
plt.show()
