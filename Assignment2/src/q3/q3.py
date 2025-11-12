import numpy as np
import pandas as pd

genes = ["G1", "G2", "G3", "G4"]
MI = pd.DataFrame([
    [1.24, 1.01, 0.69, 0.69],
    [1.01, 1.56, 1.01, 0.78],
    [0.69, 1.01, 1.01, 0.46],
    [0.69, 0.78, 0.46, 1.01]
], index=genes, columns=genes)

# calc Zi
Zi = MI.apply(lambda x: (x - x.mean()) / x.std(ddof=0), axis=1)
Zi = Zi.clip(lower=0)

# calc Zj
Zj = MI.apply(lambda x: (x - x.mean()) / x.std(ddof=0), axis=0)
Zj = Zj.clip(lower=0)

# calc CLR mat
CLR = np.sqrt(Zi ** 2 + Zj ** 2)
CLR = pd.DataFrame(CLR, index=genes, columns=genes)

pairs = []
for i in range(len(genes)):
    for j in range(i + 1, len(genes)):
        pairs.append((genes[i], genes[j], CLR.iloc[i, j]))

top3 = sorted(pairs, key=lambda x: x[2], reverse=True)[:3]

print("Zi mat：\n", Zi.round(3))
print("\nZj mat：\n", Zj.round(3))
print("\nCLR mat：\n", CLR.round(3))
print("\nTop 3 GGI：")
for g1, g2, val in top3:
    print(f"{g1}-{g2}: {val:.3f}")
