#スペクトル取得のプログラム
import numpy as np
import matplotlib.pyplot as plt

def load_spectrum(filename):
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if line.strip() == "[Data]":
                skiprows = i + 2
                break
        else:
            raise ValueError("[Data] section not found")

    ch, count = np.loadtxt(
        filename,
        delimiter=",",
        skiprows=skiprows,
        unpack=True
    )

    return ch, count

ch, count =  load_spectrum("compton90.csv")

plt.plot(ch, count)
plt.xlabel("Gamma-ray Energy [ch]")
plt.ylabel("Counts")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.savefig("../tex/compton90.pdf", dpi=300, bbox_inches="tight")
plt.show()