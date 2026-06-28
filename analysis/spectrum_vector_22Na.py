#スペクトル取得のプログラム
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Hiragino Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

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

ch, count =  load_spectrum("data/22Na_spectrum.csv")

plt.annotate("光電吸収ピーク", xy=(145, 2900), xytext=(200, 2900), ha ="left", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:orange",linewidth=1.5), fontsize=9, color ="tab:orange")
plt.annotate("Compton端", xy=(95, 700), xytext=(170, 700), ha ="left", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:orange",linewidth=1.5), fontsize=9, color ="tab:orange")
plt.annotate("後方散乱ピーク", xy=(50, 1250), xytext=(50, 1600), ha ="center", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:orange",linewidth=1.5), fontsize=9, color ="tab:orange")

plt.annotate("光電吸収ピーク", xy=(337, 320), xytext=(370, 320), ha ="left", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:green",linewidth=1.5), fontsize=9, color ="tab:green")
plt.annotate("Compton端", xy=(280, 120), xytext=(280, 1100), ha ="center", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:green",linewidth=1.5), fontsize=9, color ="tab:green")

plt.plot(ch, count)
plt.xlabel("Gamma-ray Energy [ch]")
plt.ylabel("Counts")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.savefig("tex/22Na_spectrum_vector.pdf", dpi=300, bbox_inches="tight")
plt.show()