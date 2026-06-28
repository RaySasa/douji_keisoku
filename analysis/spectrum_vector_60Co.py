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

ch, count =  load_spectrum("data/60Co_spectrum_secondtime.csv")

plt.annotate("光電吸収ピーク", xy=(305, 2700), xytext=(305, 3000), ha ="center", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:orange",linewidth=1.5), fontsize=9, color ="tab:orange")
plt.annotate("Compton端", xy=(250, 1600), xytext=(250, 2100), ha ="center", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:orange",linewidth=1.5), fontsize=9, color ="tab:orange")
plt.annotate("後方散乱ピーク", xy=(65, 3150), xytext=(120, 3150), ha ="left", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:orange",linewidth=1.5), fontsize=9, color ="tab:orange")

plt.annotate("光電吸収ピーク", xy=(345, 1900), xytext=(400, 1900), ha ="left", va ="center", arrowprops=dict(arrowstyle="->", color ="tab:green",linewidth=1.5), fontsize=9, color ="tab:green")

plt.plot(ch, count)
plt.xlabel("Gamma-ray Energy[ch]")
plt.ylabel("Counts")
plt.grid(linestyle = "--", linewidth = 0.5)
plt.savefig("tex/60Co_spectrum_vector.pdf", dpi=300, bbox_inches="tight")
plt.show()