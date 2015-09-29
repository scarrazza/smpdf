# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 10:41:50 2015

@author: zah
"""
import matplotlib.pyplot as plt

from smpdflib.initialization import init_style
from matplotlib.ticker import MaxNLocator


init_style()

thresholds = [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]
results_comb = [49, 32, 22, 20, 18, 19]
results_nnpdf = [61, 27, 21, 15, 15, 13]
results_nnpdf2 = [62, 28, 20, 14, 15, 16]

results_combined_higgs = [26, 23, 18, 10, 10, 9]
results_nnpdf_higgs = [29, 20, 13, 9, 7, 8]

results_combined_w = [35, 27, 22, 15, 14, 14]

plt.figure(figsize=(8,6))
#plt.plot(thresholds, results_comb, marker="s", label = "MC900")
#plt.plot(thresholds, results_nnpdf, marker="s", label= "NNPDF 3.0")
#plt.plot(thresholds, results_nnpdf2, marker="s", label= "NNPDF 3.0")

plt.plot(thresholds, results_comb, marker="s", label = "Ladder")
plt.plot(thresholds, results_combined_higgs, marker="s", label = "Hoggs")
plt.plot(thresholds, results_combined_w, marker="s", label = "W")



plt.xlabel("$t$")
plt.ylabel("$N_{eig}$")
plt.title("$N_{eig}$ at fixed $T_R=5\%$ for MC900 NLO", y=1.01)
plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
plt.legend()
plt.tight_layout()
plt.savefig("thresholds.pdf")
