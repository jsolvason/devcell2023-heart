import matplotlib.pyplot as plt
import numpy as np

indices_to_mutate = [0, 1, 6, 7]
bases = ['A', 'T', 'C', 'G']

def get_mutations(orig):
    mutations = []
    for index_to_mutate in indices_to_mutate:
        for base in bases:
            if orig[index_to_mutate] != base:
                mutation = orig[:index_to_mutate] + base + orig[index_to_mutate+1:]
                mutations.append(mutation)
    return mutations

# Generate dict mapping sequence to intensity
ets2intensity={}
Ets_maxIntensity=0
for line in open("Ets1_8mers.txt").read().rstrip().split("\n"): # reference PBM data
    a=line.split("\t")
    intensity=float(a[3])
    ets2intensity[a[0]]=intensity
    ets2intensity[a[1]]=intensity
    if a[0]=="CCGGAAGT" or a[1]=="CCGGAAGT":
        Ets_maxIntensity=intensity

# Normalize intensities
ets2intensity_normalized = {}
for tfbs, intensity in ets2intensity.items():
    ets2intensity_normalized[tfbs] = float(intensity)/Ets_maxIntensity

def calc_all_fold_changes(orig):
    fold_changes = []
    mutations = get_mutations(orig)
    orig_intensity = ets2intensity_normalized[orig]
    for mutation in mutations:
        fold_change = ets2intensity_normalized[mutation]/orig_intensity
        fold_changes.append(fold_change)
    return fold_changes

# x values on plot should correspond to which base pair was mutated
ji = 1
x_vals = np.array([0, ji, ji*2, 1, 1+ji, 1+ji*2, 6, 6+ji, 6+ji*2, 7, 7+ji, 7+ji*2]) # jitter
x_offsets = [0, 16, 41, 70, 103, 120]

# Determining ETS sites from foxf enhancer in putativeEnhancers_Ci.txt
for line in open("putativeEnhancers_Ci.txt").read().rstrip().split("\n"): # input file name
    if line[0]==">":
        Name=line.split(">")[1]
        continue
    seq=line

ETS_sites_to_plot = []
for i in range(2,len(seq)-5):
    substring4mer = seq[i:i+4]
    if substring4mer == "GGAA" or substring4mer == "GGAT" or substring4mer == "TTCC" or substring4mer == "ATCC":
        ETS_sites_to_plot.append(seq[i-2:i+6])

all_x_vals = np.array([])
all_y_vals = []
for i, ETS_site in enumerate(ETS_sites_to_plot):
    these_x_vals = x_offsets[i] + x_vals
    all_x_vals = np.concatenate((all_x_vals, these_x_vals))
    all_y_vals += calc_all_fold_changes(ETS_site)

plt.figure(figsize=(3.3, 3.6))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.plot([-5, 133], [3, 3], 'b--')
plt.plot([-5, 133], [1, 1], 'c--')
plt.plot(all_x_vals, all_y_vals, 'ko',markersize=2)
plt.ylim([-0.2, 8.5])
plt.xlim([-5, 133])
plt.yticks(ticks=list(range(0, 9)))
plt.ylabel("Affinity fold change")
plt.tight_layout()
plt.savefig("affinity_fold_changes.png",dpi=199)
plt.show()