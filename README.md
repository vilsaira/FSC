# FSC
Functio-Structural Current (FSC) analysis

This script contains all functions needed for the Functio-Structural Current
(FSC) analysis. FSC can be used to combine functional and structural
connectivity matrices[1].

Nothing prevents from using FSC on networks beyond neuroscience. For example, path finding can be done by setting voltage source between start/end and resistances according to the network's edge weights. The shortest path is found by following max current in FSC. As a bonus, this gives all paths.

Author: Viljami Sairanen

[1] Sairanen, Viljami. ‘Combining Function and Structure in a Single Macro-Scale Connectivity Model of the Human Brain’. bioRxiv, 6 March 2024. https://doi.org/10.1101/2024.03.03.583186.

# A Python example of usage with the Human Connectome Project from the ENIGMA TOOLBOX.
You might need to install enigmatoolbox and nilearn first.

```
from enigmatoolbox.datasets import load_sc, load_fc
from nilearn import plotting
from fsc import fsc

# Load cortico-cortical functional connectivity data
fc_ctx, fc_ctx_labels, _, _ = load_fc() 
# Load cortico-cortical structural connectivity data
sc_ctx, sc_ctx_labels, _, _ = load_sc() 
# Calculate cortico-cortical 'functio-structural current' connectivity data
fscc_ctx = fsc(V=fc_ctx, R=sc_ctx).get_I()

# Plot cortico-cortical connectivity matrices
fc_plot = plotting.plot_matrix(fc_ctx, figure=(9, 9), labels=fc_ctx_labels, vmax=0.8, vmin=0, cmap='Reds')

sc_plot = plotting.plot_matrix(sc_ctx, figure=(9, 9), labels=sc_ctx_labels, vmax=10, vmin=0, cmap='Blues')

fsc_plot = plotting.plot_matrix(fscc_ctx, figure=(9, 9), labels=sc_ctx_labels, vmax=fscc_ctx.max(), vmin=fscc_ctx.min(), cmap='Greens')
```
