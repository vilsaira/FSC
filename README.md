# FSC
Functio-Structural Current (FSC) analysis

This script contains all functions needed for the Functio-Structural Current
(FSC) analysis. FSC can be used to combine functional and structural
connectivity matrices[1].

Author: Viljami Sairanen

[1] Viljami Sairanen, Combining function and structure in a single macro-scale
connectivity model of the human brain, bioarxiv, 2024

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
