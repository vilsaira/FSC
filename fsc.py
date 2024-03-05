#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script contains all functions needed for the Functio-Structural Current
(FSC) analysis. FSC can be used to combine functional and structural
connectivity matrices[1].

Author: Viljami Sairanen

[1] Viljami Sairanen, Combining function and structure in a single macro-scale
connectivity model of the human brain, bioarxiv


Example of usage with the Human Connectome Project from the ENIGMA TOOLBOX.

from enigmatoolbox.datasets import load_sc, load_fc
from nilearn import plotting
from fsc import fsc

# Load cortico-cortical functional connectivity data
fc_ctx, fc_ctx_labels, _, _ = load_fc() 
# Load cortico-cortical structural connectivity data
sc_ctx, sc_ctx_labels, _, _ = load_sc() 
# Calculate cortico-cortical 'functio-structural current' connectivity data
fscc_ctx = fsc(V=fc_ctx, R=sc_ctx).get_I()

#  # Plot cortico-cortical connectivity matrices
fc_plot = plotting.plot_matrix(fc_ctx, figure=(9, 9), labels=fc_ctx_labels, vmax=0.8, vmin=0, cmap='Reds')
sc_plot = plotting.plot_matrix(sc_ctx, figure=(9, 9), labels=sc_ctx_labels, vmax=10, vmin=0, cmap='Blues')
fsc_plot = plotting.plot_matrix(fscc_ctx, figure=(9, 9), labels=sc_ctx_labels, vmax=fscc_ctx.max(), vmin=fscc_ctx.min(), cmap='Greens')

"""

import numpy as np
import scipy.linalg as linalg
import sys

class fsc(object):

    def __init__(self, V, R):
        self._V = V
        self._R = R
        self.get_I()

    def _validateSymmetry(self, matrix, matrix_name):
        if (np.abs(matrix - matrix.T).max() > 1e-8):
            sys.exit(f"Error in setting {matrix_name}: {matrix_name} is not symmetric along diagonal.")
        return 0
    
    def _validateDimension(self, A, B, A_name, B_name):
        if (A is not None) & (B is not None):
            if (A.shape != B.shape):
                sys.exit(f"Error in setting {A_name}: dimensions don't match with {B_name}.")

    def get_I(self):
        voltage_sources = np.argwhere(self._V != 0)
        M = len(voltage_sources)
        N = self._R.shape[0] # number of nodes or vertices

        G = np.zeros((N,N)) # Conductance matrix G: 
        # G_ii is the sum of all conductances connected to node i
        # G_ij is the negative of the sum of all conductances connected from i to j
        for i in range(N):
            for j in range(N):
                if self._R[i,j] > 0:
                    G[i,i] += 1 / self._R[i,j]
                    G[i,j] -= 1 / self._R[i,j]

        B = np.zeros((N,M))
        # If the positive/negative terminal of the i-th voltage source is connected to node j, then the element (i,j) in the B matrix is 1 / -1.
        for i, inds in enumerate(voltage_sources):
            volt = self._V[inds[0], inds[1]]
            if volt != 0:
                B[inds[0], i] = np.sign(volt)
                B[inds[1], i] = -np.sign(volt)

        C = B.T # If and only if there exist no dependent sources! In our analogy such are not present and C is the transpose of B.

        D = np.zeros((M,M)) # If and only if there exist no dependent sources!

        # A = [[G, B], [C, D]]
        A = np.concatenate((G,B), axis=1)
        tmp = np.concatenate((C,D), axis=1)
        A = np.concatenate((A,tmp), axis=0)

        z = np.zeros((N+M,)) # There are no current sources thus sum of entering/leaving currents in each node must be zero
        for i in range(M):
            volt = self._V[voltage_sources[i][0], voltage_sources[i][1]]
            z[i+N] = volt

        # Set the first node as the ground and remove corresponding elements from A and z
        A1 = A[1:,1:]
        z1 = z[1:]

        x = linalg.lstsq(A1, z1)[0]

        # Calculate voltage differences. The first node has voltage of 0 as it is ground.
        node_voltages = np.zeros((1,))
        node_voltages = np.concatenate((node_voltages, x[0:N-1]), axis=0)

        U = np.zeros((N,N)) # Voltage difference matrix
        for i in range(N-1):
            for j in range(i+1, N):
                if self._R[i,j] > 0:
                    U[i,j] = np.abs(node_voltages[j] - node_voltages[i])

        # Solve the current matrix
        I = U / self._R
        np.fill_diagonal(I, 0)
        I[np.isnan(I)] = 0.0

        return I + I.T
    
    @property
    def _V(self):
        return self.__V
    
    @_V.setter
    def _V(self, value):
        self._validateSymmetry(value, "V")
        # self._validateDimension(value, self._R, "V", "R")
        self.__V = np.triu(value)

    @property
    def _R(self):
        return self.__R
    
    @_R.setter
    def _R(self, value):
        self._validateSymmetry(value, "R")
        self._validateDimension(value, self._V, "R", "V")
        self.__R = value

    if __name__ == "__main__":
        print("No main functions implemented yet")
