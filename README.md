# FreeElectronSpinPolarizerCode
The included matlab code are used to simulate and produce the figures in the main text of, DOI: https://doi.org/10.1103/3c1m-d3hh. This code performs implementations of the equations provided in the text based on the Pauli equation of a free electron interacting with two stages of near-fields. This work explores the spin and momentum characteristics of the resulting free electron wavepacket.

SpinExpectation.m plots the spatial spin probability density from the simplified two stage interaction. User inputs include interaction strength, phase, and lengths. The total spin expectation value is also printed as a ratio corresponding to 1 being fulling spin polarized positive and -1 being fully spin polarized negative.

PauliSplitStep.m performs a split step analysis of the section stage interaction of the main text. This analysis includes both ponderomotive effects in the second stage, as well as dispersion effects resulting from the interaction of the second stage and the momentum spectrum of the first stage. The input wavefunction is assumed to be that of the first interaction stage provided in the main text.

HeatMapG.m plots the heat maps provided in Figure 3. Similarly, the user inputs include interaction strength and phase of the two stages, pixel number, and what parameters are to be plotted.
