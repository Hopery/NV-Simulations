This repository contains simulations of some experiments with NV-Centers. It's divided in two folders, **julia** and **matlab**, each one has code in its programme language.

The folder Julia contains two different simulations:
1. Three level spin system simulations for the following experiments in a simple way i.e. not taking into account carbon13 nucleus:
    * Observe Spins under magnetic field
    * ODMR experiment
    * Rabi experiment 
2. Taking into account the interactions between the nv-center and the carbon13:
    * Ramsey experiment
    * SpinEcho experiment with XYXYYX sequence
    * SpinEcho experiment with XYXYYX sequence and fixed τ

The matlab folder contains simulations taking into account the interactions between the nv-center and the carbon13, these simulations are:
* Rabi experiment 
* Ramsey experiment
* SpinEcho experiment with XYXYYX sequence
* SpinEcho experiment with XYXYYX sequence and fixed τ
* SpinEcho experiment with XYXYYX sequence and fixed τ, obtaining it by Fourier components
*  CNOT gate
* SWAP gate
* Sequence of SWAP - RF - SWAP
* ODMR of last sequence