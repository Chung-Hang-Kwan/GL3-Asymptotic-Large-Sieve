# GL3-Asymptotic-Large-Sieve

Hello World!

This repository contains the Mathematica notebooks used for the paper "Critical Zeros and Unconditional Mean Value Theorems for twisted $\hbox{PGL}(2)$ and $\hbox{PGL}(3)$ $\mathrm{L}$-functions", written by Brian Conrey, Yongxiao Lin, Caroline Turnage-Butterbaugh, and me.

We use the computer package gln.m written by Kevin A. Broughan (as part of Goldfeld's GL(n) book). Please download the gln.m file along with all of our notebooks, and place them in the same folder. To load the gln.m file to a session, use

SetDirectory[NotebookDirectory[]]; 
           << gln.m

## Notebooks

### ShiftNaRS.nb

This notebook computes the "naive" Rankin--Selberg with shifts in the "Preliminary" section  (Section 3) of the paper.

### MainEuExp.nb

This notebook computes the asymptotic of an Euler factor in the diagonal term; see Section 5 of the paper.

### U2pnmid.nb

This notebook computes the Euler factors of the $\mathcal{U}^{2}(h,k)$-sum of the paper (Section 10).

### UrFinExact.nb

This computes the finite Euler product $\mathcal{E}$ over $\mathfrak{k}'$ in the sum $\mathcal{U}^{(r)}(h,k)$; see Section 12.4 of the paper.

### UrInf.nb

This computes and bounds the infinite Euler product of $\mathfrak{A}_{p}$'s in the sum $\mathcal{U}^{(r)}(h,k)$; see Section 12.5 of the paper.


### GL3Ratio.nb

This carries out the "ratios"/mollification calculations in Section 13 of the paper.

Enjoy!




