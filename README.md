# RaFiNe
Code to model the mechanics of RAndom FIlament NEtworks.
This repository is tightly related to a paper published in International Journal of Solids and Structures:

Róbert K. Németh, András Á. Sipos, Layth S. Al-Rukaibawi, Lili E. Hlavicka-Laczák, Flórián Kovács, György Károlyi:</br>
Local shear stiffness switches between bending- and stretch-dominated regimes of a random net of elastic filaments.</br>
International Journal of Solids and Structures 321 (2025) 113576.</br>
https://doi.org/10.1016/j.ijsolstr.2025.113576

Please, find all relevant information in the paper.

This code was written by:</br>
György Károlyi</br>
Róbert K. Németh</br>
Flórián Kovács</br>

Files:</br>
rafine.c - Source code, computes the equlibrium configuration of a random network of elastic filaments.</br>
results.dat - Data processed from the output of many runs of rafine.c with varied parameters, not including anizotropy.</br>
results-anizo.dat - Data processed from the output of many runs of rafine.c with varied parameters, including the anizotropy parameter.

These are movie files that correspond to Fig. 4 of the paper:</br>
stiffness-total.gif - Fig. 4a -  Total shear stiffness of the network as a function of local shear (σ) and bend (ρ) moduli.</br>
stiffness-stretch.gif - Fig. 4b -  Stretch component of the total shear stiffness of the network as a function of local shear (σ) and bend (ρ) moduli.</br>
stiffness-bend.gif - Fig. 4c -  Bending component of the total shear stiffness of the network as a function of local shear (σ) and bend (ρ) moduli.</br>
stiffness-shear.gif - Fig. 4d -  Shear component of the total shear stiffness of the network as a function of local shear (σ) and bend (ρ) moduli.
