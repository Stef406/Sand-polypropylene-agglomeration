# Sand-Polypropylene Agglomeration in Fluidized Beds

Copyright 2024, All Rights Reserved

Code by Stefano Iannello
For Paper, "The behaviour of plastic particles during pyrolysis in 
        bubbling fluidized bed reactors: Incipient agglomeration and 
        axial segregation"
Powder Technology, 441, 119846, 2024,
https://doi.org/10.1016/j.powtec.2024.119846,
by S. Iannello, A. Sebastiani, M. Errigo, M. Materazzi.

This project provides a Monte Carlo model in MATLAB to simulate the formation of a sand-polypropylene agglomerate within a fluidized bed reactor during pyrolysis.

The main algorithm is PP_MC.m.

Useful functions are in the directory "funcs":
  - bed.m is to compute the properties of the fluidized bed reactor used in the study.
  - devol.m is to compute the devolatilization rate constant of polypropylene in pyrolysis conditions (for more information about this topic, check https://doi.org/10.1016/j.cej.2021.133807).

Agglomeration.avi shows an animation of the agglomeration between sand particles and molten PP.
