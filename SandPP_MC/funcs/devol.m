function k = devol(T, x_O2, dp0)
% Copyright 2022, All Rights Reserved
% Code and kinetics data by Stefano Iannello
% For Paper, "Investigation of single particle devolatilization in 
%        fluidized bed reactors by X-ray imaging techniques"
% by S. Iannello, P. U. Fiscolo, M. Materazzi

if x_O2 == 0

    % N2 kin parameters (exp) PP
    Ar = 3.3;                                                               % Reference pre-exp factor [1/s]
    psi = 0.87;                                                             % Fitting parameter
    E = 27.1 * 10^3;                                                        % Activation energy [J/mol]
else

    % AIR kin parameters (exp) PP
    Ar = 7.5;
    psi = 0.8;
    E = 32 * 10^3;
end

R = 8.314;                                                                  % Universal gas constant [J/K mol]
A = Ar * (8 * 10^-3 / dp0)^psi;                                             % Actual pre-exponential factor [1/s]
k = A * exp(-E / (R * T));                                                  % Kinetic rate constant pseudo-first order [1/s]
end