% =======================================================================
% Parametros numericos del algoritmo SIMPLER
% =======================================================================

% Datos numericos para el ciclo SIMPLER en cada IT ----------------------
Error=100; nitera=0; nmax=15; % Nsim=1;
omeu=0.5; omev=0.5; omep=1.5;
%     Tol=1; TolP=Tol/1000; nmaxP=Nx*Ny; % For Criterium a
Tol=0.01; 

% Parameters for the solution of Poisson eqs ----------------------------
nmaxP=5*Nx*Ny; % For Criterium c
% TolP=Tol*min(min(abs(0.5.*(b(2:Nx,2:Ny+1)+b(3:Nx+1,2:Ny+1)).*(dy_un(2:Nx,2:Ny+1)+...
%     dy_su(2:Nx,2:Ny+1)))))*1e-1;
TolPE = Tol * 1e-3;     TolPC = Tol * 1e-2;
% ------------------------------------------------------------------------

% Parameters for the solution of pseudo-velocities ----------------------
fcPV=0; fdPV=0; fSuPV=0; fSvPV=0; % Time Eval Weight,
% -----------------------------------------------------------------------
% Parameters for the solution of the Momentum Eq. -----------------------
fcM=1; fdM=1; fSuM=1; fSvM=1; % Time Eval Weight,
TolM=Tol; 
% ----------------------------------------------------------------------
        % fc=0.5; fd=0.5; fSu=0.5; fSv=0.5; % Time Eval Weight,
        % 0->Full explicit, 1-> Full implicit, 0.5-> ClankNicolson
        % fc -> For Convective terms, fd -> For difusion terms
% -----------------------------------------------------------------------
