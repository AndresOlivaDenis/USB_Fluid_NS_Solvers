% ======================================================================= %
% Datos de entrada
% =======================================================================

% Longuitud y Ubicacion de las fronteras ----------------------------
Lx = 1; X1 = 0; X2 = X1 + Lx;   % [m] Longuitud y Ubicacion de las fronteras
Ly = 1; Y1 = 0; Y2 = Y1 + Ly;   % [m] Longuitud y Ubicacion de las fronteras
B = 1; % [m] Ancho o espesor
% -------------------------------------------------------------------

% Caracteristicas del mallado ---------------------------------------
Ny = 32; Nx = 32;           % Numero de nodos
% -------------------------------------------------------------------

% Propiedades Fisicas -----------------------------------------------
% MIU= 1 * 1 / 400.0;         % [Ns/M2] Viscosidad
% RO=1; roref=RO;             % [Kg/m3] Densidad
% Nt=20000;                    % Numero de avances de paso de tiempo
% dt=80.0/Nt;                 % [s] Paso de tiempo

MIU=1; % [Ns/M2] Viscosidad
RO=1;   roref=RO; % [Kg/m3] Densidad
Nt=20000; % Numero de avances de paso de tiempo
tf=80;
Re=1000;
dt=1e-2/Re; % [s] Paso de tiempo

% -------------------------------------------------------------------

% Termino fuente ----------------------------------------------------
SU=0; SV=0;
% -------------------------------------------------------------------


% Condiciones de borde ----------------------------------------------
uCBnorth=100; uCBsouth=0; uCBwest=0; uCBeast=0; 
vCBnorth=0; vCBsouth=0; vCBwest=0; vCBeast=0;
% Tipo de condicion de borde 0-> Dirichtlet, 1-> Newman 
CuCB=[0,0,0,0]; CvCB=[0,0,0,0]; % 1->West, 2-> North, 3-> East, 4-> South
% -------------------------------------------------------------------

% Condiciones iniciales ---------------------------------------------
U0=0; V0=0;
% -------------------------------------------------------------------

% Plot mallado -----------------------------------------------------
% 1-> Plot of malla, otherwise -> no plot
PlotM = 1;
% ------------------------------------------------------------------

% Plot Avance Iteraciones ------------------------------------------
% 1-> Velocity Field + Error, 2-> Velocity Field, otherwise -> no plot
Plot1=1; 
% ------------------------------------------------------------------

