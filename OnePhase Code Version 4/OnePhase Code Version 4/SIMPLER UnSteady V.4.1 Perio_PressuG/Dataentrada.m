% ======================================================================= %
% Datos de entrada
% =======================================================================

% Longuitud y Ubicacion de las fronteras ----------------------------
Lx = 1; X1 = 0; X2 = X1 + Lx;   % [m] Longuitud y Ubicacion de las fronteras
Ly = 0.2; Y1 = 0; Y2 = Y1 + Ly;   % [m] Longuitud y Ubicacion de las fronteras
B = 1; % [m] Ancho o espesor
% -------------------------------------------------------------------

% Caracteristicas del mallado ---------------------------------------
Ny = 10; Nx = 32;           % Numero de nodos
% -------------------------------------------------------------------

% Propiedades Fisicas -----------------------------------------------
MIU= 1 * 1 / 400.0;         % [Ns/M2] Viscosidad
RO=1; roref=RO;             % [Kg/m3] Densidad
Nt=20000;                    % Numero de avances de paso de tiempo
dt=80.0/Nt;                 % [s] Paso de tiempo
% -------------------------------------------------------------------

% Termino fuente ----------------------------------------------------
SU=0; SV=0;
% -------------------------------------------------------------------


% Condiciones de borde ----------------------------------------------
%-> Inlet (y=0), Outlet (y=Ly), Walls (:,:)
uCBnorth=0; uCBsouth=0; uCBwest=0; uCBeast=0; 
vCBnorth=0; vCBsouth=0; vCBwest=0; vCBeast=0;
% Tipo de condicion de borde 0-> Dirichtlet, 1-> Newman 
CuCB=[0,0,0,0]; CvCB=[0,0,0,0]; % 1->West, 2-> North, 3-> East, 4-> South
% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
% Periodicidad en y (North & South)

% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
% -------------------------------------------------------------------

% Gradiente de Presion en y -----------------------------------------
umax = 0.003;
dpdy = umax*2*MIU/(0.5*Lx)^2; % Gradiente de presion longitudinal
% -------------------------------------------------------------------


% Condiciones iniciales ---------------------------------------------
U0=0; V0=0;
% -------------------------------------------------------------------

% Delete this! its only for test -------------------------
Nt = 300;                          % Numero de Avances de paso de tiempo
dt = 0.2e-1;                         % Paso de tiempo
dt = (0.5)/umax*(max(max(Ly/Ny)));
% dpdy = 0;
% Sv = 10;
% --------------------------------------------------------

% Plot mallado -----------------------------------------------------
% 1-> Plot of malla, otherwise -> no plot
PlotM = 1;
% ------------------------------------------------------------------

% Plot Avance Iteraciones ------------------------------------------
% 1-> Velocity Field + Error, 2-> Velocity Field, otherwise -> no plot
Plot1=1; 
% ------------------------------------------------------------------

