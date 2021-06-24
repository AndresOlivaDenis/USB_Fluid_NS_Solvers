% ======================================================================= %
% Datos de entrada
% =======================================================================

% Longuitud y Ubicacion de las fronteras ----------------------------
Lx = 1; X1 = 0; X2 = X1 + Lx;   % [m] Longuitud y Ubicacion de las fronteras
Ly = 2; Y1 = 0; Y2 = Y1 + Ly;   % [m] Longuitud y Ubicacion de las fronteras
B = 1; % [m] Ancho o espesor
% -------------------------------------------------------------------

% Caracteristicas del mallado ---------------------------------------
h = 40;
Ny = 2*h; Nx = h;           % Numero de nodos
% -------------------------------------------------------------------

% Propiedades Fisicas -----------------------------------------------
sigma = 24.5; % [??] Tension Superficial 
MIU1 = 10;      MIU2 = 1; % [Ns/M2] Viscosidad 
RO1 = 1000;     RO2 = 100;      roref=RO1; % [Kg/m3] Densidad 
dt = (1/h)/16*(1/4); % [s] Paso de tiempo
Tfinal=3; Nt=Tfinal/dt; % Numero de avances de paso de tiempo
% -------------------------------------------------------------------

% Gota --------------------------------------------------------------
rdrop=0.25;xdrop=0.5;ydrop=0.5; % Initial drop size and location
% -------------------------------------------------------------------

% Termino fuente ----------------------------------------------------
gx=0; gy=-0.98;
% -------------------------------------------------------------------


% Condiciones de borde ----------------------------------------------
%-> Inlet (y=0), Outlet (y=Ly), Walls (:,:)
uCBnorth=0; uCBsouth=0; uCBwest=0; uCBeast=0; 
vCBnorth=0; vCBsouth=0; vCBwest=0; vCBeast=0;
% Tipo de condicion de borde 0-> Dirichtlet, 1-> Newman 
CuCB=[0,0,0,0]; CvCB=[1,0,1,0]; % 1->West, 2-> North, 3-> East, 4-> South
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

