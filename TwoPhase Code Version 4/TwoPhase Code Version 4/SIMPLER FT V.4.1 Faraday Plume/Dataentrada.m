% ======================================================================= %
% Datos de entrada
% =======================================================================

% Longuitud y Ubicacion de las fronteras ----------------------------
NumP = 1; LxN = 1*NumP;             % Periodicity "copy"
Lx = LxN; X1 = 0; X2 = X1 + Lx;     % [m] Longuitud y Ubicacion de las F.
Ly = 1; Y1 = 0; Y2 = Y1 + Ly;       % [m] Longuitud y Ubicacion de las F.
B = 1;                              % [m] Ancho o espesor
% -------------------------------------------------------------------

% Caracteristicas del mallado ---------------------------------------
Ny = 60; Nx = 128*NumP;             % Numero de nodos
% -------------------------------------------------------------------

% Termino fuente ----------------------------------------------------
gx = 0.0;       gy = 0.0;
% -------------------------------------------------------------------

% Propiedades Fisicas -----------------------------------------------
sigma = 7.2e-2;                                 % [??] Tension Superficial 
MIU1 = 0.2;     MIU2 = 0.2;                     % [Ns/M2] Viscosidad 
RO1 = 1000;     RO2 = 818.18;      roref=RO1;   % [Kg/m3] Densidad 

lc=(sigma/abs(gy*(RO1-RO2)))^0.5;               % Capillary Length; Ly=5*lc;

% Faraday data .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
ac = 9e-2;                              % Critical amplitude
k = (2*pi);     landa = 2*pi/k;         % wavelength (landa), wavenumber (k); 
amp = 0.01;     ysurf = 0.5;            % Initial Pertubation configuration 
ome = 2.6728e-1;        Tv = (2*pi)/ome;
% .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

dt = Tv/70*0.1;                                 % [s] Paso de tiempo
Tfinal = 5*Tv; Nt = floor(Tfinal/dt)+1;         % Time Step Things
% -------------------------------------------------------------------

% Condiciones de borde ----------------------------------------------
%-> Inlet (y=0), Outlet (y=Ly), Walls (:,:)
uCBnorth=0; uCBsouth=0; uCBwest=0; uCBeast=0; 
vCBnorth=0; vCBsouth=0; vCBwest=0; vCBeast=0;
% Tipo de condicion de borde 0-> Dirichtlet, 1-> Newman 
CuCB=[0,0,0,0]; CvCB=[0,0,0,0]; % 1->West, 2-> North, 3-> East, 4-> South
% Notaaa! Periodicidad!!
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

