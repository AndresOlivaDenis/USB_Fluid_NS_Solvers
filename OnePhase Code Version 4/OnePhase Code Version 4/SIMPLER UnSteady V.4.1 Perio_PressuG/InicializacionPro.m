% ======================================================================= %
% Inicializacion de propiedades y variables
% =======================================================================

% Imposicion de las propiedades Fisicas ---------------------------------
% Propiedades y terminos fuentes constantes:
ro=RO*ones(Nx+2,Ny+2); miu=MIU*ones(Nx+2,Ny+2); b=B*ones(Nx+2,Ny+2);
Su=SU*ones(Nx+1,Ny+2); Sv=SV*ones(Nx+2,Ny+1); ro0=RO*ones(Nx+2,Ny+2);
% -----------------------------------------------------------------------

% Condicion Inicial -----------------------------------------------------
% Velocidades y Presion
u0=U0*ones(Nx+1,Ny+2); v0=V0*ones(Nx+2,Ny+1);
u=u0; v=v0; p=zeros(Nx+2,Ny+2); pc=zeros(Nx+2,Ny+2);
us=u0; vs=v0; ut=u0; vt=v0;
% -----------------------------------------------------------------------
