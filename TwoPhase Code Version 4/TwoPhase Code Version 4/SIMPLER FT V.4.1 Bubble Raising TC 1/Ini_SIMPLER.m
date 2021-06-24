% =======================================================================
% Inicializacion SIMPLER
% =======================================================================

% =======================================================================
% Inicializacion --------------------------------------------------------
it=0;

ErrorNitPE=zeros(Nt,1); niteraNitPE=zeros(Nt,1);
ErrorNitPC=zeros(Nt,1); niteraNitPC=zeros(Nt,1);

% Division de Intervalos (Para guardar la data) -------------------------
NInter=1; ErrorNit=zeros(Nt,1); Nit=0:Nt; niteraNit=zeros(Nt,1);
% -----------------------------------------------------------------------


% =======================================================================
% % %%%% Continue from Previous
% load('lastresults.mat'); dt=dt*50; 
% Nt=Nt+2; Nit=0:Nt;
% ErrorNit0=ErrorNit; niteraNit0=niteraNit;
% ErrorNitPE0=ErrorNitPE; niteraNitPE0=niteraNitPE;
% ErrorNitPC0=ErrorNitPC; niteraNitPC0=ErrorNitPC;
% uvels10=uvels1; vvels10=vvels1;
% 
% ErrorNit=zeros(Nt,1); niteraNit=zeros(Nt,1);
% ErrorNitPE=zeros(Nt,1); niteraNitPE=zeros(Nt,1);
% ErrorNitPC=zeros(Nt,1); niteraNitPC=zeros(Nt,1);
% uvels1=zeros(Nt,1); vvels1=zeros(Nt,1);
% =======================================================================