% =======================================================================
% Solution of 2D Unsteady momentum eq.
% SIMPLER Method
% Time Disc -> General ClankNicolson    Space -> Center Diference
% =======================================================================
clc; clear; close all

% ======================================================================= 
% Datos de entrada
run Dataentrada.m
run ParametersSIMPLER
% =======================================================================

% =======================================================================
% A Staggered Grid
% Control Volume for u & v -> Forward U,V Volume Control
[xP,yP,dx_PE,dx_Pe,dy_PN,dy_Pn,dx_WP,dx_wP,dy_SP,dy_sP,...
    xu,yu,dx_uE,dx_ue,dy_uN,dy_un,dx_Wu,dx_wu,dy_Su,dy_su,...
    xv,yv,dx_vE,dx_ve,dy_vN,dy_vn,dx_Wv,dx_wv,dy_Sv,dy_sv,...
    FMx,FMy,xi,yj,psi,psj]...
    =StaggeredGrid_3_asFunction(X1,X2,Y1,Y2,Nx,Ny);

run PlotMalla.m;    
% =======================================================================

% ======================================================================= 
% Valores Iniciales
run InicializacionPro.m
% =======================================================================

% =======================================================================
% Condiciones de borde para el campo de velocidades inicial
[u0,v0] = Velocities_BC(uCBeast, uCBwest, vCBsouth, vCBnorth,...
    uCBsouth, uCBnorth, vCBwest, vCBeast, CuCB, CvCB,...
    dx_uE, dx_Wu, dy_vN, dy_Sv, dy_uN, dx_Wv, dy_Su, dx_vE, ...
    Nx,Ny,u0,v0);
u=u0; v=v0;         us=u0; vs=v0;         ut=u0; vt=v0;
% =======================================================================
  
% =======================================================================
% Initialization of Benchmark quantities
XcBubble=zeros(Nt+1,2); CirBubble=zeros(Nt+1,1); VelBubble=zeros(Nt+1,2);
XcBubble(1,1)=xdrop; XcBubble(1,2)=ydrop;  CirBubble(1,1)=1;
VelBubble2=zeros(Nt+1,2);
% =======================================================================

% =======================================================================
% Start of The SIMPLER algorithm

% Inicializacion de las variables ---------------------------------------
run Ini_SIMPLER.m
% -----------------------------------------------------------------------


for iinte=1:NInter
for it=((iinte-1)*Nt/NInter+1):((iinte)*Nt/NInter)
% for it=1:Nt+1

    % Information for Print --------------------------------------------
    disp('------------------------------------------------------')
    info=strcat('TimeStep # = ',num2str(it),...
            '_ t = ',num2str(it*dt));               disp(info)
    % ------------------------------------------------------------------
    
    % ==================================================================
    % Last Time Step Values --------------------------------------------
    u0=u; v0=v; ro0=ro; miu0=miu; MarkF0=MarkF;
    % Front Variables -----
    NFront0=NFront; xFront0=xFront; yFront0=yFront;
    xmapFront0=xmapFront; ymapFront0=ymapFront;
    % uFront0=uFront; vFront0=vFront; 
    rodxFront0=rodxFront; rodyFront0=rodyFront;
    umapFront0=umapFront; vmapFront0=vmapFront;
    sftx0=sftx; sfty0=sfty;
    % ------------------------------------------------------------------
    
    % Local Time step iteration ----------------------------------------
    Error=100; nitera=1; 
    
    while (Error>=Tol*2)&&(nitera<=nmax)
        % Information for print ----------------------------------------
        disp('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-')
        info=strcat('SIMPLER Method #Iteracion = ',num2str(nitera),...
            '_ Error = ',num2str(Error));
        disp(info)
        % --------------------------------------------------------------
        
        nitera=nitera+1; ulast=u; vlast=v; plast=p;
        
        % Calculo de las Pseudo-Velocidades ----------------------------
        [us,vs,aPus,aPvs] = Mtum2PseudoMiu2_STFg_asFunction(Nx,Ny,...
            dx_uE,dx_ue,dy_uN,dy_un,dx_Wu,dx_wu,dy_Su,dy_su,...
            dx_vE,dx_ve,dy_vN,dy_vn,dx_Wv,dx_wv,dy_Sv,dy_sv,...
            ro,miu,b,Su,Sv,roref,...
            u0,v0,ro0,miu0,dt,...
            uCBnorth,uCBsouth,uCBwest,uCBeast,vCBnorth,vCBsouth,vCBwest,vCBeast,...
            CuCB,CvCB,fcPV,fdPV,fSuPV,fSvPV,...
            sftx,sfty,sftx0,sfty0,TolM,...
            us,vs);
        % --------------------------------------------------------------
        
        % The Pressure Equation ----------------------------------------
        [p,niteraPE,ErrorPE,aPpc] = PressureRo2_3_asFunction(Nx,Ny,...
            dx_Pe,dy_Pn,dx_wP,dy_sP,dy_un,dy_su,dx_ve,dx_wv,...
            ro,us,vs,aPus,aPvs,p,dt,ro0,nmaxP,TolPE);
        % --------------------------------------------------------------

        % Solution of Unsteady Momentum Eqs 2D -------------------------
        % ut = u; vt = v;
        [ut,vt,aPu,aPv] = Momentum2Miu2_STFg_asFunction(Nx,Ny,...
            dx_uE,dx_ue,dy_uN,dy_un,dx_Wu,dx_wu,dy_Su,dy_su,...
            dx_vE,dx_ve,dy_vN,dy_vn,dx_Wv,dx_wv,dy_Sv,dy_sv,...
            ro,miu,b,Su,Sv,roref,...
            p,u0,v0,ro0,miu0,dt,...
            uCBnorth,uCBsouth,uCBwest,uCBeast,vCBnorth,vCBsouth,vCBwest,vCBeast,...
            CuCB,CvCB,fcM,fdM,fSuM,fSvM,...
            sftx,sfty,sftx0,sfty0,TolM,...
            ut,vt);
        % --------------------------------------------------------------
        
        % The Pressure Correction Equation ------------------------------
        [pc,niteraPC,ErrorPC] = PressureCorrecRo2_3_asFunction(Nx,Ny,...
            dx_Pe,dy_Pn,dx_wP,dy_sP,dy_un,dy_su,dx_ve,dx_wv,...
            ro,ut,vt,aPu,aPv,pc,dt,ro0,nmaxP,TolPC);
        % --------------------------------------------------------------
        
        % Correct Velocities -------------------------------------------     
        u(2:Nx,2:Ny+1)=ut(2:Nx,2:Ny+1)+...
            0.5.*(b(2:Nx,2:Ny+1)+b(3:Nx+1,2:Ny+1)).*(dy_un(2:Nx,2:Ny+1)+...
            dy_su(2:Nx,2:Ny+1))...
            ./(aPu(2:Nx,2:Ny+1)).*(pc(2:Nx,2:Ny+1)-pc(3:Nx+1,2:Ny+1));
        v(2:Nx+1,2:Ny)=vt(2:Nx+1,2:Ny)+...
            0.5.*(b(2:Nx+1,2:Ny)+b(2:Nx+1,3:Ny+1)).*(dx_ve(2:Nx+1,2:Ny)+...
            dx_wv(2:Nx+1,2:Ny))...
            ./aPv(2:Nx+1,2:Ny).*(pc(2:Nx+1,2:Ny)-pc(2:Nx+1,3:Ny+1));
        % --------------------------------------------------------------
        
        % Imposicion de las condiciones de borde -----------------------
        [u,v] = Velocities_BC(uCBeast, uCBwest, vCBsouth, vCBnorth,...
            uCBsouth, uCBnorth, vCBwest, vCBeast, CuCB, CvCB,...
            dx_uE, dx_Wu, dy_vN, dy_Sv, dy_uN, dx_Wv, dy_Su, dx_vE, ...
            Nx,Ny,u,v);
        % --------------------------------------------------------------
          
        % ------ Error L2 ----------------------------------------------
        Erroru=(sum((u-ut).^2)/sum(ut.^2))^0.5;
        Errorv=(sum((v-vt).^2)/sum(vt.^2))^0.5;
        Error=max(Erroru,Errorv)*100;
        % --------------------------------------------------------------  
    end
    
    % ==================================================================
    % Front Tracking Method --------------------------------------------
    disp('--- Front Tracking ---')
    % Advection of The Front -------------------------------------------
    [xFront,yFront,uFront,vFront,...
        xmapFront,ymapFront,umapFront,vmapFront]...
        =AdvectionoftheFront03_A(NFront0,xFront0,yFront0,...
        X1,X2,Y1,Y2,Nx,Ny,u,v,u0,v0,xu,yu,xv,yv,dt,...
        xmapFront0,ymapFront0,FMx,FMy);
    % ------------------------------------------------------------------
    
    % ReStructure of The Front -----------------------------------------
%     [NFront,xFront,yFront,xmapFront,ymapFront]...
%         =ReStructureoftheFront02_B(NFront0,xFront,yFront,...
%         FrontDistMax,FrontDistMin,...
%         xmapFront,ymapFront,FMx,FMy,X1,X2,Y1,Y2);
%     info=strcat('# Front Points = ',num2str(NFront)); disp(info)
    
    [xFront,yFront,xmapFront,ymapFront]...
        =ReStructureoftheFront03_B(NFront,xFront,yFront,...
        xmapFront,ymapFront,FMx,FMy,X1,X2,Y1,Y2);
    % ------------------------------------------------------------------

    % Construction of The Marker Function  -----------------------------
    [MarkF,MarkFdx,MarkFdy]...
        =ConstructMarkerFunction03_A(MarkF0,MF01,MF02,Nx,Ny,...
        xu,yu,xv,yv,dx_PE,dy_PN,dx_WP,dy_SP,...
        NFront,xFront,yFront,xmapFront,ymapFront);
    
    for i=2:Nx+1
        for j=2:Ny+1
            MarkF(i,j)=max(MarkF(i,j),0); MarkF(i,j)=min(MarkF(i,j),1);
        end
    end
    % ------------------------------------------------------------------
%     
    % Density & Viscosities --------------------------------------------
    ro=RO2.*MarkF+(1-MarkF).*RO1; miu=MIU2.*MarkF+(1-MarkF).*MIU1;
    % ------------------------------------------------------------------
    
    % -----------------------------------------------------------------
    % Surface Tension Force ------------------------------------------- 
    %                                           using Front markers --
    [sftx,sfty] = SurfaceTensionForce04_D(Nx,Ny,sigma,...
        xu,yu,xv,yv,NFront,xFront,yFront,xmapFront,ymapFront);
    % -----------------------------------------------------------------
    % =================================================================
    disp('------------------------------------------------------')
    
    % -----------------------------------------------------------------
%     % Mass Check ----------------------------------------------------
%     Mass(it)=0;
%     for i=2:Nx+1, for j=2:Ny+1 %#ok<ALIGN>
%             Mass(it)=Mass(it)+...
%                 MarkF(i,j)*(dy_un(i,j)+dy_su(i,j))*(dx_ue(i,j)+dx_wu(i,j));
%         end,end
    % -----------------------------------------------------------------
    
    % -----------------------------------------------------------------
    % BenchMark Quanties
    run BenchQuanties
    % -----------------------------------------------------------------
        
    % -----------------------------------------------------------------
    % plot results!
    % Plot1=1; % 1-> Velocity Field + Error, 2-> Velocity Field, otherwise -> no plot
    run PlotResults
    % -----------------------------------------------------------------
    
end
% Data Save --------------------------------------------------------------
infoSave=strcat('SIMPLER_V4=',num2str(it),'_t',num2str(it*dt),'.mat');
save(infoSave);
% ------------------------------------------------------------------------
end
    
% =======================================================================


% ------------------------------------------------------------------------
% % Video
% writerObj = VideoWriter('newfile2.avi'); writerObj.FrameRate = 10; 
% open(writerObj); writeVideo(writerObj,F); close(writerObj);
% ------------------------------------------------------------------------


% Post plots -------------------------------------------------------------
run PostPlots
% ------------------------------------------------------------------------
