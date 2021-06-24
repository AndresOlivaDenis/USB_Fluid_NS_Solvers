% =======================================================================
% Pressure Correction Eq.
% =======================================================================

function [p,niteraPE,ErrorPE] = PressureRo2_3_asFunction(Nx,Ny,...
    dx_Pe,dy_Pn,dx_wP,dy_sP,dy_un,dy_su,dx_ve,dx_wv,...
    ro,us,vs,aPus,aPvs,p,dt,ro0,nmaxP,TolP)

% The Pressure Correction Equation
% =======================================================================
ErrorPE=1; niteraPE=0;  omepc=1.2;

% =======================================================================
% Coeficientes aW,aE,aN,aS,..., y Termino Independiente
aEpc=zeros(Nx+2,Ny+2); aWpc=zeros(Nx+2,Ny+2); aNpc=zeros(Nx+2,Ny+2);
aSpc=zeros(Nx+2,Ny+2); aPpc=zeros(Nx+2,Ny+2); bpc=zeros(Nx+2,Ny+2);
for i=2:Nx+1
    for j=2:Ny+1
        aEpc(i,j)=(dy_Pn(i,j)+dy_sP(i,j))*(dy_un(i,j)+dy_su(i,j))/aPus(i,j);
        aWpc(i,j)=(dy_Pn(i,j)+dy_sP(i,j))*(dy_un(i-1,j)+dy_su(i-1,j))/aPus(i-1,j);
        aNpc(i,j)=(dx_Pe(i,j)+dx_wP(i,j))*(dx_ve(i,j)+dx_wv(i,j))/aPvs(i,j);
        aSpc(i,j)=(dx_Pe(i,j)+dx_wP(i,j))*(dx_ve(i,j-1)+dx_wv(i,j-1))/aPvs(i,j-1);
        aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aNpc(i,j)+aSpc(i,j);
        bpc(i,j)=-(dy_Pn(i,j)+dy_sP(i,j))*us(i,j)+...
            (dy_Pn(i,j)+dy_sP(i,j))*us(i-1,j)+...
            -(dx_Pe(i,j)+dx_wP(i,j))*vs(i,j)+...
            (dx_Pe(i,j)+dx_wP(i,j))*vs(i,j-1)+...
            (dy_Pn(i,j)+dy_sP(i,j))*(dx_Pe(i,j)+dx_wP(i,j))*...
            (-ro(i,j)+ro0(i,j))/dt;
    end
end
i=2; j=2; aPpc(i,j)=aEpc(i,j)+aNpc(i,j); % Esquina 1

i=2; % Nodos del tope izquierdo
for j=3:Ny,aPpc(i,j)=aEpc(i,j)+aNpc(i,j)+aSpc(i,j);end

i=2; j=Ny+1; aPpc(i,j)=aEpc(i,j)+aSpc(i,j); % Esquina  3

for i=3:Nx
    j=2; aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aNpc(i,j); % Nodos del tope inferior
    j=Ny+1; aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aSpc(i,j); % Nodos del tope superior
end

i=Nx+1; j=2; aPpc(i,j)=aWpc(i,j)+aNpc(i,j); % Esquina 2

i=Nx+1; % Nodos del tope derecho
for j=3:Ny, aPpc(i,j)=aWpc(i,j)+aNpc(i,j)+aSpc(i,j); end

i=Nx+1; j=Ny+1; aPpc(i,j)=aWpc(i,j)+aSpc(i,j); % Esquina  4

aWpc(2,:)=zeros(1,Ny+2); aEpc(Nx+1,:)=zeros(1,Ny+2);
aNpc(:,Ny+1)=zeros(Nx+2,1); aSpc(:,2)=zeros(Nx+2,1);
% =======================================================================

% =======================================================================
while (ErrorPE>=TolP)&&(niteraPE<=nmaxP) % -> Esto es si se usa metodo iterativo!
    niteraPE=niteraPE+1; pclast=p;
    
    % info=strcat('#Iteracion = ',num2str(nitera),'_ Error = ',num2str(Error));
    % disp(info)
    for i=2:Nx+1
        for j=2:Ny+1
            p(i,j)=omepc*1/aPpc(i,j)*...
                (aEpc(i,j)*p(i+1,j)+aWpc(i,j)*p(i-1,j)+...
                aNpc(i,j)*p(i,j+1)+aSpc(i,j)*p(i,j-1)+...
                bpc(i,j))+...
                (1-omepc)*p(i,j);
        end
    end
    
    %     ErrorPE=max(max(abs((p-pclast)./p)))*100; 
%     ErrorPE1=max(max(abs((p-pclast)./p)))*100;
%     ErrorPE2=max(max(abs((p-pclast)./mean(mean(abs(p))))))*100;
% %     ErrorPE=max(ErrorPE1,mean([ErrorPE2,ErrorPE1]));
%     ErrorPE=min(ErrorPE1,ErrorPE2*10);
%     
    ErrorPE=(sum((p-pclast).^2)/sum(p.^2))^0.5;
    ErrorPE=ErrorPE*100;
end

disp('Sol. Pressure Equation')
info=strcat('#Iteraciones realizadas = ',num2str(niteraPE),'_ Error = ',num2str(ErrorPE));
disp(info)
% =======================================================================

% i=5; j=5; disp('Coeficientes Nodos Internos')
% aEpc(i,j)/aPpc(i,j)
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/aPpc(i,j)
% 
% disp('Coeficientes Nodos izq'); i=2; j=5;
% aEpc(i,j)/(aEpc(i,j)+aNpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aNpc(i,j)+aSpc(i,j))
% 
% disp('Coeficientes Nodos der'); i=Nx+1; j=5;
% aWpc(i,j)/(aWpc(i,j)+aNpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aWpc(i,j)+aNpc(i,j)+aSpc(i,j))
%         
% disp('Coeficientes Nodos sup'); i=5; j=Ny+1;
% aEpc(i,j)/(aEpc(i,j)+aWpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aWpc(i,j)+aSpc(i,j))
%         
% disp('Coeficientes Nodos infe'); i=5; j=2;
% aEpc(i,j)/(aEpc(i,j)+aWpc(i,j)+aNpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aWpc(i,j)+aNpc(i,j))

% disp('Coeficientes Nodos esquina1'); i=2; j=2;
% aEpc(i,j)/(aEpc(i,j)+aNpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aNpc(i,j))
% 
% disp('Coeficientes Nodos esquina2'); i=Nx+1; j=2;
% aWpc(i,j)/(aWpc(i,j)+aNpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aWpc(i,j)+aNpc(i,j))
% 
% disp('Coeficientes Nodos esquina3'); i=2; j=Ny+1;
% aSpc(i,j)/(aEpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aSpc(i,j))
%     
% disp('Coeficientes Nodos esquina4'); i=Nx+1; j=Ny+1;
% aWpc(i,j)/(aWpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aWpc(i,j)+aSpc(i,j))
    
% =======================================================================