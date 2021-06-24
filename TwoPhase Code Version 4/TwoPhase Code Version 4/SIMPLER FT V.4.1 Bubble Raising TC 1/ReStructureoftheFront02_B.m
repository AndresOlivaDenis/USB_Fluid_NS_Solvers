% ========================================================================
% ReStructure of The Front
% 02 -> Adding & Deleting Points
% 03 -> Reordening
% A -> Interpolacion Lineal, B-> Legendre
% ========================================================================

function [NFront,xFront,yFront,xmapFront,ymapFront]...
    =ReStructureoftheFront02_B(NFront0,xFront,yFront,...
    FrontDistMax,FrontDistMin,...
    xmapFront,ymapFront,FMx,FMy,X1,X2,Y1,Y2)

% ------------------------------------------------------------------------
ilNew=1; xFrontLast=xFront; yFrontLast=yFront;
xmapFrontLast=xmapFront; ymapFrontLast=ymapFront;

% 1ers Point ------------------------------------------------------------ 
ilNew=ilNew+1; il=2;
distanceNFront=((xFrontLast(il)-xFront(ilNew-1))^2+...
    (yFrontLast(il)-yFront(ilNew-1))^2)^0.5;

if (distanceNFront<=FrontDistMax)&&(distanceNFront>=FrontDistMin)
    xmapFront(ilNew)=xmapFrontLast(il); ymapFront(ilNew)=ymapFrontLast(il);
    xFront(ilNew)=X1+(X2-X1)*FMx(xmapFront(ilNew)); 
    yFront(ilNew)=Y1+(Y2-Y1)*FMy(ymapFront(ilNew));
    
elseif (distanceNFront>FrontDistMax) % Add Point
    % Legendre Interpolation
    xmapFront(ilNew)=0.0625*(-xmapFront(NFront0)-xmapFrontLast(il+1)+...
        9*(xmapFront(ilNew-1)+xmapFrontLast(il)));
    ymapFront(ilNew)=0.0625*(-ymapFront(NFront0)-ymapFrontLast(il+1)+...
        9*(ymapFront(ilNew-1)+ymapFrontLast(il)));
    xFront(ilNew)=X1+(X2-X1)*FMx(xmapFront(ilNew)); 
    yFront(ilNew)=Y1+(Y2-Y1)*FMy(ymapFront(ilNew));
    
    ilNew=ilNew+1;
    xmapFront(ilNew)=xmapFrontLast(il); ymapFront(ilNew)=ymapFrontLast(il);
    xFront(ilNew)=X1+(X2-X1)*FMx(xmapFront(ilNew));
    yFront(ilNew)=Y1+(Y2-Y1)*FMy(ymapFront(ilNew));
    
elseif (distanceNFront<FrontDistMin) % Delete Point
    ilNew=ilNew-1;
end
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
for il=3:NFront0+1
    ilNew=ilNew+1;
    distanceNFront=((xFrontLast(il)-xFront(ilNew-1))^2+...
        (yFrontLast(il)-yFront(ilNew-1))^2)^0.5;
    
    if (distanceNFront<=FrontDistMax)&&(distanceNFront>=FrontDistMin)
        xmapFront(ilNew)=xmapFrontLast(il); ymapFront(ilNew)=ymapFrontLast(il);
        xFront(ilNew)=X1+(X2-X1)*FMx(xmapFront(ilNew));
        yFront(ilNew)=Y1+(Y2-Y1)*FMy(ymapFront(ilNew));
        
    elseif (distanceNFront>FrontDistMax) % Add Point ---------------------
        % Legendre Interpolation
        xmapFront(ilNew)=0.0625*(-xmapFront(ilNew-2)-xmapFrontLast(il+1)+...
            9*(xmapFront(ilNew-1)+xmapFrontLast(il)));
        ymapFront(ilNew)=0.0625*(-ymapFront(ilNew-2)-ymapFrontLast(il+1)+...
            9*(ymapFront(ilNew-1)+ymapFrontLast(il)));
        xFront(ilNew)=X1+(X2-X1)*FMx(xmapFront(ilNew));
        yFront(ilNew)=Y1+(Y2-Y1)*FMy(ymapFront(ilNew));
        
        ilNew=ilNew+1; 
        xmapFront(ilNew)=xmapFrontLast(il); ymapFront(ilNew)=ymapFrontLast(il);
        xFront(ilNew)=X1+(X2-X1)*FMx(xmapFront(ilNew));
        yFront(ilNew)=Y1+(Y2-Y1)*FMy(ymapFront(ilNew));
        
    elseif (distanceNFront<FrontDistMin) % Delete Point -----------------
        ilNew=ilNew-1; 
    end
end
NFront=ilNew-1; 
xmapFront(1)=xmapFront(NFront+1); ymapFront(1)=ymapFront(NFront+1);
xFront(1)=xFront(NFront+1); yFront(1)=yFront(NFront+1);
xmapFront(NFront+2)=xmapFront(2); ymapFront(NFront+2)=ymapFront(2);
xFront(NFront+2)=xFront(2); yFront(NFront+2)=yFront(2);
% -----------------------------------------------------------------------

% Cosas de comprobacion -------------------------------------------------

% figure; plot(xFront,yFront,'o',xFront0,yFront0,'o'); legend('New','old')
% 
% Distances=((xFront(2:end)-xFront(1:end-1)).^2+...
%     (yFront(2:end)-yFront(1:end-1)).^2).^0.5;
% 
% Distances0=((xFront0(2:end)-xFront0(1:end-1)).^2+...
%     (yFront0(2:end)-yFront0(1:end-1)).^2).^0.5;
% 
% figure; plot(1:NFront+2-1,Distances,'o',1:NFront0+2-1,Distances0,'o',...
%     1:NFront,FrontDistMax*ones(1,NFront),1:NFront,FrontDistMin*ones(1,NFront))
% legend('New','old','Max','Min')
