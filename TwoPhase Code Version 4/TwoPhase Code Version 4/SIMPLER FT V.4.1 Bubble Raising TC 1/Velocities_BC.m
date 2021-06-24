% =======================================================================
% Velocities Boundary Condition
% =======================================================================
function [u,v] = Velocities_BC(uCBeast, uCBwest, vCBsouth, vCBnorth,...
    uCBsouth, uCBnorth, vCBwest, vCBeast, CuCB, CvCB,...
    dx_uE, dx_Wu, dy_vN, dy_Sv, dy_uN, dx_Wv, dy_Su, dx_vE, ...
    Nx,Ny,u,v)

% =======================================================================
% Normal Velocities -----------------------------------------------------
u(1,:)=uCBwest*(1-CuCB(1))  +(u(2,:)-dx_uE(1,:)*uCBwest)*CuCB(1);
u(Nx+1,:)=uCBeast*(1-CuCB(3))  +(u(Nx,:)+dx_Wu(Nx+1,:)*uCBeast)*CuCB(3);
v(:,1)=vCBsouth*(1-CvCB(4))  +(v(:,2)-dy_vN(:,1)*vCBsouth)*CvCB(4);
v(:,Ny+1)=vCBnorth*(1-CvCB(2))  +(v(:,Ny)+dy_Sv(:,Ny+1)*vCBnorth)*CvCB(2);
% ------------------------------------------------------------------------
% Tangential Velocities --------------------------------------------------
u(:,1)=(2*uCBsouth-u(:,2))*(1-CuCB(4))+(u(:,2)-dy_uN(:,2)*uCBsouth)*CuCB(4);
u(:,Ny+2)=(2*uCBnorth-u(:,Ny+1))*(1-CuCB(2))+(u(:,Ny+1)+dy_Su(:,Ny+2)*uCBnorth)*CuCB(2);
v(1,:)=(2*vCBwest-v(2,:))*(1-CvCB(1))+(v(2,:)-dx_vE(1,:)*vCBwest)*CvCB(1);
v(Nx+2,:)=(2*vCBeast-v(Nx+1,:))*(1-CvCB(3))+(v(Nx+1,:)+dx_Wv(Nx+2,:)*vCBeast)*CvCB(3);
% =======================================================================
