clc; clear all; close all
NInter=60; 

% It = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
load('SIMPLER03_h=160_It=0_t0.mat')
% hold on
% plot(xFront(1:NFront+1),yFront(1:NFront+1),'k','linewidth',1.5);pause(0.01) 
contourf(xP,yP,MarkF,1);
mapAndres=[0.65,0.75,0.90;
    0.99,0.99,0.99];
colormap(mapAndres)
axis equal; axis([0,1,0,2]);
F(1) = getframe;

% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=3840; dt=dt*2;
for iinte=1:NInter
    iinte
    it=((iinte)*Nt/NInter);
    infoLoad=strcat('SIMPLER_CN_h=80_It=',num2str(it),'_t',num2str(it*dt),'.mat');
    load(infoLoad);
    
    contourf(xP,yP,MarkF,1);
    mapAndres=[0.65,0.75,0.90;
        0.99,0.99,0.99];
    colormap(mapAndres)
    axis equal; axis([0,1,0,2]);


    F(iinte+1) = getframe;
end

writerObj = VideoWriter('Bubble_Rising_Case1.avi'); writerObj.FrameRate = 20; 
open(writerObj); writeVideo(writerObj,F); close(writerObj);
