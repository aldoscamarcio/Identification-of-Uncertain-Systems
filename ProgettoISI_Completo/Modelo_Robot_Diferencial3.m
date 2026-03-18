% CURSOS ONLINE (MI PLATAFORMA RECOMENDADO)
% 
% Robótica con Matlab y Arduino: Modelo y Simulación: http://bit.ly/2MHWYPb
% Robótica con Matlab y Arduino: Diseño de controladores: http://bit.ly/2SXaDW3
% 
% clc
% clear
% close all


t=0:t_simulation/dt;

u=0.1*ones(1,length(t));
w=0.08*ones(1,length(t));


xr(1)=log_vars.X(1); 
yr(1)=log_vars.Y(1);  
phi(1)=log_vars.theta(1);



pasos=1;  fig=figure;
set(fig,'position',[10 60 980 600]);
axis square; cameratoolbar
axis([-20 20 -20 20 -0.5 1]); grid on
MobileRobot;
M1=MobilePlot(xr(1),yr(1),phi(1));hold on
M2=plot(xr(1),yr(1),'b','LineWidth',2);


for i=1:pasos:length(t)

   
    delete (M1)
    delete (M2)
    M1=MobilePlot(log_vars.X(1,i),log_vars.Y(1,i),log_vars.theta(1,i)); hold on
    M2=plot(log_vars.X(1:i),log_vars.Y(1:i),'blue','LineWidth',2);
  
    pause(dt)

end

