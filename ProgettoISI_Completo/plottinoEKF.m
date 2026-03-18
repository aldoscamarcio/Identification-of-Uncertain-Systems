
log_vars.x = double(log_EKF.x_hat(:,1)');
log_vars.y = double(log_EKF.x_hat(:,2)');
log_vars.thetas = double(log_EKF.x_hat(:,3)');
log_vars.phi_dots=double(log_EKF.x_hat(:,4)');
log_vars.psis=double(log_EKF.x_hat(:,5)');
log_vars.psi_dots=double(log_EKF.x_hat(:,6)');

figure(1)
clf
grid on
box on
xlabel('tempo [s]')
ylabel('[m]')
title('position estimation error [m]')
hold on
plot(log_vars.t, vecnorm([log_vars.X ; log_vars.Y] - [log_vars.x ; log_vars.y]));


figure(2)
x = log_vars.t;


tiledlayout(3,2)
nexttile
y1= log_vars.X;
y2 = log_vars.x;
hold on
plot(x,y1)
plot(x,y2)
xlabel('tempo [s]')
ylabel('x [m]')
legend('posizione x vera','posizione x predetta')
hold off

nexttile
y7=log_vars.phi_dot;
y8=log_vars.phi_dots;
hold on
plot(x,y7)
plot(x,y8)
xlabel('tempo [s]')
ylabel('phi_dot [rad/s]')
legend('posizione phi_dot vera','posizione phi_dot predetta')
hold off

nexttile
y3=log_vars.Y;
y4=log_vars.y;
hold on
plot(x,y3)
plot(x,y4)
xlabel('tempo [s]')
ylabel('y [m]')
legend('posizione y vera','posizione y predetta')
hold off

nexttile
y9=log_vars.psi;
y10=log_vars.psis;
hold on
plot(x,y9)
plot(x,y10)
xlabel('tempo [s]')
ylabel('psi [rad]')
legend('posizione psi vera','posizione psi predetta')
hold off

nexttile
y5=log_vars.theta;
y6=log_vars.thetas;
hold on
plot(x,y5)
plot(x,y6)
xlabel('tempo [s]')
ylabel('theta [rad]')
legend('posizione theta vera','posizione theta predetta')
hold off


nexttile
y11=log_vars.psi_dot;
y12=log_vars.psi_dots;
hold on
plot(x,y11)
plot(x,y12)
xlabel('tempo [s]')
ylabel('psi_dot [rad/s]')
legend('posizione psi_dot vera','posizione psi_dot predetta')
hold off

figure(3)


tiledlayout(3,2)
nexttile
y1= log_vars.X;
y2 = log_vars.x;
hold on
plot(x,y1-y2)
xlabel('tempo [s]')
ylabel('x [m]')
title('X vera - X predetta')
hold off

nexttile
y7=log_vars.phi_dot;
y8=log_vars.phi_dots;
hold on
plot(x,y7-y8)
xlabel('tempo [s]')
ylabel('phi [rad]')
title('phi vero - phi predetto')
hold off

nexttile
y3=log_vars.Y;
y4=log_vars.y;
hold on
plot(x,y3-y4)
xlabel('tempo [s]')
ylabel('y [m]')
title('Y vera - Y predetta')
hold off

nexttile
y9=log_vars.psi;
y10=log_vars.psis;
hold on
plot(x,y9-y10)
xlabel('tempo [s]')
ylabel('psi [rad]')
title('psi vero - psi predetto')
hold off

nexttile
y5=log_vars.theta;
y6=log_vars.thetas;
hold on
plot(x,y5-y6)
xlabel('tempo [s]')
ylabel('theta [rad]')
title('theta vero - theta predetto')
hold off


nexttile
y11=log_vars.psi_dot;
y12=log_vars.psi_dots;
hold on
plot(x,y11-y12)
title('psi_dot vero - psi_dot predetto')
xlabel('tempo [s]')
ylabel('psi_dot [rad/s]')
hold off