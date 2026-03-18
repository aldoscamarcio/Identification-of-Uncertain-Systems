%Abbiamo i dati che fuoriescono dai sensori che noi preendiamo da Simulink.
%Questi dati li prendiamo e li mettiamo in un'unica struttura che ci
%servirà per inizializzare l'EKF


t_simulation=20;
out=sim('isiprojectsimulink',t_simulation);
log_vars = [];
log_vars.X = (out.X)';
log_vars.Y = (out.Y)';
log_vars.dx_sens = (out.dx_sens)';
log_vars.phi_dot = (out.phi_dot)';
log_vars.phi_dot_sens = (out.phi_dot_sens)';
log_vars.psi = (out.psi)';
log_vars.psi_dot= (out.psi_dot)';
log_vars.psi_sens = (out.psi_sens)';
log_vars.tau_phi = (out.tau_phi)';
log_vars.tau_psi = (out.tau_psi)';
log_vars.theta = (out.theta)';
log_vars.y_alpha = (out.y_alpha)';
log_vars.x_hat = x_hat;
log_vars.P = P;

 
k=0;
 for p=0:dt:t_simulation
     k=k+1;
 log_vars.t(k)=p;
end
 
 log_vars.dt = t_sample;
 log_vars.t_simulation = t_simulation;

 %Inizializzazione matlabFunction per assegnare valori alle variabili
 %simboliche
 ht = matlabFunction(F,'Optimize',false);
 ht1 = matlabFunction(D,'Optimize',false);
 ht2= matlabFunction(X_dot,'Optimize',false);
 ht3 = matlabFunction(H,'Optimize',false);
 ht4 = matlabFunction(h,'Optimize',false);
 log_EKF = [];

 %Dopo aver creato la struttura proseguo col costrire il ciclo for in cui
 %inserirò sia la parte di predizione che di correzione
 
 k = 0;
 for t = 0:dt:t_simulation
     k = k + 1;
 
     [x_hat,P] = predict_EKF(x_hat,P,Q,dt,log_vars.tau_phi(k),log_vars.tau_psi(k),W,F,D,X_dot,ht,ht1,ht2); %predizione
     [x_hat,P,e] = correct_EKF(x_hat, P, R,log_vars.dx_sens(k),log_vars.phi_dot_sens(k),log_vars.psi_sens(k),log_vars.y_alpha(k),ht3,ht4,Mk,h,H); %correzione
  log_EKF.x_hat(k,1:6) = x_hat;
  log_EKF.P(k,:,:) = P;
  log_EKF.e(k,1:4) = e;

 end
 
 
 %Predizione
 
 function  [x_hat,P] = predict_EKF(x_hat,P,Q,dt,tau_phi,tau_psi,W,F,D,X_dot,ht,ht1,ht2)
%Attraverso funzione feval e MatlabFunction assegno valori alle variabili
%simboliche
%Calcolo F, D e X_dot in x predetto al passo k/k e rumore nullo quindi 0,0
 F = feval(ht,x_hat(3),x_hat(4),x_hat(5),x_hat(6),tau_phi,tau_psi,0,0);
 D = feval(ht1,x_hat(5));
 X_dot = feval(ht2,x_hat(3),x_hat(4),x_hat(5),x_hat(6),tau_phi,tau_psi,0,0);

  x_hat = x_hat + X_dot * dt ; %discretizzazione
    P = F*P*F' + D*Q*D'; %Calolo matrice P al passo k+1/k
 end

%Correzione

 function [x_hat,P,e] = correct_EKF(x_hat, P, R,dx_sens,phi_dot_sens,psi_sens,y_alpha,ht3,ht4,Mk,h,H)
 %Calcolo H e h in x predetto di k / k-1 e rumore nullo
 H = feval(ht3,x_hat(1),x_hat(2));
 h = feval(ht4,x_hat(1),x_hat(2),x_hat(4),x_hat(5),0,0,0,0);
 e = [psi_sens;phi_dot_sens;dx_sens;y_alpha] - h';
 e(1) = atan2(sin(e(1)), cos(e(2)));
 e(2) = atan2(sin(e(2)), cos(e(2)));
 e(4) = atan2(sin(e(4)),cos(e(4))); %la differenza di angoli deve essere compresa tra -pi e pi quindi sfrutto atan2
 
 %Applico stima BLUE
 S = vpa(H*P*H' + Mk*R*Mk',3); %comando vpa per effettuare rapporti poichè variabili simboliche in forma razionale
 G = vpa(P*H'*inv(S),3); %G = guadagno di correzione
 
 x_hat = x_hat + G*e
 P = (eye(6) - G*H)*P*(eye(6) - G*H)' + G*Mk*R*Mk'*G'
  end