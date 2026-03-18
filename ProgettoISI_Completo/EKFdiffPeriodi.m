%In questo script verrà implementato il filtro di Kalman però cosiderando i
%periodi di campionamento dei sensori diversi
%I dati in uscita da simulink sono di tipo Timeseries

t_simulation=20;
out=sim('isiprojectsimulink2',t_simulation);
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

 ht = matlabFunction(F);
 ht1 = matlabFunction(D);
 ht2= matlabFunction(X_dot);
 ht3 = matlabFunction(H);
 ht4 = matlabFunction(h);
 log_EKF = [];
 
 %Inizializzo vettore delle misure utilizzate che conterrà tutti i tempi
 %delle misure utilizzate essendo univoci

 misure_utilizzate = [];
 misure_psi = 0;
 misure_phi_dot = 0;
 misure_dx = 0;

 %Sto aggiungendo delle misure poichè ci saranno alcuni sensori che avranno
 %più misure e altri sensori che avranno meno misure, queste misure aggiunte 
 % non sono veritiere perciò l'errore alla fine aumenterà.
prova1 = size(log_vars.y_alpha.Time,1) - size(log_vars.phi_dot_sens.Time,1);
prova2 = size(log_vars.y_alpha.Time,1) - size(log_vars.dx_sens.Time,1);
prova3 = size(log_vars.y_alpha.Time,1) - size(log_vars.psi_sens.Time,1);

for t = 1 : 1 : prova1
    log_vars.phi_dot_sens = addsample(log_vars.phi_dot_sens,'Data',0,'Time',10);
end

for t = 1 : 1 : prova2
    log_vars.dx_sens = addsample(log_vars.dx_sens,'Data',0,'Time',10);
end

for t = 1 : 1 : prova3
    log_vars.psi_sens = addsample(log_vars.psi_sens,'Data',0,'Time',10);
end

%Inizializzo matrice dei flag che mi servirà per andare a togliere eventuali 
%colonne o righe alle matrici

flag = zeros(size(log_vars.t,2),4);
  k = 0;
 for t = 0:dt:t_simulation
     k = k + 1;
     T(k) = t;
     [x_hat,P] = predict_EKF(x_hat,P,Q,dt,log_vars.tau_phi(k),log_vars.tau_psi(k),W,F,D,X_dot,ht,ht1,ht2); %predizione
     
      if T(k) < log_vars.psi_sens.Time(k) 
          if ismember(misure_utilizzate.psi,log_vars.psi_sens.time(k-1)) ~= 0 %mi dice se la misura precedente è stata utilizzata poichè va a controllare se il tempo della misura precedente appartiene al vettore delle misure utilizzate 
           flag(k,1) = 1;
           misure_psi = log_vars.psi_sens.time(k-1);
      else
          misure_psi = log_vars.psi_sens.time(k-1);%se non è stata utilizzata allora la utilizzo
      %   log_vars.psi_sens.Data(k) = log_vars.psi_sens.Data(k-1);
          flag(k,1) = 0;
          end
      end

     if t < log_vars.phi_dot_sens.Time(k)
         if ismember(misure_utilizzate.phi_dot,log_vars.phi_dot_sens.time(k-1)) ~= 0   
           flag(k,2) = 1;
           misure_phi_dot = log_vars.phi_dot_sens.time(k-1);
      else
          misure_phi_dot = log_vars.phi_dot_sens.time(k-1);
          flag(k,2) = 0;
       %   log_vars.phi_dot_sens.Data(k) = log_vars.phi_dot_sens.Data(k-1);
         end
     end
         
     if t < log_vars.dx_sens.Time(k) 
         if ismember(misure_utilizzate.dx,log_vars.dx_sens.Time(k-1)) ~= 0  
           flag(k,4) = 1;
           misure_dx= log_vars.dx_sens.Time(k-1);
      else
         misure_dx = log_vars.dx_sens.Time(k-1);
         flag(k,4) = 0;
     %    log_vars.dx_sens.Data(k) = log_vars.dx_sens.Data(k-1);
      end
     end
  
    [x_hat,P,e] = correct_EKF(x_hat, P, R,log_vars.dx_sens.Data(k),log_vars.phi_dot_sens.Data(k),log_vars.psi_sens.Data(k),log_vars.y_alpha.Data(k),ht3,ht4,Mk,h,H,flag,k); %correzione
  log_EKF.x_hat(k,1:6) = x_hat;
  log_EKF.P(k,:,:) = P;
  misure_utilizzate.psi(k) = misure_psi;
  misure_utilizzate.phi_dot(k) = misure_phi_dot;
  misure_utilizzate.dx(k) = misure_dx;

 end
 
 
 %predizione
 
 function  [x_hat,P] = predict_EKF(x_hat,P,Q,dt,tau_phi,tau_psi,W,F,D,X_dot,ht,ht1,ht2)

 F = feval(ht,x_hat(3),x_hat(4),x_hat(5),x_hat(6),tau_phi,tau_psi,0,0);
 D = feval(ht1,x_hat(5));
 X_dot = feval(ht2,x_hat(3),x_hat(4),x_hat(5),x_hat(6),tau_phi,tau_psi,0,0);
  x_hat = x_hat + X_dot * dt ;
    P = F*P*F' + D*Q*D';
 end

%Correzione
%Quando flag = 0 le misure le prende invece quando flag = 1 le misure non
%le prende, se la misura è stata utilizzata oppure non è arrivata allora le
%matrici perderanno righe o/e colonne come posssiamo notare dagli if
 function [x_hat,P,e] = correct_EKF(x_hat, P, R,dx_sens,phi_dot_sens,psi_sens,y_alpha,ht3,ht4,Mk,h,H,flag,k)
 H = feval(ht3,x_hat(1),x_hat(2));
 h = feval(ht4,x_hat(1),x_hat(2),x_hat(4),x_hat(5),0,0,0,0);
 misure = [psi_sens ; phi_dot_sens ; dx_sens ; y_alpha];
 if flag(k,1) == 1
     H(1,:) = [];
     h(1) = [];
     R(1,:) = [];
     R(:,1) = [];
     Mk(1,:) = [];
     Mk(:,1) = [];
     misure(1) = [];
 end
  if flag(k,2) == 1
     H(2,:) = [];
     h(2) = [];
     R(2,:) = [];
     R(:,2) = [];
     Mk(2,:) = [];
     Mk(:,2) = [];
     misure(2) = [];
  end
  if flag(k,3) == 1
     H(3,:) = [];
     h(3) = [];
     R(3,:) = [];
     R(:,3) = [];
     Mk(3,:) = [];
     Mk(:,3) = [];
     misure(3) = [];
 end
 e = misure - h';
 S = vpa(H*P*H' + Mk*R*Mk',3);
 G = vpa(P*H'*inv(S),3);
 x_hat = x_hat + G*e
 P = (eye(6) - G*H)*P*(eye(6) - G*H)' + G*Mk*R*Mk'*G'
 end
 
