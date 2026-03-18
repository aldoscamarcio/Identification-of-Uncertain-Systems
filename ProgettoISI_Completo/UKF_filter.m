%Filtro UKF

t_simulation=10;
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

 ht= matlabFunction(X_dot);
 ht1 = matlabFunction(h);

 k = 0;
 for t = 0:dt:t_simulation
     k = k + 1;
 
      [x_hat,P] = predict_UKF(x_hat,P,Q,dt,log_vars.tau_phi(k),log_vars.tau_psi(k),ht); %predizione
      [x_hat,P,e] = correct_UKF(x_hat, P, R,log_vars.dx_sens(k),log_vars.phi_dot_sens(k),log_vars.psi_sens(k),log_vars.y_alpha(k),ht1); %correzione
   log_UKF.x_hat(k,1:6) = x_hat;
   log_UKF.P(k,1:6,1:6) = P;
   log_UKF.e(k,1:4) = e;
 
  end
 
 
 function [x_hat,P] = predict_UKF(x_hat,P,Q,dt,tau_phi,tau_psi,ht)
 %In input alla trasformata Unscented devono andare le medie, le matrici di
 %covarianze e il modello non lineare. Le medie dei rumori sono nulle
    [x_hat,P] = UnscentedTransform_F_Function([x_hat; 0; 0], blkdiag(P, Q), dt, tau_phi, tau_psi,ht);
 end

function [x_hat,P,e] = correct_UKF(x_hat, P, R,dx_sens,phi_dot_sens,psi_sens,y_alpha,ht1)
 [hat_y, S, Pxy] = UnscentedTransform_H_Function(x_hat, P,ht1);
 %essendo il rumore additivo non va nelle medie in input alla trasformata
 %Unscented poichè verrrà sommato dopo
    y = [psi_sens; phi_dot_sens; dx_sens ;y_alpha];
    S = S + R; %sommiamo la matrice delle varianze dei rumori di misura
    e = y - hat_y;
   e(1) = atan2(sin(e(1)),cos(e(1)));
   e(2) = atan2(sin(e(2)), cos(e(2)));
   e(4) = atan2(sin(e(4)),cos(e(4))); %la differenza di angoli deve essere compresa tra -pi e pi 
    L = Pxy*S^(-1);
    
    x_hat = x_hat + L*e;
    P = P - L*S*L';
end

 function [propagated_mean, propagated_cov] = UnscentedTransform_F_Function(misura_priori, cov_priori, dt, tau_phi, tau_psi,ht)
    % definizione parametri
    alpha = 1;            
    beta = 2;             % ottimo nel caso gaussiano
    kappa = 0;            % lambda = 0 (peso centrale dei sigma point)
    
    enne = size(misura_priori,1); %n è il vettore iniziale composto dalle variabili di stato e la media dei rumori che è nulla(nel nostro caso = 8)
    
    lambda = alpha^2*(enne + kappa) - enne;
    
    % Calcolo dei pesi
    w0 = lambda/(lambda + enne);
    wc(1) = w0 + 1 - alpha^2 + beta;
    wc(2:2*enne+1) = 1/2/(enne+lambda);
    wm = wc;
    wm(1) = w0;
    
    % Fattorizzazione matrice delle covarianze
    % SVD
    [U,S,~] = svd(cov_priori);
    GAMMA = U*S^(1/2); %le colonne ci indicano le direzioni per andare a generare i sigma point
    
    % generazione degli 2n+1 sigma point (saranno 17 sigma point)
    sigma_points = misura_priori; % il primo sigma point centrale coicide con la media
    for i = 1:size(GAMMA,2) %definisco gli altri 2n sigma point simmetrici centrati sulla media
        sigma_points(:,i+1)      = misura_priori + sqrt(enne+lambda)*GAMMA(:,i);
        sigma_points(:,i+1+enne) = misura_priori - sqrt(enne+lambda)*GAMMA(:,i);
    end
    
    %Per ognuno di questi sigma-point applichiamo la funzione f in modo
    %tale da ottenere i sigma-point propagati con il modello del sistema
    for i = 1:size(sigma_points,2)
        X_dot = feval(ht,sigma_points(3,i),sigma_points(4,i),sigma_points(5,i),sigma_points(6,i),tau_phi,tau_psi,sigma_points(7,i),sigma_points(8,i));
        propagated_sigma_points(1:6,i) = sigma_points(1:6,i) + dt*X_dot;
        
    end

    %Media pesata dai wn (fare attenzione alla media di angoli: si sfrutta
    %l'atan2
   propagated_mean(1,:) = propagated_sigma_points(1,:)*wm';
   propagated_mean(2,:) = propagated_sigma_points(2,:)*wm';
   propagated_mean(3,:) = atan2(sin(propagated_sigma_points(3,:)*wm'),cos(propagated_sigma_points(3,:)*wm'));
   propagated_mean(4,:) = atan2(sin(propagated_sigma_points(4,:)*wm'),cos(propagated_sigma_points(4,:)*wm'));
   propagated_mean(6,:) = atan2(sin(propagated_sigma_points(6,:)*wm'),cos(propagated_sigma_points(6,:)*wm'));
   propagated_mean(5,:) = atan2(sin(propagated_sigma_points(5,:)*wm'),cos(propagated_sigma_points(5,:)*wm'));

     % Tilde propagato (differenza rispetto alla media che ci servirà per
     % la matrice di covarianza campionaria)
     propagated_tilde(1,:) = propagated_sigma_points(1,:) - propagated_mean(1,:);
     propagated_tilde(2,:) = propagated_sigma_points(2,:) - propagated_mean(2,:);
     propagated_tilde(3,:) = atan2(sin(propagated_sigma_points(3,:) - propagated_mean(3,:)),cos(propagated_sigma_points(3,:) - propagated_mean(3,:)));
     propagated_tilde(4,:) = atan2(sin(propagated_sigma_points(4,:) - propagated_mean(4,:)),cos(propagated_sigma_points(4,:) - propagated_mean(4,:)));
     propagated_tilde(5,:) = atan2(sin(propagated_sigma_points(5,:) - propagated_mean(5,:)),cos(propagated_sigma_points(5,:) - propagated_mean(5,:)));
     propagated_tilde(6,:) = atan2(sin(propagated_sigma_points(6,:) - propagated_mean(6,:)),cos(propagated_sigma_points(6,:) - propagated_mean(6,:)));


    % Momento incrociato non ci serve
    propagated_cov = zeros(6,6);
    for i = 1:size(sigma_points,2)
        propagated_cov = propagated_cov + wc(i)*propagated_tilde(:,i)*propagated_tilde(:,i)';
    end
 end



function [virtual_measurements_mean, virtual_measurements_cov, cross_cov] = UnscentedTransform_H_Function(prior_mean, prior_cov,ht1)
   
    alpha = 1;            
    beta = 2;          
    kappa = 0;          
    
    enne = size(prior_mean,1); % n = 6
    
    lambda = alpha^2*(enne + kappa) - enne;
    
    w0 = lambda/(lambda + enne);
    wc(1) = w0 + 1 - alpha^2 + beta;
    wc(2:2*enne+1) = 1/2/(enne+lambda);
    
    wm = wc;
    wm(1) = w0;
    
    [U,S,~] = svd(prior_cov);
    GAMMA = U*S^(1/2);
    
    sigma_points = prior_mean;
    for i = 1:size(GAMMA,2)
        sigma_points(:,i+1)      = prior_mean + sqrt(enne+lambda)*GAMMA(:,i);
        sigma_points(:,i+1+enne) = prior_mean - sqrt(enne+lambda)*GAMMA(:,i);
    end
    
    % In questo caso bisogna propagare i sigma points con il modello di
    % osservazione
    for i = 1:size(sigma_points,2)
        propagated_sigma_points(1:4,i) = feval(ht1,sigma_points(1,i),sigma_points(2,i),sigma_points(4,i),sigma_points(5,i),0,0,0,0);
    end
    
   
    virtual_measurements_mean(1,:) = atan2(sin(propagated_sigma_points(1,:)*wm'),cos(propagated_sigma_points(1,:)*wm'));
    virtual_measurements_mean(2,:) = atan2(sin(propagated_sigma_points(2,:)*wm'),cos(propagated_sigma_points(2,:)*wm'));
    virtual_measurements_mean(3,:) = propagated_sigma_points(3,:)*wm';
    virtual_measurements_mean(4,:) = atan2(sin(propagated_sigma_points(4,:)*wm'),cos(propagated_sigma_points(4,:)*wm'));

    virtual_measurements_tilde(1,:) = atan2(sin(propagated_sigma_points(1,:) - virtual_measurements_mean(1,:)),cos(propagated_sigma_points(1,:) - virtual_measurements_mean(1,:)));
    virtual_measurements_tilde(2,:) = atan2(sin(propagated_sigma_points(2,:) - virtual_measurements_mean(2,:)),cos(propagated_sigma_points(2,:) - virtual_measurements_mean(2,:)));
    virtual_measurements_tilde(3,:) = propagated_sigma_points(3,:) - virtual_measurements_mean(3,:);
    virtual_measurements_tilde(4,:) = atan2(sin(propagated_sigma_points(4,:) - virtual_measurements_mean(4,:)),cos(propagated_sigma_points(4,:) - virtual_measurements_mean(4,:)));
    
    %In questo caso il momento incrociato ci serve
    virtual_measurements_cov = zeros(4,4);
    cross_cov = zeros(6,4);
    for i = 1:size(sigma_points,2)
        virtual_measurements_cov = virtual_measurements_cov + wc(i)*virtual_measurements_tilde(:,i)*virtual_measurements_tilde(:,i)';
        cross_cov = cross_cov + wc(i)*(sigma_points(:,i) - prior_mean)*virtual_measurements_tilde(:,i)';
    end
end



 









