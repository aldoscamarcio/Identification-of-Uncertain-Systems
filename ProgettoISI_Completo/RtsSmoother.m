%La fase in avanti è stata fatta con l'EKF, in questo script ci occuperemo
%della fase all'indietro sfruttando le formule viste a teoria
%Dopo aver ottenuto tutte le x_hat e P dal filtro nella fase in avanti, 
% esse serviranno per la fase a ritroso dove si sfrutteranno le formule conosciute a teorie. 
% A questo punto si calcolano le nuove x_hat e P a ritroso.

log_EKF.x_hat_correction = double(log_EKF.x_hat_correction);

%Inizializzazione ultime colonne uguali alle ultime x_hat e P in uscita dal filtro esteso.

log_EKF.x_hat_smoothed(size(log_EKF.x_hat_prediction,1),:) = log_EKF.x_hat_correction(size(log_EKF.x_hat_prediction,1),:);
log_EKF.P_smoothed(size(log_EKF.x_hat_prediction,1),:,:) = [log_EKF.P_correction(size(log_EKF.x_hat_prediction,1),:,1);log_EKF.P_correction(size(log_EKF.x_hat_prediction,1),:,2);log_EKF.P_correction(size(log_EKF.x_hat_prediction,1),:,3);log_EKF.P_correction(size(log_EKF.x_hat_prediction,1),:,4);log_EKF.P_correction(size(log_EKF.x_hat_prediction,1),:,5);log_EKF.P_correction(size(log_EKF.x_hat_prediction,1),:,6)];


   for k = size(log_EKF.x_hat_prediction,1)-1:-1:1 %k va da n-1 fino a 1 con step di 1
    
    % Ck = P_k|k * F_k+1' * P_k+1|k^(-1)
    Ck = [log_EKF.P_correction(k,:,1);log_EKF.P_correction(k,:,2);log_EKF.P_correction(k,:,3);log_EKF.P_correction(k,:,4);log_EKF.P_correction(k,:,5); ...
        log_EKF.P_correction(k,:,6)] ...
        * [log_EKF.F_matrix(k+1,:,1); log_EKF.F_matrix(k+1,:,2);log_EKF.F_matrix(k+1,:,3);log_EKF.F_matrix(k+1,:,4);log_EKF.F_matrix(k+1,:,5); ...
        log_EKF.F_matrix(k+1,:,6);]' / ...
        [log_EKF.P_prediction(k+1,:,1);log_EKF.P_prediction(k+1,:,2);log_EKF.P_prediction(k+1,:,3);log_EKF.P_prediction(k+1,:,4); ...
        log_EKF.P_prediction(k+1,:,5);log_EKF.P_prediction(k+1,:,6)];
    
    % x_k|n = x_k|k + Ck*(x_k+1|n - x_k+1|k)
	log_EKF.x_hat_smoothed(k,:) = log_EKF.x_hat_correction(k,:)' + Ck*(log_EKF.x_hat_smoothed(k+1,:)' - log_EKF.x_hat_prediction(k+1,:)');
    
    
     % P_k|n = P_k|k + Ck*(P_k+1|n - P_k+1|k)*Ck'
    log_EKF.P_smoothed(k,:,:) = [log_EKF.P_correction(k,:,1);log_EKF.P_correction(k,:,2);log_EKF.P_correction(k,:,3); ...
        log_EKF.P_correction(k,:,4);log_EKF.P_correction(k,:,5);log_EKF.P_correction(k,:,6)] ...
         + Ck*([log_EKF.P_smoothed(k+1,:,1);log_EKF.P_smoothed(k+1,:,2);log_EKF.P_smoothed(k+1,:,3); ...
         log_EKF.P_smoothed(k+1,:,4);log_EKF.P_smoothed(k+1,:,5);log_EKF.P_smoothed(k+1,:,6)] - [log_EKF.P_prediction(k+1,:,1); ...
         log_EKF.P_prediction(k+1,:,2);log_EKF.P_prediction(k+1,:,3);log_EKF.P_prediction(k+1,:,4);log_EKF.P_prediction(k+1,:,5); ...
         log_EKF.P_prediction(k+1,:,6)])*Ck';
  end








