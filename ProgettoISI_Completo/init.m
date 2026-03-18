clear all
close all
clc

x_start = [0 0 0]';

x_meas = [0 0 0]';

L = 1.5; %m 
rp = 0.2; %m 
l=0.8; %m 
IPy=0.045;
IPz=3.5;
Mv=80; %kg
mp=2; %kg
ma=2; %kg
IAy=0.025;
ra=0.20; %m
d=0.4; %m
IGz=7.5;
IAz=7;
a = (Mv + mp + 2*ma + (2*IAy/(ra^2)));
b = 1/2*(Mv*(l^2) + mp*(L^2) + 2*ma*(d^2) +IGz + IPz + 2*IAz + 2*IAy*((d/ra)^2));
c = Mv*l + mp*L;


t_sample = 0.1;
sensor_var = 0;

%Costruzione modello su Matlab

%Vettore di stato -> x1 = x ; x2 = y ; x3 = theta ; x4 = phi_dot ; x5 = psi
% x6 = psi_dot

% tau = [tau_phi tau_psi] %vettore degli ingressi
tau= sym('tau',[ 2 1 ],'real');

% Vettore di stato in forma simbolica

X = sym('X',[6 1 ],'real');


%Inizializzazione variabili ( non conoscendo la posizione iniziale la
%ipotizzo random

limits = 10;

x_0 = limits * rand(1); % posizione iniziale x 

y_0 = limits * rand(1); %posizione iniziale di y 

theta_0 = 2*pi*rand(1); %imposto theta angolo che varia tra 0 e 360 gradi dovremmo inserire 0,1 intervallo

phi_dot_0 = 0; %suppondendo che la velocità iniziale della ruota posteriore è nulla

psi_dot_0 = 0; %supponendo che la velocità iniziale sterzante nulla

psi_0 = 0; %2*pi*rand(1); %supponogo random la posizione dello sterzo iniziale anche qui 0,1


%Non sappiamo dove si trova l'AGV ma stiamo ipotizzando che si trovi in un
%campo 10mx10m con disco di incertezza 0.3 inoltre suppongo incerto anche
%l'angolo di sterzo invece le velocità iniziali nulle

std_dev_init_x_y = 0.5 ;%supponendo disco di incertezza 0.5
std_dev_init_angoli = 0.1 ; %supponendo varianza di 0.6 gradi
std_dev_init_vel_ang = 0.05 ; %supponendo varianza di 0.15 rad/s

xr_0 = x_0 + std_dev_init_x_y*randn(1); %posizione x iniziale sporcata  
yr_0 = y_0 + std_dev_init_x_y*randn(1); %posizione y iniziale sporcata
thetar_0 = theta_0 + std_dev_init_angoli* randn(1); %posizione theta iniziale sporcata
psir_0 = psi_0 + std_dev_init_angoli*randn(1); %posizione psi iniziale sporcata

x_hat = [xr_0 yr_0 thetar_0 0 psir_0 0]';
P = blkdiag(std_dev_init_x_y , std_dev_init_x_y, std_dev_init_angoli , std_dev_init_vel_ang ,std_dev_init_angoli ,std_dev_init_vel_ang)^2; 
%suppongo che la matrice delle covarianze abbia queste varianze(che posso modificare) il che significa 
%che le misure si discostano di tot dal valore vero.


        %Caratteristiche sensori
V = sym('v',[4 1],'real'); %rumore aggiunto ai sensori ci sarà utile nella h

dev_psi = 0.0017; %circa 0,1 gradi.
        % La precisione dell'encoder è una misura dell'errore tra il valore letto dall'encoder 
        % e il valore fisico effettivo misurato. 
        % La precisione dell'encoder è misurata in minuti d'arco o secondi d'arco con 20 minuti d'arco (0,3 gradi) 
        % o meglio generalmente considerato un encoder ad alta precisione 
        % con alcuni dispositivi di precisione dell'ordine di 5 secondi d'arco (0,0014 gradi).
dev_phi_dot = 0.01; %rad/s  

dev_dx = 0.02; %m
dev_yalpha= 0.0017; %rad

%Anche gli ingressi possono essere affetti da rumore come ipotizzato anche nel modello simulink
W = sym('w',[2 1],'real'); %w indica il rumore in ingresso, ci sarà utile nella f

dev_w1 = 0.04;
dev_w2 = 0.04;

%Calcolo di h. Come sappiamo l'h si riferisce al modello di osservazione
%quindi al modello del sensore. In questo caso è molto semplice calcolarlo
%poichè dobbiamo solo considerare le 4 uscite dei sensori che esprimiamo
%dalle variabili di stato ed inoltre aggiungiamo rumore poichè sappiamo che
%alle misure aggiungiamo rumore;

y_alpha = atan2(X(2),X(1));
h = [X(5) + V(1) , X(4) + V(2) , X(1) + V(3) ,  y_alpha + V(4) ];

% con x1 = x ; x2 = y ; x3 = theta ; x4 = phi_dot ; x5 = psi
% x6 = psi_dot


%Ora definiamo la f (modello del sistema) la quale ci sarà utile nell'implementazione del filtro
%di Kalman

%La f è associata al modello del sistema e quindi all'equazioni della
%dinamica e della cinematica

%Quello che bisogna fare per scrivere f è per prima cosa come fatto nello
%schema simulink scrivere tutto nella forma B*q_ddot + C(q,q_dot) = tau
%sostituendo le variabili di stato costruite come variabili simpoliche
%(utile per caloclare jacobiano dopo). Nella f inoltre dobbiamo anche tener
%conto del rumore del segnale di ingresso ecco perchè alla tao sommiamo W
%che è il rumore( tutte variabili simboliche poichè dopo costruiremo le
%matrici a partire da queste con le derivate).(rendendo così tau=tau_w) 
% Come sappiamo per il calcolo della f bisogna calcolare x_dot (vettore di stato), discretizzare 
% e formula X_k+1 = X_k + dt* X_dot con X_k condizione iniziale

B = [a*rp^2*cos(X(5))^2+b*(rp/L)^2*sin(X(5))^2+IPy, -IPz*(rp/L)*sin(X(5)); -IPz*(rp/L)*sin(X(5)), IPz];


C = [(b*(rp/L)^2-a*rp^2)*X(4)*X(6)*sin(X(5))*cos(X(5)); -IPz*(rp/L)*X(4)*X(6)*cos(X(5))];
tau_w = [tau(1)+ W(1) tau(2)+W(2)] ;

M = (inv(B)*((tau_w)' - C ));

%Calcolo X_dot

x1_dot=rp*X(4)*cos(X(3))*cos(X(5));
x2_dot=rp*X(4)*sin(X(3))*cos(X(5));
x3_dot=-(rp*X(4)*sin(X(5)))/L;
x4_dot=M(1); %phi_ddot
x5_dot=X(6); %psi_dot
x6_dot=M(2); %psi_ddot

X_dot=[x1_dot x2_dot x3_dot x4_dot x5_dot x6_dot]';

%Ora si passa a implementare il filtro di Kalman esteso con annesse matrici
%F Q R D H M, Sappiamo dalla teoria che il filtro di Kalman esteso è per
%sistemi non lineari e pertanto si dovrà effettuare una linearizzazion al
%primo termine per poi discretizzare e derivare.


%Si sono definite le 6 matrici d'interesse R Q F H D M, in particolare, F D tramite lo Jacobiano della
%funzione f definita :f=X+dt*X_dot come già prima enunciato f è associata al modello del sistema 
% e quindi all'equazioni della dinamica e della cinematica e H M tramite lo
% jacobiano della funzione h (si riferisce al modello di osservazione
%quindi al modello del sensore sporcato dal rumore additivo. 


dt=t_sample;% tempo di campionamento
f=X+dt*X_dot;

R=blkdiag(dev_psi^2,dev_phi_dot^2,dev_dx^2,dev_yalpha^2); %R è la matrice in cui sono presenti sulla diagonale le varianze (deviazioni al quadrato) dei sensori

Q=blkdiag(dev_w1^2,dev_w2^2); %Q è la matrice in cui sono presenti sulla diagonale le varianze (deviazioni al quadrato) degli ingressi

F=jacobian(f,X); %la F si calcola derivando la f per gli stati quindi sarà una matrice 6x6

D=jacobian(f,W); % la D si calcola derivando la f per il rumore applicato agli ingressi quindi 
% uniche due righe diverse da zero saranno la 4 e la 6

H=jacobian(h,X); %la H si calcola derivando la h per gli stati quindi sarà una matrice 4x6

Mk=jacobian(h,V); % la M si calcola derivando la h per gli ingressi rumorosi, siccome gli ingressi rumorosi sono additivi allora la matrice sarà identità 4x4


%Adesso implemento il filtro con i diversi periodi di campionamento dei sensori 
t_sample_psi = 0.10012;
t_sample_phi_dot = 0.10010;
t_sample_dx = 0.10010;
t_sample_y_alpha = 0.1;









