clc
clear
close all
if 1
% return
addpath(genpath(pwd))
% return
load('LIMITI_per_il_gradiente_proiettato.mat')
% Definizione di dati
DATA.tipo_u = 3;         DATA.tipo_mesh = 10; % 13; ;      
DATA.tipo_base = 'hex'; 'corone circolari';'rettangoli'; 
DATA.raggio_ext = 0.6;
DATA.raggio_osta = 1/6;     DATA.raggio_mant = DATA.raggio_osta * 1.575;
% if ~exist("model",'var')
%     model = mphload('Mesh_concavo.mph');
% end
% DATA.model = model;
% DATA.fem = 'P1';
DATA = Dati_cloak(DATA);
DATA.hh_mesh = 1500/max(DATA.omega/2/pi) / 5.5*2;
DATA.mesh = 'mesh_cloak';
% Creo MESH e FEM
[MESH, FE_SPACE] = crea_MESH_e_FEM(DATA);
% disegno MESH e normali
disegno_mesh(DATA,MESH)
% Definizione basi
DATA.lato_hex = 2*5*10^(-3);    % 5 mm
DATA = def_basi(DATA,MESH);
end
% load('Concavo soft.mat')
% Info
fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Numero di vertici       = %d \n',MESH.numVertices);
fprintf(' * Numero di nodi          = %d \n',MESH.numNodes);
fprintf(' * Numero di elementi      = %d \n',MESH.numElem);
fprintf(' * Numero di contorni      = %d \n',MESH.n_contorni);
fprintf(' * Numero di domini        = %d \n',MESH.n_domini);
fprintf(' * Numero di basi di controllo  = %d \n\n',DATA.n_basi);

%% Generazione matrici dominio
quali = {'A_0','A_j', 'B_0','B_j', 'C_e','D  ','E  ','B_a'};
        [ A_0,  A_j,   B_0,  B_j,   C_e,  D,    E,    B_a] = Genera_matrici_parfor(DATA, MESH, FE_SPACE,quali);
% Dirichlet
gamma_in = trova_gamma_in(DATA,MESH,'s');
%% Inizializziamo

% Y_inc nei vertici
Y_inc = zeros(MESH.numNodes,DATA.n_frq);
for ff = 1:DATA.n_frq
    Y_inc(:,ff) = DATA.y_inc_fun({MESH.nodes(1,:), MESH.nodes(2,:)}, DATA.k_amb(ff),DATA.direz(ff)).';
end

% Inizializzo i vettori Y, P e U
Y = zeros(MESH.numNodes,DATA.n_frq);
P = Y;
Y(gamma_in,:) = -Y_inc(gamma_in,:);
U = zeros(DATA.n_basi,1)+0.1*0;
U_old = U;
V = U;
V_old = V;

dJ_U = zeros(DATA.n_basi,1);
dJ_V = dJ_U; 
% return
%% Cose per il costo
[A_costo,N_costo] = calcola_vicini(DATA);

%% Iteriamo
% Costanti di iterazione
tau =  15*0;
MaxIt =  3000;
lmb_u = 0*1e-6;  % costo di U -> Bulk
lmb_v = 0*1e-6;  % costo di V -> rho

tol = 1e-3;
disegna_ogni = 100;
uscita = false;
figure(4)
set(4,'Position',[400  1  1136 780])

ii = 1; if tau == 0,  MaxIt = 2;     end
J = zeros(MaxIt,1); J_costo = J;
J(1) = 1e10;
dJ_norma = zeros(MaxIt,2);
tempo_medio = 0; contatore = 0; fattore = 1.2;
testo = fprintf('Ciclo %d in corso',ii);
%%
while ii < MaxIt && ~uscita
tic
dJ_U(:) = 0;
dJ_V(:) = 0;
%% Definiamo Bu e Av        
    Av = ttv(A_j,exp(-V)-1,3);
    Av = (spmatrix(Av) + A_0) / DATA.rho_0;

    Bu = ttv(B_j,exp(-U)-1,3);
    Bu = (spmatrix(Bu) + B_0) / DATA.B_0;

for ff = 1:DATA.n_frq
    
    omega = DATA.omega(ff);
    omega2 = omega^2;
    k_0 = DATA.k_amb(ff);
    % Definiamo Ay e Fy              
    Atot = Av - Bu * omega2 + (1i*k_0 + 1/2/DATA.raggio_ext)/DATA.rho_0 * C_e;

    Fy =  D{ff} * ( (exp(-U) -1) /DATA.B_0 * omega2 ) ...   
         - E{ff} * ((exp(-V) -1) /DATA.rho_0);
    Fy = Fy - Atot(:,gamma_in) * Y(gamma_in,ff);
    % Impongo Dirichlet (omogeneo su p_tot) e risolvo in Y
    Y(~gamma_in,ff) = Atot(~gamma_in,~gamma_in) \ Fy(~gamma_in);

    
    % Definiamo Ap (= Atot) e Fp    
    Fp = B_a * Y(:,ff);
    % Impongo Dirichlet omogeneo e risolvo in P
    P(~gamma_in,ff) = conj(Atot(~gamma_in,~gamma_in))  \  Fp(~gamma_in);
 
    %% Calcolo del costo e dei gradienti per la freq ff  
    J(ii+1) = J(ii+1) + Y(:,ff)' * B_a * Y(:,ff);

    dJ_loc = real( P(:,ff)' *   ( E{ff} + spmatrix(ttv(A_j, Y(:,ff),2)) ));
    dJ_loc = dJ_loc(:) ./ exp(V) /DATA.rho_0 ;
    dJ_V = dJ_V + dJ_loc;

    dJ_loc = real( P(:,ff)' *   ( D{ff} + spmatrix(ttv(B_j, Y(:,ff), 2) ) ));
    dJ_loc = -omega2/DATA.B_0 * dJ_loc(:) ./ exp(U);
    dJ_U = dJ_U + dJ_loc;

end

% Costo sul gradiente del controllo
J_costo(ii+1) = + 2*lmb_u* U' * (diag(N_costo)-A_costo) * U + 2*lmb_v* V' * (N_costo-A_costo) * V;
J(ii+1) = J(ii+1) + J_costo(ii+1);
dJ_V = dJ_V + 4*lmb_v*(diag(N_costo) - A_costo)*V;
dJ_U = dJ_U + 4*lmb_u*(diag(N_costo) - A_costo)*U;
dJ_norma(ii+1,:) = [norm(dJ_V), norm(dJ_U)];

if ii == 1
    dJ0 = norm(dJ_U);
end

%% Stopping criteria    
if real(J(ii+1)) > real(J(ii))
    if tau > 1e-4
        tau = tau/1.4;
        U = U_old;
        V = V_old;
        J(ii+1) = 0;
        ritento = true;
        contatore = 0;
        fattore = fattore^0.9;
    else
        fprintf(repmat('\b',1,9));        fprintf('                Concluso in %.3f s.\n',toc);
        fprintf('Ottimizzazione terminata dopo %i step: \n',ii)
        fprintf('    anche con step piccoli, il costo non diminuisce,   tau = %5.4f \n',tau)
        uscita = 1;
    end
elseif norm(dJ_U) > tol*dJ0 - eps
    U_old = U;
    V_old = V;
    U = U - tau * dJ_U;
    V = V - tau * dJ_V;
%     [U,V] = proiettato(U,V,rho_lim,bulk_lim);
    
    ritento = false;
    ii = ii+1;
    if contatore == 10
        tau = tau * fattore;
        contatore = 0;
    else
        contatore = contatore + 1;
    end
        
else
    fprintf(repmat('\b',1,9));        fprintf('                Concluso in %.3f s.\n',toc);
    fprintf('\nOttimizzazione terminata con successo dopo %i step. \n',ii)
    fprintf('    Proseguendo, il costo diminuirebbe troppo poco. \n            |dJ|/|dJ0| = %3.2f*1e-4 < tol\n',1e4*norm(dJ_U)/dJ0)
    uscita = 1;
end

%% Disegnamo          
if mod(ii-2,disegna_ogni) == 0  || uscita || ii == MaxIt || ii <= 10 || ii == 50
    for ff = 1:DATA.n_frq
        subplot(2,3,ff)
        disegna_risultato(MESH,real(Y(1:MESH.numVertices,ff)) , true)
        title(['f = ' num2str(DATA.omega(ff)/2/pi,2) ' Hz,' ...
                    ' Direzione: (' num2str(cosdir(DATA.direz(ff))',2) ')'])
        subplot(2,3,ff+1)
        disegna_risultato(MESH,real(Y(1:MESH.numVertices,ff)+Y_inc(1:MESH.numVertices,ff)) , true)
        title(['f = ' num2str(DATA.omega(ff)/2/pi,2) ' Hz,' ...
                    ' Direzione: (' num2str(cosdir(DATA.direz(ff))',2) ')'])
    end

    subplot(235)
    cla
    disegna_u(exp(U),DATA)
    title('\kappa')
    xlabel('x'); ylabel('y'); view(2)
    colorbar

    subplot(236)
    cla
    disegna_u(exp(V),DATA)
    title('\rho')
    xlabel('x'); ylabel('y'); view(2)
    colorbar
    
    colormap('jet')
    
    subplot(234)
    cla
    plot(exp(V),exp(U),'x')
    xlabel('\rho');ylabel('\kappa');
    drawnow
    
end
    subplot(234)
    cla
    plot(exp(V),exp(U),'x')
%     hold on; plot(rho_lim,bulk_lim,'r'); hold off
    xlabel('\rho'); ylabel('\kappa');
%     xlim([0.5 1.5]); ylim([0.5 1.5])
    drawnow
%% Informazioni ciclo
tempo_medio = (tempo_medio*(ii-2) + toc)/(ii-1);
if ~uscita
    fprintf(repmat('\b',1,testo))
    testo = fprintf('Ciclo %d                Concluso in %.3f s.\n Media: %.2f  J = %.2f  tau = %.3f  T_tot: %.1f s\n',ii,toc,tempo_medio,J(ii),tau,tempo_medio*ii);
    if ii < MaxIt
        testo = testo + fprintf('Ciclo %d in corso',ii+1);
    end
end
end % Fine del while

if ii == MaxIt
    fprintf('Ottimizzazione terminata dopo %i step\n',ii)
    fprintf('    Il numero massimo di iterazioni è stato raggiunto\n')
end
% Disp costo
J = J(3:min(ii+1,MaxIt));
J_costo = J_costo(3:min(ii+1,MaxIt));
dJ_norma = dJ_norma(3:min(ii+1,MaxIt),:);
if all(abs(imag(J)) < abs(real(J)) * 1e-6)
    J = real(J);
    J_costo = real(J_costo);
else
    warning('J è complesso')
end
figure
plot(J); hold on
plot(dJ_norma(:,1)); plot(dJ_norma(:,2))
legend('J','dJ_V','dJ_U')
title('Costo VS iterazioni')
set(gca,'YScale','log')
%% Salvataggi
tempo_tot = tempo_medio*ii;
if tau ~=0 
    U = U_old;
    V = V_old;
    if lmb_u
        save(['Backup ' char(date) ' - soft - si vincoli'],'Y','U','V','J','MESH','DATA','tempo_tot','FE_SPACE','B_0','J_costo','dJ_norma')
    else
        save(['Backup ' char(date) ' - soft - no vincoli'],'Y','U','V','J','MESH','DATA','tempo_tot','FE_SPACE','B_0','J_costo','dJ_norma')
    end
end
beep
%% Plot poco utile
figure(4)
subplot(234)
plot(exp(V),exp(U),'x')
hold on
plot(rho_lim,bulk_lim,'r')
xlim([0 3]); ylim([0 3])
grid on
