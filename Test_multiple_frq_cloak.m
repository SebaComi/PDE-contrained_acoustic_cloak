
load('Circ soft - no vincoli - Direz ortog.mat')    

DATA = Dati_cloak(DATA);
DATA.psi = @(xy,nn,vertices) esagono(xy,nn,DATA);
DATA.n_frq = 51;
% DATA.omega = 2000*pi*linspace(10,16,DATA.n_frq);
DATA.omega = 2*pi*13000*ones(DATA.n_frq,1);
DATA.k_amb = DATA.omega / DATA.c_0;
DATA.direz = linspace(0, pi/2, DATA.n_frq);
MESH.bound_tipo_cont = 5;
quali = {'A_0','A_j', 'B_0','B_j', 'C_e','D  ','E  ','B_a'};
        [ A_0,  A_j,   B_0,  B_j,   C_e,  D,    E,    B_a] = Genera_matrici_parfor(DATA, MESH, FE_SPACE,quali);

Av = ttv(A_j,exp(-V)-1,3);
Av = (spmatrix(Av) + A_0) / DATA.rho_0;
Bu = ttv(B_j,exp(-U)-1,3);
Bu = (spmatrix(Bu) + B_0) / DATA.B_0;

% Dirichlet
gamma_in = trova_gamma_in(DATA,MESH,'s');


for ff = 1:DATA.n_frq
    Y(gamma_in,ff) = -DATA.y_inc_fun({MESH.nodes(1,gamma_in),MESH.nodes(2,gamma_in)},DATA.k_amb(ff),DATA.direz(ff));

    omega = DATA.omega(ff);
    omega2 = omega^2;
    k_0 = DATA.k_amb(ff);
    % Definiamo Ay e Fy              
    Atot = Av - Bu * omega2 + (1i*k_0 + 1/2/DATA.raggio_ext)/DATA.rho_0 * C_e;

    Fy =  D{ff} * ( (exp(-U) -1) /DATA.B_0 * omega2 )  - E{ff} * ((exp(-V) -1) /DATA.rho_0);
    Fy = Fy - Atot(:,gamma_in) * Y(gamma_in,ff);
    % Impongo Dirichlet (omogeneo su p_tot) e risolvo in Y
    Y(~gamma_in,ff) = Atot(~gamma_in,~gamma_in) \ Fy(~gamma_in);
    
    if ff < 50
        subplot(7,7,ff)
        disegna_risultato(MESH,real(Y(:,ff)) , true)
    end
    J_ff(ff) = real(Y(:,ff)' * B_a * Y(:,ff));
    
end

return
%%
theta = [flip(-DATA.direz), DATA.direz];
theta = [theta, pi-flip(theta)];
rrr = [flip(J_ff), J_ff];
rrr = [rrr, rrr];
figure
polarplot(theta,rrr )
return
%%
figure
plot(DATA.omega/2/pi/1000,J_ff)

%%
% save('Circ soft - si vincoli - 02-Nov-2021- fine - 1 direz - direzioni.mat','Y','J_ff')
save(' Circ soft - no vincoli - Direz ortog - direzioni.mat','Y','J_ff')
