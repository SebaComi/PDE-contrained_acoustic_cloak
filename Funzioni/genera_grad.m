function [dJ_U, dJ_V, J] = genera_grad(U,V,DATA,MESH)

global A_0 A_j A_j_conc B_0 B_j B_j_conc  C D E F  M_D_a lunghezze lmb_u lmb_v

% Y_inc nei nodi
Y_inc = zeros(MESH.n_vertices,DATA.n_frq);
for ff = 1:DATA.n_frq
    Y_inc(:,ff) = DATA.y_inc_fun(MESH.vertices(1,:),MESH.vertices(2,:),DATA.k_amb(ff),DATA.direz(:,ff)).';
end

Y = zeros(MESH.n_vertices,DATA.n_frq);
P = Y;

dJ_U = zeros(DATA.n_basi,1);
dJ_V = dJ_U;
J = 0;
%% Definiamo Bu e Av        
array_A = cat(1,A_j{:});
array_B = cat(1,B_j{:});

    for jj = 1:DATA.n_basi        
        array_A(lunghezze(jj)+1:lunghezze(jj+1),3) = (exp(-V(jj))-1) * array_A(lunghezze(jj)+1:lunghezze(jj+1),3);        % D1_rho = (1/rho - 1/rho_0)
        array_B(lunghezze(jj)+1:lunghezze(jj+1),3) = (exp(-U(jj))-1) * array_B(lunghezze(jj)+1:lunghezze(jj+1),3);     % DB = (1/B - 1/B_0)
    end
    Av = sparse(array_A(:,1),array_A(:,2),array_A(:,3),MESH.n_vertices,MESH.n_vertices);
    Av = (Av + A_0)/DATA.rho_0;  % Cioè (1/rho - 1/rho_0) + 1/rho_0 = 1/rho
    Bu = sparse(array_B(:,1),array_B(:,2),array_B(:,3),MESH.n_vertices,MESH.n_vertices);
    Bu = (Bu + B_0) / DATA.B_0;  % Cioè - delta B  - B_0 = - B
%%
for ff = 1:DATA.n_frq
    
%     k_amb = DATA.k_amb(ff);
    omega = DATA.omega(ff);
    omega2 = omega^2;
    k_0 = DATA.k_amb(ff);
    %% Definiamo Ay e Fy              
    Atot = Av - Bu * omega2 + (1i*k_0 + 1/DATA.lato_q)/DATA.rho_0 * C;
    Fy = D{ff} * ( (exp(-U) -1) /DATA.B_0 * omega2 ) ...   
         - E{ff} * ((exp(-V) -1) /DATA.rho_0) ...
         - F{ff} / DATA.rho_0;
     
    %% Solve for Y (Questo solving forse si potrebbe fare meglio)
    Y(:,ff) = Atot \ Fy;
    
    %% Definiamo Ap (= Atot) e Fp    
    Fp = 2* M_D_a * conj(Y(:,ff));   % derivo y*T M y ---> 2y*T M traspongo 2M y*

    %% Solve for P          
    P(:,ff) = Atot  \  Fp;
        
    %% Calcolo del costo e dei gradienti per la freq ff  
    J = J + Y(:,ff)' * M_D_a * Y(:,ff);
    
    dJ_loc = real( sparse(transpose(P(:,ff))) * (E{ff} +  ...
                   reshape(A_j_conc * sparse(Y(:,ff)) ,MESH.n_vertices,DATA.n_basi) ...
                  ));                                                                                 % non dovrebbe essere trasposto? Si, infatti alla riga dopo lo uso come colonna
    dJ_loc = dJ_loc(:) .* exp(-V) /DATA.rho_0  ;
    dJ_V = dJ_V + dJ_loc;
    
    dJ_loc = real( sparse(transpose(P(:,ff))) * ( D{ff} +  ...
                   reshape(B_j_conc * sparse(Y(:,ff)) ,MESH.n_vertices,DATA.n_basi) ...
                  ));
    dJ_loc = -omega2/DATA.B_0 * dJ_loc(:) ./ exp(U);
    dJ_U = dJ_U + dJ_loc;

end
%% Aggiungo il costo e i gradienti dovuti a U e V
J = J + lmb_u/2 * (U)'*(DATA.area .*U) + lmb_v/2 * (V)'*(DATA.area .*V);    % Perché la matrice di massa del controllo è una matrice diagonale
dJ_U = dJ_U + lmb_u * (DATA.area .*U);
dJ_V = dJ_V + lmb_v * (DATA.area .*V);

end
