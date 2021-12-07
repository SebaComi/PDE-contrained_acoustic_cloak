function [dP_n, dP_x_mesh, dP_y_mesh] = derivata_contorno(P, boundaries, DATA_dP, MESH_dP, FE_SPACE_dP )
% Questa funzione crea le derivate direzionali lungo la normale al dominio
% in tutti i punti del contorno

DATA_dP.force = @(x,y,t,param) 0.*x;
DATA_dP.diffusion = @(x,y,t,param) 0.*x;

% Inizializzo dP_n
dP_n = zeros(size(boundaries,2),1);
% [1 2 3 4] sono rispettivamente le rette [y = 0, x = 1, y = 1, x = 0] a
% cui corrispondono le seguenti direzioni normali: [0 -1],  [1 0],  [0 1],  [-1 0]

%% Calcolo dP_n lungo la direzione y := [0 1] per i lati y = 1 e y = 0
DATA_dP.transport{1} = @(x,y,t,param) 0 + 0.*x;
DATA_dP.transport{2} = @(x,y,t,param) 1 + 0.*x;

[A_dP, ~, M_dP] = ADR_Assembler(MESH_dP,DATA_dP,FE_SPACE_dP);
dP_y_mesh = M_dP\A_dP*P;     % Derivata direzionale su tutti i punti della mesh
posizione = boundaries(5,:) == 3;       % Mi fornisce la posizione su boundaries
dP_n(posizione) =  dP_y_mesh(boundaries(1,posizione));
posizione = boundaries(5,:) == 1;
dP_n(posizione) = -dP_y_mesh(boundaries(1,posizione));

%% Calcolo dP_n lungo la direzione x := [1 0] per i lati x = 1 e x = 0
DATA_dP.transport{1} = @(x,y,t,param) 1 + 0.*x;
DATA_dP.transport{2} = @(x,y,t,param) 0 + 0.*x;

[A_dP, ~, M_dP] = ADR_Assembler(MESH_dP,DATA_dP,FE_SPACE_dP);
dP_x_mesh = M_dP\A_dP*P;     % Derivata direzionale su tutti i punti della mesh
posizione = boundaries(5,:) == 2;
dP_n(posizione) =  dP_x_mesh(boundaries(1,posizione));
posizione = boundaries(5,:) == 4;
dP_n(posizione) = -dP_x_mesh(boundaries(1,posizione));

%% Ci manca la dP_n nei 4 vertici, perciò componiamo i due risultati che abbiamo già

vertice_00 = MESH_dP.vertices(1,boundaries(1,:)) == 0 & MESH_dP.vertices(2,boundaries(1,:)) == 0;   % Prendo la posizione su boundaries dei vertici
vertice_10 = MESH_dP.vertices(1,boundaries(1,:)) == 1 & MESH_dP.vertices(2,boundaries(1,:)) == 0;
vertice_01 = MESH_dP.vertices(1,boundaries(1,:)) == 0 & MESH_dP.vertices(2,boundaries(1,:)) == 1;
vertice_11 = MESH_dP.vertices(1,boundaries(1,:)) == 1 & MESH_dP.vertices(2,boundaries(1,:)) == 1;

posizione_00 = boundaries(1,vertice_00);       % Prendo la posizione sulla mesh dei vertici
posizione_10 = boundaries(1,vertice_10);
posizione_01 = boundaries(1,vertice_01);
posizione_11 = boundaries(1,vertice_11);

dP_n(vertice_00) = sqrt(2) * (- dP_x_mesh(posizione_00) - dP_y_mesh(posizione_00) );
dP_n(vertice_10) = sqrt(2) * (+ dP_x_mesh(posizione_10) - dP_y_mesh(posizione_10) );
dP_n(vertice_01) = sqrt(2) * (- dP_x_mesh(posizione_01) + dP_y_mesh(posizione_01) );
dP_n(vertice_11) = sqrt(2) * (+ dP_x_mesh(posizione_11) + dP_y_mesh(posizione_11) );
end