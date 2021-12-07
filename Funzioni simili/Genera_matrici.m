function [A_0, A_j, A_j_conc, B_0, B_j,B_j_conc, C,D,E,F, M_D_a] = Genera_matrici(DATA, MESH, FE_SPACE)
%ADR_ASSEMBLER assembler for 2D ADR equations
% https://it.mathworks.com/help/matlab/cc-mx-matrix-library.html?s_tid=CRUX_topnav

%% Assemblo A_0 e B_0
[Arows, Acols, Acoef, Mcoef] = Assembla_A_e_M(MESH.dim, MESH.elements, FE_SPACE.numElemDof, FE_SPACE.quad_weights, ...
                                              MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);

A_0    = sparse(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);     % Matrice di rigidezza
B_0    = sparse(Arows,Acols,Mcoef,MESH.numNodes,MESH.numNodes);     % Matrice di massa

%% Creo i nodi di quadratura
coord_ref = MESH.chi;
x = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
y = x;
for j = 1 : 3
    i = MESH.elements(j,:);
    vtemp = MESH.vertices(1,i);
    x = x + vtemp'*coord_ref(j,:);
    vtemp = MESH.vertices(2,i);
    y = y + vtemp'*coord_ref(j,:);
end

%% Assemblo A_j
% Evaluation of coefficients in the quadrature nodes
A_j{DATA.n_basi,1} = [];
A_j_conc = [];

% mu  = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
% f = mu;
% b = [mu mu];
% si = mu;

for nn = 1:DATA.n_basi

    mu  = DATA.psi(x,y,nn,MESH.vertices);
    [Arows, Acols, Acoef] = Assembla_solo_mu(MESH.dim, MESH.elements, FE_SPACE.numElemDof, mu, FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.dphi_ref);
%     mu2 = DATA.psi(x,y,2,[]);
%     [Arows, Acols, Acoef, ~, ~, ~] = ADR_assembler_C_omp(MESH.dim, 'all', [10 10], 10, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,  ...
%         FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
    [i,j,v] = find(sparse(Arows,Acols,Acoef)); %,  MESH.numNodes,MESH.numNodes));
    A_j{nn} = [i,j,v];
   
    A_j_conc = [A_j_conc; sparse(Arows,Acols,Acoef,  MESH.numNodes,MESH.numNodes)];
end

%% Assemblo B_j
B_j{DATA.n_basi,1} = [];
B_j_conc = [];
% mu(:) = 0;
for nn = 1:DATA.n_basi

    si  = DATA.psi(x,y,nn,MESH.nodes);
    [Brows, Bcols, Bcoef] = Assembla_solo_si(MESH.dim, MESH.elements, FE_SPACE.numElemDof, si, FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi);
    
%     [Brows, Bcols, Bcoef, ~, ~, ~] = ADR_assembler_C_omp(MESH.dim, 'all', [10 10], 10, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,  ...
%         FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
    [i,j,v] = find(sparse(Brows,Bcols,Bcoef)); %,  MESH.numNodes,MESH.numNodes));
    B_j{nn} = [i,j,v];
    
    B_j_conc = [B_j_conc; sparse(Brows,Bcols,Bcoef,  MESH.numNodes,MESH.numNodes)];
end

%% Assemblo la forzante D

D{DATA.n_frq} = [];
D_loc{DATA.n_basi} = [];

for ff = 1:DATA.n_frq
    f_ff = DATA.y_inc_fun(x,y,DATA.k_amb(ff),DATA.direz(:,ff));
    switch DATA.tipo_base
        case {'rettangoli', 'corone circolari','hex'}
            for nn = 1:DATA.n_basi
                % complex() serve per evitare il caso in cui f sia reale -> Assembla_D sfolla
                f =  complex(f_ff .* DATA.psi(x,y,nn,MESH.vertices) );
                [Frows, Fcoef] = Assembla_D(MESH.dim, MESH.elements, FE_SPACE.numElemDof,f,FE_SPACE.quad_weights,MESH.jac,FE_SPACE.phi);

                D_loc{nn} = sparse(Frows,1,Fcoef, MESH.numNodes,1);
            end
            D{ff} = cat(2,D_loc{:});
        otherwise
            error('DATA.tipo_base non riconosciuto')
    end
end

%% Assemblo la forzante E

E{DATA.n_frq} = [];
E_loc{DATA.n_basi} = [];

for ff = 1:DATA.n_frq
    f_ff_x = DATA.dy_inc_dx_fun(x,y,DATA.k_amb(ff),DATA.direz(:,ff));
    f_ff_y = DATA.dy_inc_dy_fun(x,y,DATA.k_amb(ff),DATA.direz(:,ff));
    switch DATA.tipo_base
        case {'rettangoli', 'corone circolari','hex'}
            for nn = 1:DATA.n_basi
                % complex() serve per evitare il caso in cui f sia reale -> Assembla_E sfolla
                valore_base = DATA.psi(x,y,nn,MESH.vertices);
                f_x = complex(f_ff_x .* valore_base);
                f_y = complex(f_ff_y .* valore_base );
                [Frows, Fcoef] = Assembla_E(MESH.dim, MESH.elements, FE_SPACE.numElemDof, f_x, f_y, ...
                                            FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);

                E_loc{nn} = sparse(Frows,1,Fcoef, MESH.numNodes,1);
            end
            E{ff} = cat(2,E_loc{:});
        otherwise
            error('DATA.tipo_base non riconosciuto')
    end
end

%% Assemblo la matrice di massa del dominio ambiente
mu  = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
f_real = mu;
b = [mu mu];
 
[idx_bordo_mant, idx_bordo_osta] = Calcola_cont_mant_e_ost(MESH.boundaries,MESH.nodi_boundaries,DATA.tipo_contorno,MESH.normali,DATA.tipo_dominio);
si = ~isdentro(x,y, idx_bordo_mant, idx_bordo_osta,MESH.vertices) *1;


[Arows, Acols, Acoef, ~, ~, ~] = ADR_assembler_C_omp(MESH.dim, 'all', [10 10], 10, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f_real,  ...
                                                 FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);

M_D_a = sparse(Arows,Acols,Acoef,  MESH.numNodes,MESH.numNodes);

%% Definiamo la lista di nodi del contorno esterno e interno
gamma_ext = find( MESH.nodi_boundaries(:,1) & ~MESH.nodi_boundaries(:,2)); % Appartengono al dominio 1 (=fluido), ma non al 2 (=mantello)

gamma_in =  MESH.nodi_boundaries(:,2) & ~MESH.nodi_boundaries(:,1); % Appartengono al dominio 2, ma non all' 1

lato_sopra = abs(MESH.vertices(2,MESH.boundaries(1,:)))==0.5 & abs(MESH.vertices(2,MESH.boundaries(2,:)))==0.5;
gamma_in = gamma_in | lato_sopra(:);
gamma_ext = setdiff(gamma_ext, find(lato_sopra(:)));

%% Assemblo C
% Matrice di massa sul contorno esterno           
[csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
[phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);

nof         = length(gamma_ext);
nbn         = MESH.numBoundaryDof;

Crows       = zeros(nbn*nbn*nof,1);
Ccols       = Crows;
Ccoef       = Crows;

[rows,cols] = meshgrid(1:nbn,1:nbn);
rows        = rows(:);
cols        = cols(:);


x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma_ext));
y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma_ext));

side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);

coso = ( repmat(wi,nbn,1) .* phi )* phi';
for ll = 1 : nof
    face = gamma_ext(ll);

    Crows(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(rows,face);
    Ccols(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(cols,face);
    Ccoef(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  side_length(ll) * coso;

end

C = sparse(Crows,Ccols,Ccoef,MESH.numNodes,MESH.numNodes);

%% Creo le normali ai lati della mesh

x_1    =  MESH.vertices(1,MESH.boundaries(1, gamma_in));
x_2    =  MESH.vertices(1,MESH.boundaries(2, gamma_in));
y_1    =  MESH.vertices(2,MESH.boundaries(1, gamma_in));
y_2    =  MESH.vertices(2,MESH.boundaries(2, gamma_in));

ha_a_dx_lo_0 = MESH.boundaries(7,gamma_in) == 0;
Dx = x_2 - x_1;
Dy = y_2 - y_1;
direzione = [Dx; Dy];

direzione(:,ha_a_dx_lo_0) =  [0 1; -1 0] * direzione(:,ha_a_dx_lo_0);
direzione(:,~ha_a_dx_lo_0) = [0 -1; 1 0] * direzione(:,~ha_a_dx_lo_0);

direzione = direzione ./ repmat(vecnorm(direzione),2,1);

%% Assemblo F
% forzante su gamma_in

[csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
[phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
eta            =  1 - csi;
nqn            =  length(csi);


nof         = sum(gamma_in);
nbn         = MESH.numBoundaryDof;

Rrows       = zeros(nbn*nof,1);
Rcoef       = Rrows;


xlt = zeros(nof,nqn); ylt = xlt;
coord_ref = [eta; csi];
for j = 1 : 2
    dof = MESH.boundaries(j,gamma_in);
    vtemp = MESH.vertices(1,dof);
    xlt = xlt + vtemp'*coord_ref(j,:);
    vtemp = MESH.vertices(2,dof);
    ylt = ylt + vtemp'*coord_ref(j,:);

end

x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma_in));
y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma_in));

side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);

F{DATA.n_frq} = [];
gamma_in = find(gamma_in);
for ff = 1:DATA.n_frq
    u_Robin_x = DATA.dy_inc_dx_fun(xlt,ylt,DATA.k_amb(ff),DATA.direz(:,ff));
    u_Robin_y = DATA.dy_inc_dy_fun(xlt,ylt,DATA.k_amb(ff),DATA.direz(:,ff));

    u_Robin = u_Robin_x .* repmat(direzione(1,:)',1,nqn) + u_Robin_y .* repmat(direzione(2,:)',1,nqn);

    for ll = 1 : nof
        face = gamma_in(ll);

        u_Robin_loc  = u_Robin(ll,:).*wi;
%         u_Robin_loc  = u_Robin_loc(1,:).';

        Rrows(1+(ll-1)*nbn:ll*nbn)    = MESH.boundaries(1:nbn,face);
        Rcoef(1+(ll-1)*nbn:ll*nbn)    = side_length(ll)*phi*(u_Robin_loc).';

    end

    F{ff} = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
end
    