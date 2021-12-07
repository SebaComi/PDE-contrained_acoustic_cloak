function [A, B, M, C_x, C_y, D, B_conc, A_rho, A_rho_conc] = MIO_ADR_Assembler(MESH, FE_SPACE, nodi_dominio,DATA,richiesta,flag_B)
%ADR_ASSEMBLER assembler for 2D/3D ADR equations
% https://it.mathworks.com/help/matlab/cc-mx-matrix-library.html?s_tid=CRUX_topnav
% Costruiamo le matrici di massa, rigidezza e la cella di matrici B
% A matrice di rigidezza
% B è un array di celle contenenti sparse double
% M è la matrice di massa

if ~exist('flag_B','var')
    flag_B = 'lista';
end

%% Assembly A and M
if any(richiesta == 1) || any(richiesta == 3)

    [Arows, Acols, Acoef, Mcoef] = Assembla_A_e_M(MESH.dim, MESH.elements, FE_SPACE.numElemDof, FE_SPACE.quad_weights, ...
                                                  MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
    A    = sparse(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);
    M    = sparse(Arows,Acols,Mcoef,MESH.numNodes,MESH.numNodes);
else
    A = [];
    M = [];
end
%% Assembly A_var in caso di controllo su rho
if any(richiesta == 8)
    switch DATA.tipo_base
        case {'rettangoli', 'corone circolari','hex'}
            coord_ref = MESH.chi;
            x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x;
            for j = 1 : 3
                i = MESH.elements(j,:);
                vtemp = MESH.vertices(1,i);
                x = x + vtemp'*coord_ref(j,:);
                vtemp = MESH.vertices(2,i);
                y = y + vtemp'*coord_ref(j,:);
            end
            % Evaluation of coefficients in the quadrature nodes
            si  = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
            f_real = si;
            b = [si si];
            A_rho{DATA.n_basi,1} = [];
            A_rho_conc = [];
            for nn = 1:DATA.n_basi

                mu  = DATA.psi(x,y,nn,MESH.vertices);
                [Arows, Acols, Acoef, ~, ~, ~] = ADR_assembler_C_omp(MESH.dim, 'all', [10 10], 10, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f_real,  ...
                                                                 FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
                [i,j,v] = find(sparse(Arows,Acols,Acoef)); %,  MESH.numNodes,MESH.numNodes));
                A_rho{nn} = [i,j,v];
                if any(richiesta == 9)
                    A_rho_conc = [A_rho_conc; sparse(Arows,Acols,Acoef,  MESH.numNodes,MESH.numNodes)];
                end
            end
        otherwise
            error('A_rho non disponibile')
    end
else
    A_rho = [];
end
if ~exist('A_rho_conc','var')
    A_rho_conc = [];
end

%% Assembly B, a 3D sparse matrix, as a rectangular matrix
if any(richiesta == 2)
switch DATA.tipo_base
    case 'vertici'
        [Brows, Bcols, Balte, Bcoef] = Assembla_B( MESH.dim, MESH.elements, FE_SPACE.numElemDof, ...
                                                   FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi );
        switch flag_B
            case 'cell'
                B{MESH.numNodes} = [];
                for ii = 1:MESH.numNodes
                    dove = Balte == ii;
                    B{ii} = sparse(Brows(dove),Bcols(dove),Bcoef(dove),  MESH.numNodes,MESH.numNodes);
                end
            case 'conc'
                B = [];
                for ii = 1:MESH.numNodes
                    dove = Balte == ii;
                    B = [B; sparse(Brows(dove),Bcols(dove),Bcoef(dove),  MESH.numNodes,MESH.numNodes)];
                end
            case 'lista'
                B{MESH.numNodes,1} = [];
                for ii = find(nodi_dominio(:)') % 1:MESH.numNodes % con il find ci risparmiamo un po' di iterazioni inutili
                    dove = Balte == ii;
                    [i,j,v] = find(sparse(Brows(dove),Bcols(dove),Bcoef(dove)));%,  MESH.numNodes,MESH.numNodes));
                    B{ii} = [i,j,v];
                end
            otherwise
                error(['DATA.tipo_base = ''' DATA.tipo_base 'non riconosciuto'])
        end
        D = [];
    case {'rettangoli', 'corone circolari','hex'}
        coord_ref = MESH.chi;
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x;
        for j = 1 : 3
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
        end
        % Evaluation of coefficients in the quadrature nodes
        mu  = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
        f_real = mu;
        b = [mu mu];
        B{DATA.n_basi,1} = [];
        B_conc = [];
        for nn = 1:DATA.n_basi

            si  = DATA.psi(x,y,nn,MESH.vertices);
            [Brows, Bcols, Bcoef, ~, ~, ~] = ADR_assembler_C_omp(MESH.dim, 'all', [10 10], 10, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f_real,  ...
                                                             FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
            [i,j,v] = find(sparse(Brows,Bcols,Bcoef)); %,  MESH.numNodes,MESH.numNodes));
            B{nn} = [i,j,v];
            if any(richiesta == 7)
                B_conc = [B_conc; sparse(Brows,Bcols,Bcoef,  MESH.numNodes,MESH.numNodes)];
            end
        end
otherwise
    error('tipo_base non riconosciuto')
end
else
    B = [];
end
if ~exist('B_conc','var')
    B_conc = [];
end


%% Assembly of transport matyrices C
% C.x e C.y sono le matrici dove il trasporto è (1,0) e (0,1) rispettivamente
if any(richiesta == 4)
    [Crows, Ccols, Ccoef_x, Ccoef_y] = Assembla_C(MESH.dim, MESH.elements, FE_SPACE.numElemDof, FE_SPACE.quad_weights, ...
                                                    MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);

    C_x     = sparse(Crows,Ccols,Ccoef_x,MESH.numNodes,MESH.numNodes);
    C_y     = sparse(Crows,Ccols,Ccoef_y,MESH.numNodes,MESH.numNodes);
else
    C_x = [];
    C_y = [];
end

%% Assemblo la forzante D
if any(richiesta == 6)
    D{DATA.n_frq} = [];
            coord_ref = MESH.chi;
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x;

        for j = 1 : 3
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
        end
        D_loc{DATA.n_basi} = [];
    for ff = 1:DATA.n_frq
        switch DATA.tipo_base
            case {'rettangoli', 'corone circolari','hex'}
                for nn = 1:DATA.n_basi
                    % complex() serve per evitare il caso in cui f sia reale -> Assembla_D sfolla
                    f = complex(DATA.y_inc_fun(x,y,DATA.k_amb(ff),DATA.direz(:,ff)) .* DATA.psi(x,y,nn,MESH.vertices) );
                    [Frows, Fcoef] = Assembla_D(MESH.dim, MESH.elements, FE_SPACE.numElemDof,f,FE_SPACE.quad_weights,MESH.jac,FE_SPACE.phi);
                    
                    D_loc{nn} = sparse(Frows,1,Fcoef, MESH.n_vertices,1);
                end
                D{ff} = cat(2,D_loc{:});
        end
    end
else
    D = [];
end
    
    
    