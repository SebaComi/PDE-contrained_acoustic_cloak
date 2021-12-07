function [out1, out2,  out3, out4, out5, out6, out7, out8, out9, out10, out11, out12, out13, out14, out15, out16, out17] = Genera_matrici_parfor(DATA, MESH, FE_SPACE,quali)
%ADR_ASSEMBLER assembler for 2D ADR equations
% https://it.mathworks.com/help/matlab/cc-mx-matrix-library.html?s_tid=CRUX_topnav
fprintf('Generazione delle matrici... '); tic
quali_vec = cat(1,quali{:});
%% Assemblo A_0 e B_0      
% matrici di rigidezza e massa dell'intera mesh
if any(strcmp(quali,'A_0')) || any(strcmp(quali,'B_0'))
    [Rows, Cols, CoefA, CoefB] = Assembla_A_e_M(MESH.dim, MESH.elements, FE_SPACE.numElemDof, FE_SPACE.quad_weights, ...
                                                  MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
    
    A_0    = sparse(Rows,Cols,CoefA,MESH.numNodes,MESH.numNodes);     % Matrice di rigidezza
    B_0    = sparse(Rows,Cols,CoefB,MESH.numNodes,MESH.numNodes);     % Matrice di massa
end
%% Creo i nodi di quadratura
coord_ref = MESH.chi;
switch MESH.dim
    case 2
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes);      y = x;
        for j = 1 : 3
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
        end
        nodi_quad = {x,y};
    case 3
            x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x; z = x;
        for j = 1:4
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,i);
            z = z + vtemp'*coord_ref(j,:);
        end
        nodi_quad = {x,y,z};
end
%% Assemblo A_j             
if any(strcmp(quali,'A_j'))
    A_j{DATA.n_basi} = [];      k = A_j;
    parfor nn = 1:DATA.n_basi
        mu  = DATA.psi(nodi_quad,nn,MESH.vertices);
        [Rows, Cols, Coef] = Assembla_solo_mu(MESH.dim, MESH.elements, FE_SPACE.numElemDof, mu, FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.dphi_ref);
        [i,j,v] = find(sparse(Rows,Cols,Coef)); %,  MESH.numNodes,MESH.numNodes));
        A_j{nn} = [i,j,v];
        k{nn} = nn * ones(length(i),1);
    end
    A_j = cat(1,A_j{:});     k = cat(1,k{:});
    A_j = sptensor([A_j(:,1:2), k], A_j(:,3), [[1 1]*MESH.numNodes   DATA.n_basi]);
end
%% Assemblo B_j             
if any(strcmp(quali,'B_j'))
    B_j{DATA.n_basi} = [];    k = B_j;
    parfor nn = 1:DATA.n_basi
        si  = DATA.psi(nodi_quad,nn,MESH.vertices);
        [Rows, Cols, Coef] = Assembla_solo_si(MESH.dim, MESH.elements, FE_SPACE.numElemDof, si, FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi);
        [i,j,v] = find(sparse(Rows,Cols,Coef)); %,  MESH.numNodes,MESH.numNodes));
        B_j{nn} = [i,j,v];
        k{nn} = nn * ones(length(i),1);
    end
    B_j = cat(1,B_j{:});     k = cat(1,k{:});
    B_j = sptensor([B_j(:,1:2), k], B_j(:,3), [[1 1]*MESH.numNodes   DATA.n_basi]);
end

%% Assemblo la matrice di massa del dominio ambiente B_a (vecchia M_a o M_D_a)  
% Matrice tale che [B_a]_{k,i} = \int_{D_a} \phi_k phi_i d\Omega
dove = find(quali_vec(:,1)' == 'B' & ismember(quali_vec(:,3)', 'abc'));
for jj = dove(:)'
    lettera = quali_vec(jj,3);
    chi_dominio = find(DATA.tipo_dominio == lettera,1);
    si = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
    si(MESH.elements(end,:) == chi_dominio, :) = 1;
    [Rows, Cols, Coef] = Assembla_solo_si(MESH.dim, MESH.elements, FE_SPACE.numElemDof, si, FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi);
    eval(['B_' lettera ' = sparse(Rows,Cols,Coef,  MESH.numNodes,MESH.numNodes);'])
end

%% Assemblo la forzante D       
% Vettore tale che [D]_{k,j} = \int_\Omega  \phi_k \psi_j  p_{inc} d\Omega
if any(strcmp(quali,'D  '))
    D{DATA.n_frq} = [];
    parfor ff = 1:DATA.n_frq
        D_loc = cell(DATA.n_basi,1);
        f_ff = DATA.y_inc_fun(nodi_quad,DATA.k_amb(ff),DATA.direz(ff));
        for nn = 1:DATA.n_basi
            % complex() serve per evitare il caso in cui f sia reale -> Assembla_D sfolla
            f =  complex(f_ff .* DATA.psi(nodi_quad,nn,MESH.vertices) );
            [Frows, Fcoef] = Assembla_D(MESH.dim, MESH.elements, FE_SPACE.numElemDof,f,FE_SPACE.quad_weights,MESH.jac,FE_SPACE.phi);

            D_loc{nn} = sparse(Frows,1,Fcoef, MESH.numNodes,1);
        end
        D{ff} = cat(2,D_loc{:});
    end
end
%% Assemblo la forzante E       
if any(strcmp(quali,'E  '))
    E{DATA.n_frq} = [];
    parfor ff = 1:DATA.n_frq
        E_loc = cell(DATA.n_basi,1);
        f_ff_x = DATA.dy_inc_dx_fun(nodi_quad,DATA.k_amb(ff),DATA.direz(ff));
        f_ff_y = DATA.dy_inc_dy_fun(nodi_quad,DATA.k_amb(ff),DATA.direz(ff));
                for nn = 1:DATA.n_basi
                    % complex() serve per evitare il caso in cui f sia reale -> Assembla_E sfolla
                    valore_base = DATA.psi(nodi_quad,nn,MESH.vertices);
                    f_x = complex(f_ff_x .* valore_base);
                    f_y = complex(f_ff_y .* valore_base );
                    [Frows, Fcoef] = Assembla_E(MESH.dim, MESH.elements, FE_SPACE.numElemDof, f_x, f_y, ...
                                                FE_SPACE.quad_weights, MESH.invjac, MESH.jac, FE_SPACE.phi, FE_SPACE.dphi_ref);
                    E_loc{nn} = sparse(Frows,1,Fcoef, MESH.numNodes,1);
                end
                E{ff} = cat(2,E_loc{:});
    end
end

%% Assemblo la forzante del secondo dominio ambiente F_b    -  Questa è D_b, quindi va aggiornata facendo D_x per ogni x
% Vettore tale che [F_b]_k = \int_{D_b} \phi_k p_{inc} d\Omega
if any(strcmp(quali,'F_b'))
    F_b{DATA.n_frq} = [];
    chi_dominio = find(DATA.tipo_dominio == 'b');
%     FaceID = nearestFace(MESH.gm,[x(:),y(:)]);
%     dove = reshape(FaceID == chi_dominio, size(x));
%     dove = reshape(trova_dominio([x(:),y(:)]', MESH, chi_dominio),size(x));
    dove = false(size(x));
    dove(MESH.elements(end,:) == chi_dominio, :) = true;
    for ff = 1:DATA.n_frq
        f_ff = complex( DATA.y_inc_fun(nodi_quad,DATA.k_amb(ff),DATA.direz(ff))  );     % complex() serve per evitare il caso in cui f sia reale -> Assembla_D sfolla
        f_ff(~dove) = 0;
        [Frows, Fcoef] = Assembla_D(MESH.dim, MESH.elements, FE_SPACE.numElemDof,f_ff,FE_SPACE.quad_weights,MESH.jac,FE_SPACE.phi);
        F_b{ff} = sparse(Frows,1,Fcoef, MESH.numNodes,1);
    end
    clearvars f_ff
end
%% Definiamo la lista di nodi del contorno esterno (e), rigido (r), di controllo (c) e altri
for lettera = intersect(quali_vec(:,3)', 'ercxzgf')
    eval(['gamma_' lettera '= [];'])
    for index = find(DATA.tipo_contorno == lettera)
        eval(['gamma_' lettera '= [gamma_' lettera ', find(MESH.boundaries(MESH.bound_tipo_cont,:) == index)];'])
    end
end

%% Assemblo C_x - Matrice di massa sui contorni x
% [C_x]_{i,k} = \int_{\Gamma_x} { \phi_i \phi_k  d\Gamma }
dove = find(quali_vec(:,1)' == 'C');
for jj = dove(:)'
    lettera = quali_vec(jj,3);
    if ~strcmp(lettera, 'j' )
        eval(['C_' lettera '= massa_contorno_C(DATA,MESH,FE_SPACE,gamma_' lettera ');'])
    end
end
%% Assemblo F_r (bordo rigido) e F_z (bordo impedente) e altre
% [F_x]_k = \int_{\Gamma_x} { \nabla p_inc \cdot n \phi_k d\Gamma }
dove = find(quali_vec(:,1)' == 'F' & ~ismember(quali_vec(:,3)', 'b'));
for jj = dove(:)'
    lettera = quali_vec(jj,3);
    eval(['F_' lettera '= forzante_F(DATA,MESH,FE_SPACE,gamma_' lettera ');'])
end
%% Assemblo G_r (bordo rigido) e G_z (bordo impedente) e altre
% [G_x]_k = \int_{\Gamma_x} { p_inc \phi_k d\Gamma } 
dove = find(quali_vec(:,1)' == 'G' & ~ismember(quali_vec(:,3)','b'));
for jj = dove(:)'
    lettera = quali_vec(jj,3);
    eval(['G_' lettera '= forzante_G(DATA,MESH,FE_SPACE,gamma_' lettera ');'])
end
%% Assemblo H_g (bordo galleggiante) e altre
% [H_x]_k = \sum_{dd = 1}^{dim} \int_{\Gamma_x} b_dd \phi_k d\Gamma  
dove = find(quali_vec(:,1)' == 'H' & ~ismember(quali_vec(:,3)','b'));
for jj = dove(:)'
    lettera = quali_vec(jj,3);
    eval(['H_' lettera '= matrice_H(DATA,MESH,FE_SPACE,gamma_' lettera ', DATA.x_bari);'])
end
%% Assemblo il tensore C_j, matrici di massa sul contorno di controllo \Gamma_j
% [C_j]_{i,k} = \int_{\Gamma_c} {\psi_j \phi_i \phi_k  d\Gamma}
%%% NON è detto che gamma_c sia stata definita, quindi in caso dà errore !!!
if any(strcmp(quali,'C_j'))
    C_j = massa_contorno_C_j(DATA,MESH,FE_SPACE,gamma_c);
end
%% Riordino output
for jj = 1:20
    eval(['out' num2str(jj) '= [];'])
end

for jj = 1:length(quali)
    eval(['out' num2str(jj) '=' quali{jj} ';'])
end

fprintf('Fatto in %.3f s\n',toc)
end


%% Nested
function C_j = massa_contorno_C_j(DATA,MESH,FE_SPACE,gamma)
switch MESH.dim
    case 2
        [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
        eta            =  1 - csi;
        nqn            =  length(csi);
    
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
            
        [rows,cols] = meshgrid(1:nbn,1:nbn);
        rows        = rows(:);
        cols        = cols(:);
        
        % Calcolo i valori del parametro nei punti di quadratura
        xlt = zeros(nof,nqn); ylt = xlt;
        coord_ref = [eta; csi];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
        end
        C_j{DATA.n_basi} = [];  k = C_j;
    for nn = 1:DATA.n_basi
        Crows       = zeros(nbn*nbn*nof,1);
        Ccols       = Crows;
        Ccoef       = Crows;
    
        alphaR  = DATA.psi({xlt,ylt},nn,MESH.vertices);
        
        one       = ones(nof,nqn);
        alphaR    = alphaR.*one;
    
        % trovo le lunghezze dei lati
        x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma));
        
        side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
        
        for ll = 1 : nof
            face = gamma(ll);
            alphaR_loc   = repmat(alphaR(ll,:).*wi, nbn,1);
    
            Crows(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(rows,face);
            Ccols(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(cols,face);
            Ccoef(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  side_length(ll) *(alphaR_loc.*phi)*phi';
        end
        C_j{nn} = [Crows,Ccols,Ccoef];
        k{nn} = nn * ones(nbn*nbn*nof,1);
    end
        C_j = cat(1,C_j{:});     k = cat(1,k{:});
        C_j = sptensor([C_j(:,1:2), k], C_j(:,3), [[1 1]*MESH.numNodes   DATA.n_basi]);
    


    case 3
        [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
        csi = quad_points(1,:);
        eta = quad_points(2,:);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
        eta1           =  1-csi-eta;
        nqn            =  length(wi);

        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
        coord_ref = [eta1; csi; eta];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,dof);
            zlt = zlt + vtemp'*coord_ref(j,:);
        end

        x    =  MESH.vertices(1,MESH.boundaries(1:3, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:3, gamma));
        z    =  MESH.vertices(3,MESH.boundaries(1:3, gamma));
        areav = cross( [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                       [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);

        [rows,cols] = meshgrid(1:nbn,1:nbn);
        rows        = rows(:);
        cols        = cols(:);
        
        Crows       = zeros(nbn*nbn*nof,1);
        Ccols       = Crows;
        Ccoef       = Crows;

        C_j{DATA.n_basi} = [];  k = C_j;

    for nn = 1:DATA.n_basi
        alphaR  = DATA.psi({xlt,ylt,zlt},nn,MESH.vertices);

        for ll = 1 : nof
            face = gamma(ll);
            alphaR_loc   = repmat(alphaR(ll,:).*wi, nbn,1);
            area   = 0.5*norm(areav(:,ll));
            detjac = 2*area;

            Crows(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(rows,face);
            Ccols(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(cols,face);
            Ccoef(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  detjac *(alphaR_loc.*phi)*phi';
        end
        C_j{nn} = [Crows,Ccols,Ccoef];
        k{nn} = nn * ones(nbn*nbn*nof,1);
    end
        C_j = cat(1,C_j{:});     k = cat(1,k{:});
        C_j = sptensor([C_j(:,1:2), k], C_j(:,3), [[1 1]*MESH.numNodes   DATA.n_basi]);
%         C = sparse(Crows,Ccols,Ccoef,MESH.numNodes,MESH.numNodes);

end
end


function C = massa_contorno_C(~,MESH,FE_SPACE,gamma)
switch MESH.dim
    case 2
        [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Crows       = zeros(nbn*nbn*nof,1);
        Ccols       = Crows;
        Ccoef       = Crows;
        
        [rows,cols] = meshgrid(1:nbn,1:nbn);
        rows        = rows(:);
        cols        = cols(:);
        
        x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma));
        side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
        
        coso = ( repmat(wi,nbn,1) .* phi )* phi';
        for ll = 1 : nof
            face = gamma(ll);
            Crows(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(rows,face);
            Ccols(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(cols,face);
            Ccoef(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  side_length(ll) * coso;
        end
        C = sparse(Crows,Ccols,Ccoef,MESH.numNodes,MESH.numNodes);

    case 3
        [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
        csi = quad_points(1,:);
        eta = quad_points(2,:);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
%         eta1           =  1-csi-eta;
        nqn            =  length(wi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Crows       = zeros(nbn*nbn*nof,1);
        Ccols       = Crows;
        Ccoef       = Crows;

        [rows,cols] = meshgrid(1:nbn,1:nbn);
        rows        = rows(:);
        cols        = cols(:);
                
        x    =  MESH.vertices(1,MESH.boundaries(1:3, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:3, gamma));
        z    =  MESH.vertices(3,MESH.boundaries(1:3, gamma));
        areav = cross( [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                       [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
        
        coso = ( repmat(wi,nbn,1) .* phi )* phi';
        for ll = 1 : nof
            face = gamma(ll);
            area   = 0.5*norm(areav(:,ll));
            detjac = 2*area;

            Crows(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(rows,face);
            Ccols(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  MESH.boundaries(cols,face);
            Ccoef(1+(ll-1)*nbn*nbn:ll*nbn*nbn)   =  detjac * coso;
        end
        C = sparse(Crows,Ccols,Ccoef,MESH.numNodes,MESH.numNodes);


%         u_Neumann = DATA.y_inc_fun({xlt,ylt,zlt},DATA.k_amb(ff),DATA.direz(ff));
%         for ll = 1 : nof
%             
%             area   = 0.5*norm(areav(:,ll));
%             detjac = 2*area;
%             
%             face = MESH.Neumann_side(ll);
%             
%             u_Neumann_loc  = u_Neumann(ll,:).*wi;
%             u_Neumann_loc  = u_Neumann_loc(1,:)';
%             
%             Rrows(1+(ll-1)*nbn:ll*nbn)    = MESH.boundaries(1:nbn,face);
%             Rcoef(1+(ll-1)*nbn:ll*nbn)    = detjac*phi*u_Neumann_loc;
%         end
%         
end
end


function F = forzante_F(DATA,MESH,FE_SPACE,gamma)
switch MESH.dim
    case 2
        normali = MESH.Normal_Faces(:,gamma);
    %
        [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
        eta            =  1 - csi;
        nqn            =  length(csi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Rrows       = zeros(nbn*nof,1);
        Rcoef       = Rrows;
        
        xlt = zeros(nof,nqn); ylt = xlt;
        coord_ref = [eta; csi];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
        end
        
        x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma));
        side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
        
        F{DATA.n_frq} = [];
        for ff = 1:DATA.n_frq
            u_Robin_x = DATA.dy_inc_dx_fun({xlt,ylt},DATA.k_amb(ff),DATA.direz(ff));
            u_Robin_y = DATA.dy_inc_dy_fun({xlt,ylt},DATA.k_amb(ff),DATA.direz(ff));
        
            u_Robin = u_Robin_x .* repmat(normali(1,:)',1,nqn) + u_Robin_y .* repmat(normali(2,:)',1,nqn);
        
            for ll = 1 : nof
                face = gamma(ll);
                u_Robin_loc  = u_Robin(ll,:).*wi;
        %         u_Robin_loc  = u_Robin_loc(1,:).';
                Rrows(1+(ll-1)*nbn:ll*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(ll-1)*nbn:ll*nbn)    = side_length(ll)*phi*(u_Robin_loc).';
            end
            F{ff} = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
        end

    case 3
        [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
        csi = quad_points(1,:);
        eta = quad_points(2,:);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
        eta1           =  1-csi-eta;
        nqn            =  length(wi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Rrows       = zeros(nbn*nof,1);
        Rcoef       = Rrows;
        
        xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
        coord_ref = [eta1; csi; eta];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,dof);
            zlt = zlt + vtemp'*coord_ref(j,:);
        end
        
        x    =  MESH.vertices(1,MESH.boundaries(1:3, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:3, gamma));
        z    =  MESH.vertices(3,MESH.boundaries(1:3, gamma));
        
        areav = cross( [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                       [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
        
        normali = MESH.Normal_Faces(:,gamma);

        F{DATA.n_frq} = [];
        for ff = 1:DATA.n_frq
            u_Robin_x = DATA.dy_inc_dx_fun({xlt,ylt,zlt},DATA.k_amb(ff),DATA.direz(ff));
            u_Robin_y = DATA.dy_inc_dy_fun({xlt,ylt,zlt},DATA.k_amb(ff),DATA.direz(ff));
            u_Robin_z = DATA.dy_inc_dz_fun({xlt,ylt,zlt},DATA.k_amb(ff),DATA.direz(ff));

            u_Neumann = u_Robin_x .* repmat(normali(1,:)',1,nqn) ...
                      + u_Robin_y .* repmat(normali(2,:)',1,nqn) ...
                      + u_Robin_z .* repmat(normali(3,:)',1,nqn);
        
            for l = 1 : nof
                
                area   = 0.5*norm(areav(:,l));
                detjac = 2*area;
                
                face = gamma(l);
                
                u_Neumann_loc  = u_Neumann(l,:).*wi;
                u_Neumann_loc  = u_Neumann_loc(1,:)';
                
                Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(l-1)*nbn:l*nbn)    = detjac*phi*u_Neumann_loc;
            end
            F{ff} = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
        end
end
end


function G = forzante_G(DATA,MESH,FE_SPACE,gamma)
switch MESH.dim
    case 2
        [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
        eta            =  1 - csi;
        nqn            =  length(csi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Rrows       = zeros(nbn*nof,1);
        Rcoef       = Rrows;
        
        xlt = zeros(nof,nqn); ylt = xlt;
        coord_ref = [eta; csi];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
        end
        
        x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma));
        side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
        
        G{DATA.n_frq} = [];
        for ff = 1:DATA.n_frq
            u_Robin = DATA.y_inc_fun({xlt,ylt},DATA.k_amb(ff),DATA.direz(ff));
        
            for ll = 1 : nof
                face = gamma(ll);
                u_Robin_loc  = u_Robin(ll,:).*wi;
        %         u_Robin_loc  = u_Robin_loc(1,:).';
                Rrows(1+(ll-1)*nbn:ll*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(ll-1)*nbn:ll*nbn)    = side_length(ll)*phi*(u_Robin_loc).';
            end
            G{ff} = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
        end

    case 3

        [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
        csi = quad_points(1,:);
        eta = quad_points(2,:);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
        eta1           =  1-csi-eta;
        nqn            =  length(wi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Rrows       = zeros(nbn*nof,1);
        Rcoef       = Rrows;
        
        xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
        coord_ref = [eta1; csi; eta];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,dof);
            zlt = zlt + vtemp'*coord_ref(j,:);
        end
        
        x    =  MESH.vertices(1,MESH.boundaries(1:3, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:3, gamma));
        z    =  MESH.vertices(3,MESH.boundaries(1:3, gamma));
        areav = cross( [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                       [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
        
        G{DATA.n_frq} = [];
        for ff = 1:DATA.n_frq
            u_Neumann = DATA.y_inc_fun({xlt,ylt,zlt},DATA.k_amb(ff),DATA.direz(ff));
            if ~all(size(u_Neumann) == [nof,nqn])
                error('!!!')
            end

            for ll = 1 : nof
                
                area   = 0.5*norm(areav(:,ll));
                detjac = 2*area;
                
                face = gamma(ll);
                
                u_Neumann_loc  = u_Neumann(ll,:).*wi;
                u_Neumann_loc  = u_Neumann_loc(1,:)';
                
                Rrows(1+(ll-1)*nbn:ll*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(ll-1)*nbn:ll*nbn)    = detjac*phi*u_Neumann_loc;
            end
            G{ff} = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
        end

end
end


function H = matrice_H(~,MESH,FE_SPACE,gamma,x_G)
switch MESH.dim
    case 2
        normali = MESH.Normal_Faces(:,gamma);
    %
        [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
        eta            =  1 - csi;
        nqn            =  length(csi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Rrows       = zeros(nbn*nof,1);
        Rcoef       = Rrows;
        
        xlt = zeros(nof,nqn); ylt = xlt;
        coord_ref = [eta; csi];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
        end
        
        x    =  MESH.vertices(1,MESH.boundaries(1:2, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:2, gamma));
        side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
        
        xx_G{1} = xlt - x_G(1);
        xx_G{2} = ylt - x_G(2);

        nn{1} = repmat(normali(1,:)',1,nqn);
        nn{2} = repmat(normali(2,:)',1,nqn);

        momento = xx_G{1}.*nn{2} - xx_G{2}.*nn{1};

        u_Robin = momento;

        for ll = 1 : nof
            face = gamma(ll);
            u_Robin_loc  = u_Robin(ll,:).*wi;

            Rrows(1+(ll-1)*nbn:ll*nbn)    = MESH.boundaries(1:nbn,face);
            Rcoef(1+(ll-1)*nbn:ll*nbn)    = side_length(ll)*phi*(u_Robin_loc).';
        end
        H = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
    
    case 3
        [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
        csi = quad_points(1,:);
        eta = quad_points(2,:);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
        eta1           =  1-csi-eta;
        nqn            =  length(wi);
        
        nof         = length(gamma);
        nbn         = MESH.numBoundaryDof;
        
        Rrows       = zeros(nbn*nof,1);
        Rcoef       = Rrows;
        
        xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
        coord_ref = [eta1; csi; eta];
        for j = 1 : 2
            dof = MESH.boundaries(j,gamma);
            vtemp = MESH.vertices(1,dof);
            xlt = xlt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,dof);
            ylt = ylt + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,dof);
            zlt = zlt + vtemp'*coord_ref(j,:);
        end
        
        x    =  MESH.vertices(1,MESH.boundaries(1:3, gamma));
        y    =  MESH.vertices(2,MESH.boundaries(1:3, gamma));
        z    =  MESH.vertices(3,MESH.boundaries(1:3, gamma));
        
        areav = cross( [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                       [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
        
        normali = MESH.Normal_Faces(:,gamma);
        

        xx_G{1} = xlt - x_G(1);
        xx_G{2} = ylt - x_G(2);
        xx_G{3} = zlt - x_G(3);
        nn{1} = repmat(normali(1,:)',1,nqn);
        nn{2} = repmat(normali(2,:)',1,nqn);
        nn{3} = repmat(normali(3,:)',1,nqn);
        momento{1} = xx_G{2}.*nn{3} - xx_G{3}.*nn{2};
        momento{2} = xx_G{3}.*nn{1} - xx_G{1}.*nn{3};
        momento{3} = xx_G{1}.*nn{2} - xx_G{2}.*nn{1};

        H = sparse(MESH.numNodes,MESH.numNodes);
        for dd = 1:MESH.dim
            u_Neumann = momento{dd};
        
            for l = 1 : nof
                
                area   = 0.5*norm(areav(:,l));
                detjac = 2*area;
                
                face = gamma(l);
                
                u_Neumann_loc  = u_Neumann(l,:).*wi;
                u_Neumann_loc  = u_Neumann_loc(1,:)';
                
                Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(l-1)*nbn:l*nbn)    = detjac*phi*u_Neumann_loc;
            end
            H_loc = sparse(Rrows,1,Rcoef,MESH.numNodes,1);
            H = H + H_loc*H_loc';
        end
end
end

