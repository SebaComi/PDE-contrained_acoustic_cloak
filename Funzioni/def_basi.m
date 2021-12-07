function [DATA] = def_basi(DATA,MESH)
% Definizione basi
switch DATA.dim
    case 2
        switch DATA.tipo_base
            case 'vertici'
                DATA.n_basi = MESH.n_vertices;
                DATA.ctrl2nodi = (1:MESH.n_vertices)';
            case 'rettangoli'
                    DATA.nodi_basi = isbase(DATA,MESH);
                    DATA.n_basi = length(DATA.nodi_basi);
                    
                    DATA.psi = @(xy,n,vertices) rettangolo(xy,n,DATA); % ,idx_bordo_mant,idx_bordo_osta
                    
%                     DATA.ctrl2nodi = zeros(MESH.n_vertices,1);
%                     for jj = 1:MESH.n_vertices
%                         dove = MESH.vertices(1,jj) >= DATA.nodi_basi(1,:)               & ...
%                                MESH.vertices(1,jj) <  DATA.nodi_basi(1,:) + DATA.lato_x -100*eps & ...  % riduco leggermente i quadretti altrimenti capita che dei vertici appartengano a due quadretti
%                                MESH.vertices(2,jj) >= DATA.nodi_basi(2,:)               & ...
%                                MESH.vertices(2,jj) <  DATA.nodi_basi(2,:) + DATA.lato_y -100*eps ;
%                        if any(dove)
%                            DATA.ctrl2nodi(jj) = find(dove);
%                        else
%                            DATA.ctrl2nodi(jj) = 0;
%                        end
%                     end
                    
                    DATA.area = DATA.lato_x * DATA.lato_y;
            case 'corone circolari'
                    DATA.nodi_basi = isbase(DATA,[]);
                    DATA.n_basi = length(DATA.nodi_basi);
                    
                    DATA.psi = @(xy,n,vertices) corona_circ(xy,n,DATA);
                    
                    DATA.ctrl2nodi = zeros(MESH.numVertices,1);
                    vect = MESH.vertices(1,:) + 1i * MESH.vertices(2,:);
                    r = abs(vect);
                    theta = angle(vect);
                    theta = mod(theta,2*pi);
                    
                    for jj = 1:MESH.numVertices
                        dove = r(jj) >= DATA.nodi_basi(1,:)                     & ...
                               r(jj) <  DATA.nodi_basi(1,:) + DATA.dr -100*eps  & ...  % riduco leggermente i quadretti altrimenti capita che dei vertici appartengano a due quadretti
                               theta(jj) >= DATA.nodi_basi(2,:)                 & ...
                               theta(jj) <  DATA.nodi_basi(2,:) + DATA.dtheta -100*eps ;
                       if any(dove)
                           DATA.ctrl2nodi(jj) = find(dove);
                       else
                           DATA.ctrl2nodi(jj) = 0;
                       end
                    end
                    DATA.area = (DATA.nodi_basi(1,:)' + DATA.dr/2) * DATA.dtheta * DATA.dr;
            case 'hex'
    %         Gli esagoni hanno due lati orizzontali !!!
                    DATA.nodi_basi = isbase(DATA,MESH);
                    DATA.n_basi = length(DATA.nodi_basi);
                    
                    DATA.psi = @(xy,n,vertices) esagono(xy,n,DATA);
                    
                    %%% Per ogni vertice, definisco a quale base appartiene
            
                    DATA.ctrl2nodi = zeros(MESH.numVertices,1);
                    vert = MESH.vertices';
                    
                    for nn = 1:DATA.n_basi
                        dove = logical( esagono({vert(:,1),vert(:,2)},nn,DATA) );
                        
                        DATA.ctrl2nodi(dove) = nn;
                    end
                    
                    DATA.area = DATA.lato_hex^2 * 3*sqrt(3)/2;
            case 'lati'
                DATA.nodi_basi = [];
                dx = DATA.dx;
                for jj = find(DATA.tipo_contorno == 'c')
                    x1 = min(min(MESH.vertices(1,MESH.boundaries(1:2, MESH.boundaries(5,:) == jj))));
                    x2 = max(max(MESH.vertices(1,MESH.boundaries(1:2, MESH.boundaries(5,:) == jj))));
                    DATA.nodi_basi = [ DATA.nodi_basi, x1:dx:(x2-dx) ];
                end
                DATA.n_basi = length(DATA.nodi_basi);
                DATA.psi = @(xy,n,vertices) (xy{1} >= DATA.nodi_basi(n)) & (xy{1} < DATA.nodi_basi(n)+dx) & (abs(xy{2}) < 10*eps );
            otherwise
                error('Ciao!')
        end

    case 3
        switch DATA.tipo_base
            case 'rettangoli'
                DATA.nodi_basi = isbase(DATA,MESH);
                DATA.n_basi = length(DATA.nodi_basi);
                
                DATA.psi = @(xyz,n,vertices) rettangolo(xyz,n,DATA) .* (xyz{3} == 0); % ,idx_bordo_mant,idx_bordo_osta
                
%                 DATA.ctrl2nodi = zeros(MESH.n_vertices,1);
%                 for jj = 1:MESH.n_vertices
%                     dove = MESH.vertices(1,jj) >= DATA.nodi_basi(1,:)               & ...
%                            MESH.vertices(1,jj) <  DATA.nodi_basi(1,:) + DATA.lato_x -100*eps & ...  % riduco leggermente i quadretti altrimenti capita che dei vertici appartengano a due quadretti
%                            MESH.vertices(2,jj) >= DATA.nodi_basi(2,:)               & ...
%                            MESH.vertices(2,jj) <  DATA.nodi_basi(2,:) + DATA.lato_y -100*eps ;
%                    if any(dove)
%                        DATA.ctrl2nodi(jj) = find(dove);
%                    else
%                        DATA.ctrl2nodi(jj) = 0;
%                    end
%                 end
                DATA.area = DATA.lato_x * DATA.lato_y;

            case 'corone circolari'
                DATA.nodi_basi = isbase(DATA,[]);
                DATA.n_basi = length(DATA.nodi_basi);
                
                DATA.psi = @(xyz,n,vertices) corona_circ(xyz,n,DATA);
                
                % La parte commentata non funziona per il caso 3D !!!
                
%                 DATA.ctrl2nodi = zeros(MESH.numVertices,1);
%                 vect = MESH.vertices(1,:) + 1i * MESH.vertices(2,:);
%                 r = abs(vect);
%                 theta = angle(vect);
%                 theta = mod(theta,2*pi);
                
%                 for jj = 1:MESH.numVertices
%                     dove = r(jj) >= DATA.nodi_basi(1,:)                     & ...
%                            r(jj) <  DATA.nodi_basi(1,:) + DATA.dr -100*eps  & ...  % riduco leggermente i quadretti altrimenti capita che dei vertici appartengano a due quadretti
%                            theta(jj) >= DATA.nodi_basi(2,:)                 & ...
%                            theta(jj) <  DATA.nodi_basi(2,:) + DATA.dtheta -100*eps ;
%                    if any(dove)
%                        DATA.ctrl2nodi(jj) = find(dove);
%                    else
%                        DATA.ctrl2nodi(jj) = 0;
%                    end
%                 end
                DATA.area = (DATA.nodi_basi(1,:)' + DATA.dr/2).^2 * DATA.dtheta^2 * DATA.dr;
            otherwise
                error('Ciao!')
        end
end

