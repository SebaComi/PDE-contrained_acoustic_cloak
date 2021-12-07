function nodi_basi = isbase(DATA,MESH)
% (idx_bordo_mant,idx_bordo_osta,DATA,MESH)

switch DATA.dim
    case 2
        switch DATA.tipo_base
            case 'rettangoli'
%                 lato_q = DATA.lato_q;
%                 [x,y] = meshgrid(-lato_q:DATA.lato_x:lato_q,   -lato_q:DATA.lato_y:lato_q);
                lato_max = max(max(abs(MESH.vertices(1:2,:))));
                [x,y] = meshgrid(-lato_max:DATA.lato_x:lato_max,   -lato_max:DATA.lato_y:lato_max);
        
                x = x(:);
                y = y(:);
        
                % impongo che siano contenuti i punti (x;y), (x+lato;y), (x;y+lato) e
                % (x+lato;y+lato), in questo modo restano solo i vertici in basso a sx dei
                % quadrati interamente contenuti nel mantello
                in = true(size(x));
                quale_dominio = find(DATA.tipo_dominio == 'c');
                for D = [0 0 1 1
                         0 1 0 1]
                    xx = x+D(1)*DATA.lato_x; yy = y+D(2)*DATA.lato_y;
                    in = in &  trova_dominio([xx(:),yy(:)]',MESH,quale_dominio)'  ;
                end
                nodi_basi = [x(in), y(in)].';   % sono i nodi in basso a sx di ogni rettangolo

            case 'corone circolari'
                r = DATA.raggio_osta;
                R = DATA.raggio_mant;
                dtheta = DATA.dtheta;
                dr = DATA.dr;
                [rr, theta] = meshgrid(r:dr:R-dr,    0:dtheta:(2*pi-100*eps));
                nodi_basi = [rr(:), theta(:)].'; % sono i nodi verso il centro dell'ostacolo e ad angolo minore, quindi senso orario
        %         La base quindi si sviluppa fra [rr(ii); rr(ii)+dr] x [theta(ii); theta(ii)+dtheta] 
            case 'hex'
            %         Gli esagoni hanno due lati orizzontali !!!
                lato = DATA.lato_hex;
                range = 2*DATA.raggio_mant;
                asse_x = [flip(0 : -3*lato : -range)        3*lato : 3*lato : range];
                asse_y = [flip(0 : -sqrt(3)*lato : -range)      sqrt(3)*lato : sqrt(3)*lato : range];
                [x,y] = meshgrid(asse_x, asse_y);
                x = [x(:); x(:)+1.5*lato];
                y = [y(:); y(:)+sqrt(3)/2*lato];
                
                % impongo che siano contenuti i 6 vertici degli esagoni, in questo
                % modo restano solo i centri degli esagoni
                % interamente contenuti nel mantello
                % Ho aggiunto anche i punti medi dei lati, spero che restino solo
                % gli esagoni che servono
                
                in = true(size(x));
                vert_esagoni = lato*[1     1/2         -1/2        -1      -1/2        1/2          3/4         0           -3/4        -3/4        0           3/4
                                     0   sqrt(3)/2   sqrt(3)/2      0   -sqrt(3)/2   -sqrt(3)/2     sqrt(3)/4   sqrt(3)/2   sqrt(3)/4   -sqrt(3)/4  -sqrt(3)/2  -sqrt(3)/4];
                quale_dominio = find(DATA.tipo_dominio == 'c');
                for D = vert_esagoni
                    xx = x+D(1); yy = y+D(2);
        %             FaceID = nearestFace(MESH.gm,[xx(:),yy(:)]);
        %             in = in &  (FaceID == quale_dominio)' ;
                    in = in &  trova_dominio([xx(:),yy(:)]',MESH,quale_dominio)' ;
        %             in = in & isdentro(x+D(1),y+D(2), idx_bordo_mant,idx_bordo_osta,MESH.vertices);
                end
                nodi_basi = [x(in), y(in)].';   % sono i centri degli esagoni
    
            otherwise
                error('Ciao!')
        end

    case 3
        switch DATA.tipo_base
            case 'rettangoli'
%                 lato_q = DATA.lato_q;
                lato_max = max(max(abs(MESH.vertices(1:2,:))));
                [x,y] = meshgrid(-lato_max:DATA.lato_x:lato_max,   -lato_max:DATA.lato_y:lato_max);
        
                x = x(:);
                y = y(:);
        
                % impongo che siano contenuti i punti (x;y), (x+lato;y), (x;y+lato) e
                % (x+lato;y+lato), in questo modo restano solo i vertici in basso a sx dei
                % quadrati interamente contenuti nel mantello
                in = true(size(x));
                quale_contorno = find(DATA.tipo_contorno == 'c');
                for D = [ [0 0 1 1]*DATA.lato_x
                          [0 1 0 1]*DATA.lato_y ]
                    xx = x+D(1); yy = y+D(2);
                    in = in &  trova_dominio([xx(:),yy(:)]', MESH.MESH_2, quale_contorno)'  ;
                end
                nodi_basi = [x(in), y(in)].';   % sono i nodi in basso a sx di ogni rettangolo
            case 'corone circolari'
                r = DATA.raggio_osta;
                R = DATA.raggio_mant;
                dtheta = DATA.dtheta;
                dr = DATA.dr;
                [rr, theta1,theta2] = meshgrid(r:dr:R-dr,   0:dtheta:(2*pi-100*eps),   0:dtheta:(pi-100*eps));
                nodi_basi = [rr(:), theta1(:), theta2(:)].'; % sono i nodi verso il centro dell'ostacolo e ad angolo minore
    %         La base quindi si sviluppa fra [rr(ii); rr(ii)+dr]  x  [theta1(ii); theta1(ii)+dtheta]  x  [theta2(ii); theta2(ii)+dtheta]
            otherwise
                error('Ciao!')
        end
end
