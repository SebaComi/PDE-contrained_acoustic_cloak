function [] = disegno_mesh(DATA,MESH)
    % Disegno mesh
    fprintf('Disegno mesh. '); tic
    % Visualizziamo la mesh, i contorni e le normali
    % figure(2)
    switch DATA.dim
        case 2
            pdeplot(MESH.vertices,MESH.elements(1:3,:));
        case 3
            pdeplot3D(MESH.vertices,MESH.elements(1:4,:),'FaceAlpha',0);
    end
    hold on
    colori = get(gca,'ColorOrder'); colori = [colori;colori];
    colori(3,:) = [];
    % Domini e normali
    for dd = 1:MESH.n_domini
        if any(MESH.nodi_elements(:,dd))
            switch DATA.dim
                case 2
                    % Elementi
                    h = pdeplot(MESH.vertices,MESH.elements(1:3,MESH.nodi_elements(:,dd)));
                    h.Color = colori(dd,:);
                    % Normali
                    xx = MESH.vertices(1,MESH.boundaries(1:2,MESH.nodi_boundaries(:,dd)));
                    xx = mean(reshape(xx,2,[]),1);
                    yy = MESH.vertices(2,MESH.boundaries(1:2,MESH.nodi_boundaries(:,dd)));
                    yy = mean(reshape(yy,2,[]),1);
                    hold on
                    quiver(xx,yy,MESH.Normal_Faces(1,MESH.nodi_boundaries(:,dd)), ...
                                 MESH.Normal_Faces(2,MESH.nodi_boundaries(:,dd)), ...
                           'Color',colori(dd,:), 'LineWidth',1,'AutoScale','on')
                case 3
                    % Elementi
                    h = pdeplot3D(MESH.vertices,MESH.elements(1:4,MESH.nodi_elements(:,dd)),'FaceAlpha',0.5);
                    h.EdgeColor = colori(dd,:);
                    % Normali
                    xx = MESH.vertices(1,MESH.boundaries(1:3,MESH.nodi_boundaries(:,dd)));
                    xx = mean(reshape(xx,3,[]),1);
                    yy = MESH.vertices(2,MESH.boundaries(1:3,MESH.nodi_boundaries(:,dd)));
                    yy = mean(reshape(yy,3,[]),1);
                    zz = MESH.vertices(3,MESH.boundaries(1:3,MESH.nodi_boundaries(:,dd)));
                    zz = mean(reshape(zz,3,[]),1);
                    hold on
                    quiver3(xx,yy,zz, ...
                        MESH.Normal_Faces(1,MESH.nodi_boundaries(:,dd)), ...
                        MESH.Normal_Faces(2,MESH.nodi_boundaries(:,dd)), ...
                        MESH.Normal_Faces(3,MESH.nodi_boundaries(:,dd)), ...
                        'Color',colori(dd,:), 'LineWidth',1,'AutoScale','on')
            end
        end
    end
switch DATA.dim
    case 2
        % Contorni
    for cc = 1:MESH.n_contorni
        dove = MESH.boundaries(5,:) == cc;
        if any(dove)
            linea = MESH.boundaries(:,dove);
            [~,index] = sort(linea(4,:));   % aggiornamento 14/07/21:
                                            % guardo l'ascissa curvilinea sulla riga 4 e non 3 perché redbKit
                                            % sovrascrive la riga 3...
            linea = linea(:,index);
            plot(MESH.vertices(1,[linea(1,:) linea(2,end)]),  MESH.vertices(2,[linea(1,:) linea(2,end)]), ...
                'k','LineWidth',.8)
        end
    end
end
    drawnow
    fprintf('             Fatto in %.3f s\n',toc)
end

