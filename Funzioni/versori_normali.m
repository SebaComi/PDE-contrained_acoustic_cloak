function normali = versori_normali(MESH)
boundaries = MESH.boundaries;
vertices = MESH.vertices;
nodi_boundaries = MESH.nodi_boundaries;
% normali{dd,cc} è una matrice 3xn che ha per colonne i versori normali al
    % contorno cc del dominio dd nei suoi nodi e orientati verso l'esterno del dominio
    % la prima riga dei versori è l'indice del nodo a cui è riferito
% Se il contorno cc non è un contorno di dd, normali{dd,cc} = []
n_domini = size(nodi_boundaries,2);
n_contorni = max(boundaries(MESH.bound_tipo_cont,:));
normali{n_domini,n_contorni} = [];

for dd = 1:n_domini
    loc_bound = boundaries(:,nodi_boundaries(:,dd));
    
    cc_del_dd = unique(loc_bound(MESH.bound_tipo_cont,:));
    inizi = zeros(length(cc_del_dd),1);  % code contiene l'indice dei vertici che sono coda di contorno
        ii = 0;
    for cc = cc_del_dd(:)'
        
        ii = ii+1;
        
        loc_loc_bound = loc_bound(:,loc_bound(MESH.bound_tipo_cont,:) == cc);
        % riordino rispetto all'ascissa curvilinea
        %%% aggiornamento 14/07/21:
        % guardo l'ascissa curvilinea sulla riga 4 e non 3 perché redbKit
        % sovrascrive la riga 3 se uso l'ordine P2 degli elementi...
        [~,index] = sort(loc_loc_bound(4,:));
        loc_loc_bound = loc_loc_bound(:,index);
        
        inizi(ii) = loc_loc_bound(1,1);
        
        % Decido il verso di rotazione in base a:
        if loc_loc_bound(7,1) == dd      % se il dominio dd è a destra, ruoto a sisnistra
            ruoto = [0 -1; 1 0];
        elseif loc_loc_bound(6,1) == dd  % se il dominio dd è a sinistra, ruoto a destra
            ruoto = [0 1; -1 0];
        else
            error('C''è un problema')
        end
        
        x1 = vertices(1,loc_loc_bound(1,:));
        x2 = vertices(1,loc_loc_bound(2,:));
        y1 = vertices(2,loc_loc_bound(1,:));
        y2 = vertices(2,loc_loc_bound(2,:));
        % resta il problema del segno
        vect = [0 x2-x1 0; 0 y2-y1 0];   % Sono i vettori dei lati orlati di zeri
        vect = vect(:,1:end-1) + vect(:,2:end); % Sono i vettori delle diagonali che uniscono i nodi alterni,
                                                % a parte inizio e fine che restano i vect dei lati.
                                                % in questo modo le BC nei vertici hanno contributo
                                                % metà da un lato e metà dall'altro.
                                                % Resta il problema della normalizzazione sui vertici
        vect = ruoto * vect;
        norma = vecnorm(vect);
        vect = vect ./ [norma; norma];
        
        normali{dd,cc} = [loc_loc_bound(1,:) loc_loc_bound(2,end);   vect];
        
    end
   %% Qui si risolve il problema della normalizzazione sui vertici 
%    ********* OLD *********
%    inizi = unique(inizi);
%     for jj = code(:)'
%         vect = [];
%         for cc = cc_del_dd(:)'
%             dove = normali{dd,cc}(1,:) == jj;
%             vect = [vect, normali{dd,cc}(2:3,dove)];
%         end
%         
%         if size(vect,2) ~= 2    % Ci devono essere sempre e solo due contorni con la stessa coda jj
%             error('Errore')
%         end
%         norma = norm(sum(vect,2));
%         
%         for cc = cc_del_dd(:)'
%             dove = normali{dd,cc}(1,:) == jj;
%             normali{dd,cc}(2:3,dove) = normali{dd,cc}(2:3,dove);% / norma; % /norma affinché la norma della somma sia sempre 1
%         end
%     end
end

