function [nodi_domini,nodi_contorni,nodi_boundaries,nodi_elements] = distribuzione_domini(MESH,DATA)
vertices = MESH.vertices;
boundaries = MESH.boundaries;
elements = MESH.elements;
n_domini = MESH.n_domini;

% nodi_domini(ii,jj) = 1   sse il nodo ii appartiene al dominio jj
% elementi(ii,jj) = 1   sse l'elemento ii appartiene al dominio jj
% contorni(ii,jj) = 1   sse il contorno ii appartiene al dominio jj
n_vertices = size(vertices,2);
n_nodi_contorni = size(boundaries,2);
n_elements = size(elements,2);
n_contorni = length(DATA.tipo_contorno);

nodi_domini = false(n_vertices,n_domini);
nodi_elements = false(n_elements,n_domini);
nodi_boundaries = false(n_nodi_contorni,n_domini);

for dd = 1:n_domini
    
    a = elements(end,:) == dd;
    nodi_elements(a, dd) = true;
    
    b = elements(1:end-1, a);
    b = unique(b(:));
    nodi_domini(b,dd) = true;
    
    nodi_boundaries(:,dd) = all(ismember(boundaries(1:MESH.dim,:),b),1);
%     switch MESH.dim
%         case 2
% %             nodi_boundaries(:,dd) = boundaries(6,:) == dd | boundaries(7,:) == dd;
%         warning('Ciao!')
%                         nodi_boundaries(:,dd) = ismember(boundaries(1,:),b);
%         case 3
%             warning('Ciao!')
%             nodi_boundaries(:,dd) = ismember(boundaries(1,:),b);
%     end
end

nodi_contorni{n_contorni} = [];
% warning('Ho tolto questo riordinamento!!')
% [~,index] = sort(boundaries(3,:));  % Riordino i boundaries rispetto all'ascissa curvilinea
% boundaries = boundaries(:,index);   % In questo modo nodi_contorni è ordinato come il vettore normali
for cc = 1:n_contorni
    dove = boundaries(MESH.bound_tipo_cont,:) == cc;
    nodi_contorni{cc} = unique(boundaries(1:MESH.numBoundaryDof,dove),'stable');
end




