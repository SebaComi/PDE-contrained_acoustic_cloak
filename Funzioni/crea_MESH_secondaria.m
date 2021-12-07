function [MESH_2] = crea_MESH_secondaria(MESH,DATA,nomi_contorno,dim_eliminare)


nodi_estratti = unique(cat(1,MESH.nodi_contorni{ismember(DATA.tipo_contorno, nomi_contorno)}));
vertices = MESH.vertices(setdiff(1:3,dim_eliminare),nodi_estratti);

dove = all(ismember(MESH.boundaries(1:3,:),nodi_estratti));
elements = MESH.boundaries([1:3 MESH.bound_tipo_cont],dove);
[~,ordine3] = ismember(elements(1:3,:),nodi_estratti);
elements(1:3,:) = ordine3;

dove = all(ismember(MESH.edges(1:2,:),nodi_estratti));
boundaries = MESH.edges(:,dove);
% error('C''è da risolvere la cosa del plot dei boundaries')

[~,ordine4] = ismember(boundaries(1:2,:),nodi_estratti);
if ~all(ordine4)
    error('C''è uno zero!')
end
boundaries(1:2,:) = ordine4;

MESH_2 = buildMESH( 2, elements, vertices, boundaries, DATA.fem, DATA.quad_order,[],'CFD');

MESH_2.nodi_estratti = nodi_estratti;
