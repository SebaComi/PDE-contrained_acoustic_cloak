function [vertices, boundaries, elements, edges, DATA] = mesh_3D_Comsol(DATA)

model = DATA.model;

[~,meshdata] = mphmeshstats(model);
% info = mphxmeshinfo(model);
vertices = meshdata.vertex;

dove = strcmp(meshdata.types,'tet');
elements = [ double(meshdata.elem{dove}+1);
            meshdata.elementity{dove}'     ];
elements = double(elements);

dove = strcmp(meshdata.types,'tri');
boundaries = [ meshdata.elem{dove}+1;
               zeros(size(meshdata.elem{dove}));
               meshdata.elementity{dove}'       ];
boundaries = double(boundaries);

dove = strcmp(meshdata.types,'edg');
edges = [ meshdata.elem{dove}+1;
          zeros(2, size(meshdata.elem{dove},2));
          meshdata.elementity{dove}' ];
edges = double(edges);

ciao = mphxmeshinfo(model);
nodes = ciao.nodes.coords;
idx = zeros(size(vertices,2),1);

for jj = 1:length(idx)
    i1 = abs(vertices(1,:) - nodes(1,jj)) < 10*eps;
    i2 = abs(vertices(2,:) - nodes(2,jj)) < 10*eps;
    i3 = abs(vertices(3,:) - nodes(3,jj)) < 10*eps;
%     if jj == 9903
%         ciao
%     end
    idx(jj) = find(i1 & i2 & i3);
end

DATA.idx_vert = idx;

end