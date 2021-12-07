function gamma_in = trova_gamma_in(DATA,MESH,nome)
dirichlet = [];
for chi = find(ismember(DATA.tipo_contorno, nome))
    dove = find(MESH.boundaries(MESH.bound_tipo_cont,:) == chi);
    dirichlet = union(dirichlet, MESH.boundaries(1:MESH.numBoundaryDof,dove) );
end
gamma_in = false(MESH.numNodes,1);
gamma_in(dirichlet(:)) = true;