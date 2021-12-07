function [A_in, F_in, y_Dirichlet] = App_BC(A, F, MESH, DATA)
% Applica le BCs

y_Dirichlet = DATA.bcDir;

F_in = F(MESH.internal_dof,:)  -  A(MESH.internal_dof,MESH.Dirichlet_dof) * y_Dirichlet;

A_in = A(MESH.internal_dof,MESH.internal_dof);

end