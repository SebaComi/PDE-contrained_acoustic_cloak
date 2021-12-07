function [vertices, boundaries, elements, dl,DATA] = mesh_WW(DATA)

% https://it.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html
switch DATA.tipo_mesh
    case 1 % 
        L = 1;
        R = 0.05;
        h = DATA.h;
dl =  [ 2   2   2   2   2   2   2   2
       -L  -L   L   1   R   R  -R  -R
       -L   L   L   R   R  -R  -R  -L
        0  -h  -h   0   0  -R  -R   0
       -h  -h   0   0  -R  -R   0   0
        1   1   1   1   1   1   1   1
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0];
    dl(2:3,:) = dl(2:3,:)*4;
    case 2 % Gamma controllo 
        L = 1;
        R = 0.1;
        C = 0.2;
        h = DATA.h;
dl =  [ 2   2   2   2   2   2   2   2   2   2
       -L  -L   L   1  R+C  R   R  -R  -R -R-C
       -L   L   L  R+C  R   R  -R  -R -R-C -L
        0  -h  -h   0   0   0  -R  -R   0   0
       -h  -h   0   0   0  -R  -R   0   0   0
        1   1   1   1   1   1   1   1   1   1
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0];
    dl(2:3,:) = dl(2:3,:)*2.5;
end
% figure

pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
[vertices, boundaries, elements] = initmesh(dl,'Jiggle','minimum','Hgrad',1.05,'Hmax',DATA.hh_mesh);
% [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,find(DATA.tipo_dominio == 'c',1));
% switch DATA.tipo_mesh
%     case 7
%         [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,1);
% %     case 9
% %         [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,)
% end
% [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,2);
% if any(boundaries(5,:) == 3)
%     [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,1);
% end

DATA.gm = [];
return

%% Cosa che serve per la funzione isbase.m e Genera_matrici_parfor.m
model = createpde();
mask = false(size(dl));
mask(6:7,interno_mantello) = true;
mask = mask & (dl == 0);
dl(mask) = max(max(dl(6:7,:)))+1;
% [vertices_fake, ~, elements_fake] = initmesh(dl);%,'Hmax',DATA.hh_mesh/2);
gm = geometryFromEdges(model,dl);
% [gm,~] = geometryFromMesh(model,vertices_fake,elements_fake(1:3,:),elements_fake(4,:));
% pg = geometryFromEdges(model,dl);
DATA.gm = gm;
