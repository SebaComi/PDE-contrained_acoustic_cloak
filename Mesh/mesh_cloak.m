function [vertices, boundaries, elements, dl,DATA] = mesh_cloak(DATA)

% https://it.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html
switch DATA.tipo_mesh
    case 1  %  Dominio quadrato nuovo
%       1   2   3   4   5   6   7   8     9  10  11  12    13  14  15  16    17  18  19  20    21  22  23  24
dl =  [ 2   2   2   2   2   2   2   2     2   2   2   2     1   1   1   1     2   2   2   2     1   1   1   1
        3   0  -3  -3  -3   0   3   3     0  -3   0   2     2   0  -2   0     0  -2   0   1     1   0  -1   0
        0  -3  -3  -3   0   3   3   3     0  -2   0   3     0  -2   0   2     0  -1   0   2     0  -1   0   1
        3   3   3   0  -3  -3  -3   0     2   0  -3   0     0   2   0  -2     1   0  -2   0     0   1   0  -1
        3   3   0  -3  -3  -3   0   3     3   0  -2   0     2   0  -2   0     2   0  -1   0     1   0  -1   0
        3   2   2   1   1   4   4   3     2   2   1   3     7   6   5   8     6   6   5   7     0   0   0   0
        0   0   0   0   0   0   0   0     3   1   4   4     3   2   1   4     7   5   8   8     7   6   5   8
        0   0   0   0   0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0
        0   0   0   0   0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0
        0   0   0   0   0   0   0   0     0   0   0   0     2   2   2   2     0   0   0   0     1   1   1   1   ];

    dove3 = dl([2:5, 10],:) == 3 | dl([2:5, 10],:) == -3;
    dove = [false(1,24); dove3(1:4,:); false(4,24); dove3(5,:)];
    dl(dove) = dl(dove)/3 * DATA.lato_q;

    dove2 = dl([2:5, 10],:) == 3 | dl([2:5, 10],:) == -3  & ~dove3;
    dove = [false(1,24); dove2(1:4,:); false(4,24); dove2(5,:)];
    dl(dove) = dl(dove)/2 * DATA.raggio_mant;
    
    dove1 = dl([2:5, 10],:) == 3 | dl([2:5, 10],:) == -3  & ~dove3 & ~dove2;
    dove = [false(1,24); dove1(1:4,:); false(4,24); dove1(5,:)];
    dl(dove) = dl(dove) * DATA.raggio_osta;
    
    case 2      %  Dominio circolare nuovo
%       1   2   3   4     5   6   7   8     9  10  11  12    13  14  15  16    17  18  19  20
dl =  [ 1   1   1   1     2   2   2   2     1   1   1   1     2   2   2   2     1   1   1   1
        3   0  -3   0     0  -3   0   2     2   0  -2   0     0  -2   0   1     1   0  -1   0
        0  -3   0   3     0  -2   0   3     0  -2   0   2     0  -1   0   2     0  -1   0   1
        0   3   0  -3     2   0  -3   0     0   2   0  -2     1   0  -2   0     0   1   0  -1
        3   0  -3   0     3   0  -2   0     2   0  -2   0     2   0  -1   0     1   0  -1   0
        3   2   1   4     2   2   1   3     7   6   5   8     6   6   5   7     0   0   0   0
        0   0   0   0     3   1   4   4     3   2   1   4     7   5   8   8     7   6   5   8
        0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0
        0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0     0   0   0   0
        3   3   3   3     0   0   0   0     2   2   2   2     0   0   0   0     1   1   1   1   ];

    dove3 = dl([2:5, 10],:) == 3 | dl([2:5, 10],:) == -3;
    dove = [false(1,20); dove3(1:4,:); false(4,20); dove3(5,:)];
    dl(dove) = dl(dove)/3 * DATA.lato_q;

    dove2 = dl([2:5, 10],:) == 3 | dl([2:5, 10],:) == -3  & ~dove3;
    dove = [false(1,20); dove2(1:4,:); false(4,20); dove2(5,:)];
    dl(dove) = dl(dove)/2 * DATA.raggio_mant;
    
    dove1 = dl([2:5, 10],:) == 3 | dl([2:5, 10],:) == -3  & ~dove3 & ~dove2;
    dove = [false(1,20); dove1(1:4,:); false(4,20); dove1(5,:)];
    dl(dove) = dl(dove) * DATA.raggio_osta;
    
    case 3  % Dominio quadrato
dl =  [ 2   2   2   2       1   1   1   1       1   1   1   1
        1  -1  -1   1      -.5  0   .5  0      -.25 0   .25 0
       -1  -1   1   1       0   .5  0  -.5      0   .25 0  -.25
        1   1  -1  -1       0  -.5  0   .5      0  -.25 0   .25
        1  -1  -1   1      -.5  0   .5  0      -.25 0   .25 0
        1   1   1   1       2   2   2   2       0   0   0   0
        0   0   0   0       1   1   1   1       2   2   2   2
        0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       .5  .5  .5  .5      .25 .25 .25 .25];
    
    dl(2:5,1:4) = DATA.lato_q/2*dl(2:5,1:4);
    dl([2:5 end],5:8) = dl([2:5 end],5:8) * 2 * DATA.raggio_mant;
    dl([2:5 end],9:12) = dl([2:5 end],9:12) * 4 * DATA.raggio_osta;

    case {4,10}  %  Dominio circolare
dl =  [ 1   1   1   1       1   1   1   1       1   1   1   1
       -1   0   1   0      -.5  0   .5  0      -.25 0   .25 0
        0   1   0  -1       0   .5  0  -.5      0   .25 0  -.25
        0  -1   0   1       0  -.5  0   .5      0  -.25 0   .25
       -1   0   1   0      -.5  0   .5  0      -.25 0   .25 0
        1   1   1   1       2   2   2   2       0   0   0   0
        0   0   0   0       1   1   1   1       2   2   2   2
        0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       0   0   0   0       0   0   0   0
        1   1   1   1       .5  .5  .5  .5      .25 .25 .25 .25];
    
    dl([2:5 end],1:4) = dl([2:5 end],1:4) * DATA.raggio_ext;
    dl([2:5 end],5:8) = dl([2:5 end],5:8) * 2 * DATA.raggio_mant;
    dl([2:5 end],9:12) = dl([2:5 end],9:12) * 4 * DATA.raggio_osta;
    interno_mantello = false(size(dl,2),1);
    interno_mantello(9:12) = true;
    
        case {5,11}  %  Barca
L_barca = DATA.raggio_osta;
L_cloak = DATA.raggio_mant;
delta = (L_cloak - L_barca)/2;
dl =  [ 1   1   1   1       2    2    1        1        2          2    2    1      1     2
       -1   0   1   0      -.5  -.5   0        1.25^.5  0         -1   -1    0      3^.5  0
        0   1   0  -1      -.5   0    1.25^.5  0       -.5        -1    0    3^.5   0    -1
        0  -1   0   1       .5  -.5  -.5       0        .5         1   -1   -1      0     1
       -1   0   1   0      -.5  -.5   0        .5       .5        -1   -1    0      1     1
        1   1   1   1       0    0    0        0        0          2    2    2      2     2 
        0   0   0   0       2    2    2        2        2          1    1    1      1     1
        0   0   0   0       0    0    0        0        0          0    0    0      0     0
        0   0   0   0       0    0    1       -1        0          0    0    1     -1     0
        1   1   1   1       0    0   1.5      1.5       0          0    0    2      2     0  ];
    
    dl([2:5 end],1:4) = dl([2:5 end],1:4) *  DATA.raggio_ext;
    dl([2:5 end-2:end],5:9) = dl([2:5 end-2:end],5:9) * L_barca;
    dl([2:5 end-2:end],[10 11 14]) = dl([2:5 end-2:end],[5 6 9]) * L_cloak/L_barca;
    dl(8:9,12:13) = dl(8:9,7:8);
    dl(end,[12 13]) =  1.5*L_barca + delta;
    dl(3,12) = sqrt((1.5*L_barca+delta)^2 - L_barca^2);
    dl(2,13) = sqrt((1.5*L_barca+delta)^2 - L_barca^2);
    dl(4,12) = -L_barca/2 -delta;
    dl(5,13) = L_barca/2 + delta;
    
        case 7  %  Intensificatore
dl =  [ 1   1   1   1       1   1   1   1       1   1   1   1       1   1   1   1
       -1   0   1   0      -1   0   1   0      -1   0   1   0      -1   0   1   0
        0   1   0  -1       0   1   0  -1       0   1   0  -1       0   1   0  -1
        0  -1   0   1       0  -1   0   1       0  -1   0   1       0  -1   0   1
       -1   0   1   0      -1   0   1   0      -1   0   1   0      -1   0   1   0
        3   3   3   3       2   2   2   2       4   4   4   4       1   1   1   1
        0   0   0   0       3   3   3   3       2   2   2   2       4   4   4   4
        0   0   0   0       0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       0   0   0   0       0   0   0   0       0   0   0   0
        1   1   1   1       1   1   1   1       1   1   1   1       1   1   1   1   ];
    
    dl([2:5 end],1:4) = dl([2:5 end],1:4) * DATA.raggio_ext;
    dl([2:5 end],5:8) = dl([2:5 end],5:8) * DATA.raggio_mant;
    dl([2:5 end],9:12) = dl([2:5 end],9:12) * DATA.raggio_osta;
    dl([2:5 end],13:16) = dl([2:5 end],13:16) * DATA.raggio_fuoco;
    interno_mantello = false(size(dl,2),1);
    interno_mantello(9:12) = true;
    case 9  %  Intensificatore microfono
        L = DATA.raggio_ext;
        R = DATA.raggio_l_ext;
        r = DATA.raggio_l_in;
        s = DATA.raggio_mic;
        h = DATA.d_mic/2;
        p = @(ii) sqrt(ii^2 - h^2);
%                    External cfr and hidrophone                    cfr ext lens    cfr int lens
dl =  [ 1      2     2     2     1     1     2     2     2     1        1    1        1    1
       -L    p(L)  p(R)  p(r)   -s   p(s)   p(s)  p(r)  p(R)  p(L)      -R  p(R)      -r  p(r)
      p(L)   p(R)  p(r)  p(s)  p(s)   -s    p(r)  p(R)  p(L)   -L      p(R)  -R      p(r)  -r
        0     -h    -h    -h     0     h     h     h     h     h        0    h        0    h
       -h     -h    -h    -h    -h     0     h     h     h     0       -h    0       -h    0
        3      3     2     1     0     0     1     2     3     3        2    2        1    1
        0      0     0     0     1     1     0     0     0     0        3    3        2    2
        0      0     0     0     0     0     0     0     0     0        0    0        0    0
        0      0     0     0     0     0     0     0     0     0        0    0        0    0
        L      0     0     0     s     s     0     0     0     L        R    R        r    r   ];
    
    interno_mantello = false(size(dl,2),1);
    interno_mantello([11 3 13 14 8 12]) = true;
    case 8  %  Dominio circolare - Silent cloak
dl =  [ 1   1   1   1       1   1   1   1       1   1   1   1
       -1   0   1   0      -.5  0   .5  0      -.25 0   .25 0
        0   1   0  -1       0   .5  0  -.5      0   .25 0  -.25
        0  -1   0   1       0  -.5  0   .5      0  -.25 0   .25
       -1   0   1   0      -.5  0   .5  0      -.25 0   .25 0
        1   1   1   1       2   2   2   2       3   3   3   3
        0   0   0   0       1   1   1   1       2   2   2   2
        0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       0   0   0   0       0   0   0   0
        1   1   1   1       .5  .5  .5  .5      .25 .25 .25 .25];
    
    dl([2:5 end],1:4) = dl([2:5 end],1:4) *  DATA.lato_q;
    dl([2:5 end],5:8) = dl([2:5 end],5:8) * 2 * DATA.raggio_mant;
    dl([2:5 end],9:12) = dl([2:5 end],9:12) * 4 * DATA.raggio_osta;
    interno_mantello = false(size(dl,2),1);
    interno_mantello(9:12) = true;
    case {12, 13} % Meshlente_3D
        return
end

% figure

pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
[vertices, boundaries, elements] = initmesh(dl,'Jiggle','minimum','Hgrad',1.05,'Hmax',DATA.hh_mesh);
[vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,find(DATA.tipo_dominio == 'c',1));
[vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,find(DATA.tipo_dominio == 'c',1));
switch DATA.tipo_mesh
    case 7
        [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,1);
%     case 9
%         [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,)
end
% [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,2);
% if any(boundaries(5,:) == 3)
%     [vertices, boundaries, elements] = refinemesh(dl,vertices, boundaries, elements,1);
% end



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
