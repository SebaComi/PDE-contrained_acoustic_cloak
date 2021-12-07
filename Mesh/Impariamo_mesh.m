clc
clear
close all


%% Build & refine mesh

% gg = [1;
%       0;
%       0;
%       1];
%   
% dl = decsg(gg);

dl =  [ 2   2   2   2       1   1   1   1       1   1   1   1
        1  -1  -1   1      -.5  0   .5  0      -.25 0   .25 0
       -1  -1   1   1       0   .5  0  -.5      0   .25 0  -.25
        1   1  -1  -1       0  -.5  0   .5      0  -.25 0   .25
        1  -1  -1   1      -.5  0   .5  0      -.25 0   .25 0
        1   1   1   1       2   2   2   2       3   3   3   3
        0   0   0   0       1   1   1   1       2   2   2   2
        0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       0   0   0   0       0   0   0   0
        0   0   0   0       .5  .5  .5  .5      .25 .25 .25 .25];
pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
figure
% return
[vertices, boundaries, elements] = initmesh(dl,'Jiggle','minimum','Hgrad',1.05,'Hmax',0.07);

u = zeros(size(vertices));

 pdeplot(vertices,elements(1:3,:))
 