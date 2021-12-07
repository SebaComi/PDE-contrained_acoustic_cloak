function [vertices, boundaries, elements, dl,DATA] = Mesh_concavo(DATA)

model = DATA.model;

[ciao,meshdata] = mphmeshstats(model);
% info = mphxmeshinfo(model);
vertices = meshdata.vertex;

elements = [ double(meshdata.elem{strcmp(meshdata.types,'tri')}+1);
            meshdata.elementity{strcmp(meshdata.types,'tri')}' ];
elements = double(elements);

boundaries = [ meshdata.elem{strcmp(meshdata.types,'edg')}+1;
               zeros(size(meshdata.elem{strcmp(meshdata.types,'edg')}))
               meshdata.elementity{strcmp(meshdata.types,'edg')}' ];
boundaries = double(boundaries);

dl = [];

DATA.model = model;