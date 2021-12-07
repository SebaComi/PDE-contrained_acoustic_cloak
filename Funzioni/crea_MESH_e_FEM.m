function [MESH, FE_SPACE,DATA] = crea_MESH_e_FEM(DATA,MESH_0,DATA_0)
    fprintf('Creazione mesh. '); tic
%     if nargin < 2
%         nome = 'mesh_cloak';
%     end
    if isfield(DATA,'mesh')
        eval(['[vertices, boundaries, elements, edges ,DATA] = ' DATA.mesh '(DATA);'])
%             [vertices, boundaries, elements,~,DATA] = mesh_cloak(DATA);
%             [vertices, boundaries, elements,~,DATA] = mesh_WW(DATA);
%             [vertices, boundaries, elements,~,DATA] = Mesh_lente_3D(DATA);
%             [vertices, boundaries, elements,~,DATA] = Mesh_concavo(DATA);

        % Create and fill the MESH data structure
        MESH = buildMESH( DATA.dim, elements, vertices, boundaries, DATA.fem, DATA.quad_order,[],'CFD');
        MESH.edges = edges;
        switch MESH.dim
            case 2
                MESH.bound_tipo_cont = 5;
            case 3
                MESH.bound_tipo_cont = 7;
        end
    %     MESH.coppie = domini_contigui(MESH.boundaries);
    else
        [MESH] = crea_MESH_secondaria(MESH_0,DATA_0,DATA.tipo_dominio, 3);
        switch MESH.dim
            case 2
                MESH.bound_tipo_cont = 5;
            case 3
                MESH.bound_tipo_cont = 7;
        end
    end
    MESH.n_domini = length(DATA.tipo_dominio);
    [MESH.nodi_domini,MESH.nodi_contorni,MESH.nodi_boundaries,MESH.nodi_elements] = distribuzione_domini(MESH,DATA);
%     MESH = rmfield(MESH,'nodi_contorni');
%     Ora uso le normali di redbKIT
%     MESH.normali = versori_normali(MESH.boundaries,MESH.vertices,MESH.nodi_boundaries);
    MESH.n_contorni = length(DATA.tipo_contorno);
    
    
    % Create and fill the FE_SPACE data structure
    [ FE_SPACE ] = buildFESpace( MESH, DATA.fem, 1, DATA.quad_order );
    
    fprintf('           Fatto in %.3f s\n',toc)
end