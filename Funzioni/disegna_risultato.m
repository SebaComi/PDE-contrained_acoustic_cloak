function [] = disegna_risultato(MESH,Y,limite)
if ~exist('limite','var')
    limite = false;
end
if ~isreal(Y)
    Y = real(Y);
end

switch MESH.dim
    case 2
        pdeplot(MESH.vertices,MESH.elements(1:3,:),'xydata',Y(1:MESH.numVertices),'xystyle','interp',...
                'zdata',Y(1:MESH.numVertices),'zstyle','continuous','colorbar','on','mesh','off','ColorMap',jet);
        view(2)
        xlabel('x'); ylabel('y'); title('Real(Y)')
    case 3
        pdeplot3D(MESH.vertices,MESH.elements(1:4,:),'ColorMapData',Y(1:MESH.numVertices),'Mesh','on','FaceAlpha',1)
        xlabel('x'); ylabel('y'); zlabel('z'); title('Real(Y)')
end
%     axis equal
%     colormap(hot)
if limite
    caxis([-1 1])
end
end