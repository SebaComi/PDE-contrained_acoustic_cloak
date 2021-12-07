function out = trova_dominio(punti,MESH,dominio_n)
% Il punto appartiene al dominio_n (dominio_n può essere un vettore)?
if size(punti,1) ~= MESH.dim
    error('Ciao!')
end
n_punti = size(punti,2);
elem = MESH.elements;
if length(dominio_n) == 1
    chi_serve = elem(end,:) == dominio_n;
else
    chi_serve = any(elem(end,:) == dominio_n(:));
end
elem = elem(:,chi_serve);
% noe = size(elem,2);
vert = single(MESH.vertices);
invjac = single(MESH.invjac(chi_serve,:,:));
out = false(1,n_punti);
switch MESH.dim
    case 2
        tic
        for jj = 1:ceil(n_punti/200)
            range = (1+(jj-1)*200):min(200*jj,n_punti);
            x_x_0 = punti(1,range) - vert(1,elem(1,:))';
            y_y_0 = punti(2,range) - vert(2,elem(1,:))';
    
            s = invjac(:,1,1) .* x_x_0 + invjac(:,2,1) .* y_y_0;
            t = invjac(:,1,2) .* x_x_0 + invjac(:,2,2) .* y_y_0;
            
            dove = all(cat(3,  s<=1, s>= 0, t<=1, t>=0, s+t<=1), 3);
            out(range) = any(dove,1);
        end
        toc
    case 3
        x_x_0 = punti(1,:) - vert(1,elem(1,:))';
        y_y_0 = punti(2,:) - vert(2,elem(1,:))';
        z_z_0 = punti(3,:) - vert(3,elem(1,:))';

        s = invjac(:,1,1) .* x_x_0 + invjac(:,2,1) .* y_y_0 + invjac(:,3,1) .* z_z_0;
        t = invjac(:,1,2) .* x_x_0 + invjac(:,2,2) .* y_y_0 + invjac(:,3,2) .* z_z_0;
        r = invjac(:,1,3) .* x_x_0 + invjac(:,2,3) .* y_y_0 + invjac(:,3,3) .* z_z_0;
        
        dove = all(cat(3,  s<=1, s>= 0, t<=1, t>=0, r<=1, r>=0, s+t+r<=1),3);
end

% FaceId = zeros(length(punti),1);
% FaceId(any(dove,1)) = elem(end,mod(find(dove),noe));


