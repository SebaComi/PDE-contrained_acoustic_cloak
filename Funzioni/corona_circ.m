function [in] = corona_circ(xy,n,DATA)
dtheta = DATA.dtheta;
dr = DATA.dr;
nodi = DATA.nodi_basi;
x = xy{1};
in = ones(size(x));

switch DATA.dim
    case 2
        y = xy{2};
        vect = x + 1i*y;
        r = abs(vect);
        theta = angle(vect);
        theta = mod(theta,2*pi);
        
        in = in .* (    r >= nodi(1,n)  &      r < nodi(1,n) + dr      & ...
                  theta >= nodi(2,n)  &  theta < nodi(2,n) + dtheta );
    case 3
        y = xy{2}; z = xy{3};
        vect = x + 1i*y;
        r_xy = abs(vect);
        theta1 = angle(vect);
        theta1 = mod(theta1,2*pi);
        vect = r_xy + 1i*z;
        theta2 = angle(vect);
        theta2 = mod(theta2,2*pi);
        r = sqrt(r_xy.^2 + abs(vect).^2);
    
        in = in .* (   r >= nodi(1,n)  &       r < nodi(1,n) + dr      & ...
                  theta1 >= nodi(2,n)  &  theta1 < nodi(2,n) + dtheta  & ...
                  theta2 >= nodi(3,n)  &  theta2 < nodi(3,n) + dtheta  );
end


