function [A_costo,N_costo] = calcola_vicini(DATA)

Ax = bsxfun(@ minus, DATA.nodi_basi(1,:)', DATA.nodi_basi(1,:));
Ay = bsxfun(@ minus, DATA.nodi_basi(2,:)', DATA.nodi_basi(2,:));
switch DATA.tipo_base
    case 'hex'
        Ax = Ax.^2;
        Ay = Ay.^2;
        A_costo = (Ax + Ay) < (2*DATA.lato_hex)^2;
    case 'corone circolari'
        A_costo = abs(Ax) <= (DATA.dr + 3*eps);
        A_costo = A_costo & (abs(Ay) <= (DATA.dtheta + 3*eps) );
    case 'rettangoli'
        A_costo = abs(Ax) <= (DATA.lato_x + 3*eps);
        A_costo = A_costo & (abs(Ay) <= (DATA.lato_y + 3*eps) );
end

A_costo = A_costo & ~eye(DATA.n_basi);
N_costo = sum(A_costo,2);


end