function [] = disegna_u(U,DATA,lato_disegno_10_standard)
if nargin < 3
    lato_disegno_10_standard = 10;
end
% lato_disegno = 30;
% lato = DATA.lato;
% [aa,bb] = meshgrid([0 1],[0 1]);
% hold on
% for uu = 1:DATA.n_basi
%     x = lato*aa + DATA.nodi_basi(1,uu);
%     y = lato*bb + DATA.nodi_basi(2,uu);
%     surf(x,y,U(uu)*ones(2),U(uu));
% end

colore = jet;
z = (U-min(U)+5*eps);
z = ceil(256*z/max(z));
z = colore(z,:);
switch DATA.tipo_base
    case 'rettangoli'
        x = DATA.nodi_basi(1,:) + DATA.lato_x/2;
        y = DATA.nodi_basi(2,:) + DATA.lato_y/2;
        scatter(x,y,lato_disegno_10_standard,z,'filled','s')
    case 'corone circolari'
        vect = (DATA.nodi_basi(1,:) + DATA.dr).* exp(1i * (DATA.nodi_basi(2,:) + DATA.dtheta));
        if DATA.dtheta > pi/3
            vect = vect.' * exp(1i * DATA.dtheta * (0:149)/150);
            z = repmat(z,150,1);
        end
        x = real(vect(:));
        y = imag(vect(:));
        scatter(x,y,lato_disegno_10_standard,z,'filled','o')
    case 'hex'
        x = DATA.nodi_basi(1,:);
        y = DATA.nodi_basi(2,:);
        scatter(x,y,lato_disegno_10_standard,z,'filled','o')
end

axis equal
if max(U) - min(U) < 1e-8
    caxis([min(U)-1, max(U)+1])
else
    caxis([min(U) max(U)])
end
