function [] = disegna_u_bello(U,DATA,~)
colore = jet;
if max(U) - min(U) < 1e-3
    z = 0.5 * ones(size(U));
    z = ceil(256*z/max(z) - 128);
else
    z = (U-min(U)+5*eps);
    z = ceil(256*z/max(z));
end

z = colore(z,:);
x = DATA.nodi_basi(1,:);
y = DATA.nodi_basi(2,:);

switch DATA.tipo_base
    case 'rettangoli'
        xy = x + 1i*y + [0;   DATA.lato_x;     DATA.lato_x+1i*DATA.lato_y;    1i*DATA.lato_y];
    case 'corone circolari'
        xy = x + 1i*y + [0;   DATA.dr;     DATA.dr+1i*DATA.dtheta;    1i*DATA.dtheta];
        xy = real(xy) .* exp(1i * imag(xy));
    case 'hex'
        xy = x + 1i*y + DATA.lato_hex * exp(1i*pi/3 * (0:5)');
end

hold on
poligoni{DATA.n_basi} = []; ehi = poligoni;
for ii = 1:DATA.n_basi
    poligoni{ii} = polyshape(real(xy(:,ii)),imag(xy(:,ii)));
    ehi{ii} = plot(poligoni{ii});
    ehi{ii}.FaceColor = z(ii,:);
    ehi{ii}.FaceAlpha = 1;
%             ehi{ii}.EdgeAlpha = 0;
end

axis equal
if max(U) - min(U) < 1e-3
    caxis([min(U)-1, max(U)+1])
else
    caxis([min(U) max(U)])
end
