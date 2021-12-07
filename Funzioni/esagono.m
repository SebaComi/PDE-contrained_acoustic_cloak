function [z] = esagono(xy,n,DATA)
%     Gli esagoni hanno due lati orizzontali !!!
x = xy{1}; y = xy{2};

z = zeros(size(x));
lato = DATA.lato_hex;
latosq3 = lato*sqrt(3);

x_c = DATA.nodi_basi(1,n);
y_c = DATA.nodi_basi(2,n);

x = x - x_c;
y = y - y_c;

in =  (y >= -latosq3/2 & y < latosq3/2   & ...  -lato*sqrt(3)/2 < y < lato*sqrt(3)/2                                 Limiti _ e ‾
       -latosq3 <= (y - sqrt(3)*x)  & (y - sqrt(3)*x) < latosq3  &...  -sqrt(3)*lato < y - sqrt(3)*x < sqrt(3)*lato  Limiti / e /
       -latosq3 <= (y + sqrt(3)*x)  & (y + sqrt(3)*x) < latosq3);   % -sqrt(3)*lato < y + sqrt(3)*x < sqrt(3)*lato   Limiti \ e \
%    NB: il segno = solo da una parte evita che un punto appartenga a due esagoni

z(in) = 1;

