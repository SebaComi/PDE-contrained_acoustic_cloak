function [z] = rettangolo(xy,n,DATA) % ,idx_bordo_mant,idx_bordo_osta
% base 'n' per il controllo fatta a rettangoli di lati = 'lato_x' e 'lato_y' e con
% vertice in basso a sinistra  in (x;y) = nodi_basi(1:2,n)
x = xy{1}; y = xy{2};

z =  ones(size(x));
nodi = DATA.nodi_basi;
z = z .*(x>=nodi(1,n) & x<nodi(1,n)+DATA.lato_x & y>=nodi(2,n) & y<nodi(2,n)+DATA.lato_y);
% NB: il segno = solo da una parte evita che un punto appartenga a due rettangoli


% Questo è inutile perché tutti i quadretti sono completamente all'interno
% del mantello, può tornare utili in caso i quadretti si definissero in
% modo diverso

% xx_quadret = [nodi(1,n), nodi(1,n),      nodi(1,n)+lato,    nodi(1,n)+lato];
% yy_quadret = [nodi(2,n), nodi(2,n)+lato, nodi(2,n)+lato,    nodi(2,n)];
% if ~all(isdentro(xx_quadret,yy_quadret,idx_bordo_mant,idx_bordo_osta,vertices))
%     in = isdentro(x,y,idx_bordo_mant,idx_bordo_osta,vertices);   % se è fuori dal mantello z = 0
%     z = z.*in;
% end


end