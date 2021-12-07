function out_continuo = stimatore_dist(x,y,out_discreto, vertici_mesh)
% x e y sono vettori che contengono le coordinate di interesse

% definiamo il raggio di ricerca cercando i due punti di vertici_mesh che sono
% più vicini fra loro
% distanza = @(xy1,xy2) sqrt((xy1(1) - xy2(1)).^2 + (xy1(2) - xy2(2)).^2);
mat_dist_x = bsxfun(@ minus, vertici_mesh(1,:)', vertici_mesh(1,:));
mat_dist_x = mat_dist_x.^2;
mat_dist_y = bsxfun(@ minus, vertici_mesh(2,:)', vertici_mesh(2,:));
mat_dist_y = mat_dist_y.^2;
mat_dist = mat_dist_x + mat_dist_y;
clearvars mat_dist_x mat_dist_y
mat_dist(mat_dist == 0) = max(max(mat_dist)); % Non voglio considerare le auto-distanze
raggio = mat_dist(find(mat_dist == min(mat_dist),1));
clearvars mat_dist
raggio = sqrt(raggio);      % sqrt( distanza^2 )
raggio = raggio * sqrt(0.5) * 2; % Il raggio attorno a cui guardare è 2 volte la più piccola distanza fra i punti della mesh
% fprintf('Raggio: %f \n',raggio)
tic

out_continuo = zeros(size(x));
for ii = 1:numel(x)
    att_quad = find(vertici_mesh(1,:) < x(ii) + raggio & vertici_mesh(1,:) > x(ii) - raggio & vertici_mesh(2,:) < y(ii) + raggio & vertici_mesh(2,:) > y(ii) - raggio);
    att_quad = att_quad';
%     fprintf('Attorno quadrato: %i   ',att_quad)
%     fprintf('\n')
    distanze = sqrt((vertici_mesh(1,att_quad) - x(ii)).^2 + (vertici_mesh(2,att_quad) - y(ii)).^2);
    distanze = distanze';
%     pesi = raggio - distanze;
%     dove_no = pesi < 1e-6;    % Togliamo raggio-distanze < 0 ed evitiamo errori numerici in caso (raggio-distanze) ~ 0
%     pesi(dove_no) = [];
%     att_quad(dove_no) = [];
    pesi = exp(-10*distanze); % Questa pesatura è da tarare
    out_continuo(ii) = sum(out_discreto(att_quad) .* pesi) / sum(pesi);
end



end