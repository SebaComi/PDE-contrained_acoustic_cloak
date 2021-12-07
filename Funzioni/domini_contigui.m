function coppie = domini_contigui(boundaries)
% coppie(ii,jj) = 1 sse il dominio ii confina col dominio jj
% coppie_new(:,cc) contiene il dominio destro e sinistro del contorno cc
% n_facce = max(max(boundaries(6:7,:)));
% coppie = false(n_facce);
% 
% for ii = 1:n_facce
%     aa = boundaries(6,:) == ii;
%     bb = boundaries(7,:) == ii;
%     
%     cc = [boundaries(6,bb)  boundaries(7,aa)];
%     cc = unique(cc);
%     cc(cc == 0) = [];
%     coppie(cc,ii) = 1;
% end
[~,index] = unique(boundaries(5,:));
coppie = boundaries(6:7,index);
