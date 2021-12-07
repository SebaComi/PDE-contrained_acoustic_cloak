

figure(4)
hold on
% subplot(231)
colori = get(gca,'ColorOrder');
colori(3,:) = [];
    % Domini e normali
for dd = 1:n_domini

    % Gradienti 
    cc_del_dd = unique(boundaries(5,nodi_boundaries(:,dd)));
    for cc = cc_del_dd(:)'
        if tipo_contorno(cc) == 3
            quiver(vertices(1,nodi_contorni{cc}), vertices(2,nodi_contorni{cc}),...
                dY_cont{dd,cc}, dY_cont{dd,cc}, ...
                'Color',colori(dd,:), 'LineWidth',1,'AutoScale','off')
        end
    end
end

