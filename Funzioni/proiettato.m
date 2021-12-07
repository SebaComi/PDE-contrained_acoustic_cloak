function [U_pr,V_pr] = proiettato(U,V,rho_lim,bulk_lim)

rho = exp(V);
bulk = exp(U);

U_pr = U;
V_pr = V;

fuori = ~ inpolygon(rho(:),bulk(:),  rho_lim,bulk_lim);
if any(fuori)
    distanze = (rho(fuori) - rho_lim(:)').^2   + (bulk(fuori) - bulk_lim(:)').^2;
    switch 4    % Il 4 è il solo sicuramente corretto
        case 1
            [d1,pr1] = min(distanze,[],2);
            distanze(:,pr1) = 1e10;
            [d2,pr2] = min(distanze,[],2);

            U_pr(fuori) = log( (bulk_lim(pr1).*d2 + bulk_lim(pr2).*d1) ./ (d1+d2) );
            V_pr(fuori) = log( (rho_lim(pr1) .*d2 + rho_lim(pr2) .*d1) ./ (d1+d2) );

        case 2
%     % Questa è la priezione ortogonale sul lato AB
            [~,prA] = min(distanze,[],2);   % trovo i punti A
            distanze(:,prA) = 1e10;
            [~,prB] = min(distanze,[],2);   % trovo i punti B

            x_AP = rho(fuori)  - rho_lim(prA);
            y_AP = bulk(fuori) - bulk_lim(prA);

            x_AB = rho_lim(prB) - rho_lim(prA);
            y_AB = bulk_lim(prB)- bulk_lim(prA);

            coeff = (x_AP .* x_AB  +  y_AP .* y_AB) ./ (x_AB.^2 + y_AB.^2);
            coeff(coeff < -0.05) = -0.05;
            coeff(coeff >  1.05) = 1.05;
            V_pr(fuori) = log(   rho_lim(prA)  +   coeff .* x_AB  );
            U_pr(fuori) = log(   bulk_lim(prA) +   coeff .* y_AB  );

        case 3
            % Metodo del confronto con la bisettrice
            % trovo il punto A più vicino a P
            [~,prA] = min(distanze,[],2);
            A = [ rho_lim(prA) , bulk_lim(prA) ];

            % guardo il punto prima e quello dopo
            prB1 = mod(prA-1, length(rho_lim)) + 1;
            prB2 = mod(prA+1, length(rho_lim)) + 1;

            % ottengo i vettori relativi ad A
            AP = [rho(fuori), bulk(fuori)] - A;
            AB1 = [rho_lim(prB1) , bulk_lim(prB1)] - A;
            AB2 = [rho_lim(prB2) , bulk_lim(prB2)] - A;

            % calcolo i coseni
            cos1 = sum(AP.*AB1,2) ./ vecnorm(AP,2,2) ./ vecnorm(AB1,2,2);
            cos2 = sum(AP.*AB2,2) ./ vecnorm(AP,2,2) ./ vecnorm(AB1,2,2);

            % capisco a quale vettore è più allineato, cioè il cui cos è maggiore.
            % NB: sulla bisettrice i cos sono uguali
            interni = cos1 > cos2;

            % definisco le cose "giuste"
            AB = AB2;
            AB(interni) = AB1(interni);
            coseno = cos2;
            coseno(interni) = cos1(interni);

            % calcolo AP'
            AP_pr = repmat(coseno .* vecnorm(AP,2,2) ./ vecnorm(AB,2,2) ,1,2) .* AB;
            % assicuro che la proiezione sia dentro il lato AB
            AP_pr(coseno < 0, :) = 0;
            interni = vecnorm(AP_pr,2,2) > vecnorm(AB,2,2);
            AP_pr(interni,:) = AB(interni,:);

            P_pr = A + AP_pr;

            V_pr(fuori) = log( P_pr(:,1) );
            U_pr(fuori) = log( P_pr(:,2) );
        case 4
            distanze = sqrt(distanze);
            F = sum(fuori);
            L = length(rho_lim);
            % Calcolo dei lati
            a = [rho_lim([2:end 1])'; bulk_lim([2:end 1])'] - [rho_lim(:)'; bulk_lim(:)'] ;
            
            b_x = rho(fuori) - rho_lim(:)';
            b_y = bulk(fuori) - bulk_lim(:)';
            
            prod = b_x .* a(1,:) + b_y .* a(2,:);
            
            lunghe2 = sum(a.^2,1);% lunghezze dei lati
            criterio = prod ./ lunghe2 ;
            interni = criterio >=0 & criterio <=1; % In questo modo solo le proiezioni SUI lati sono considerate
            
            % troviamo tutte e altezze facendo il prodotto vettoriale
            h = abs( b_y .* a(1,:) - b_x .* a(2,:) ) ./ sqrt(lunghe2);
            
            h(~interni) = inf;
            [~,dove] = min([h distanze],[],2);
            
            chi_h = find(dove <= L);
            dove_h = dove(chi_h);           % Indici di chi ha come distanza minore le altezze
            chi_dist = find(dove > L);  % Indici di chi ha come dist minore le dist dai vertici
            dove_dist = dove(chi_dist) - L;
            
            c = zeros(2, F);
            if ~isempty(chi_dist)
                c(:,chi_dist) = [rho_lim(dove_dist)'; bulk_lim(dove_dist)'];
            end
            if ~isempty(chi_h)
                c(:,chi_h) = [rho_lim(dove_h)'; bulk_lim(dove_h)'] + a(:,dove_h) .* criterio(chi_h + (dove_h-1)*F)';
            end
            U_pr(fuori) = log( c(2,:) );
            V_pr(fuori) = log( c(1,:) );
            
            
    end
end
