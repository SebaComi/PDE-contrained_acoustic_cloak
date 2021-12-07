function [idx_bordo_mant, idx_bordo_osta] = Calcola_cont_mant_e_ost(boundaries,nodi_boundaries,tipo_contorno,normali,tipo_dominio)
for tipo = [3,2]
    contorno_ost = find(tipo_contorno == tipo);
    code = [];     dove = [];
    ii = 1;
    for dd = find(tipo_dominio == 2)
        cc_del_dd = unique(boundaries(5,nodi_boundaries(:,dd)));
        for cc = contorno_ost(:)'
            if any(cc_del_dd == cc)
                vect = normali{dd,cc}(1,:); % questi indici contenuti in "normali" sono già ordinati secondo l'ascissa curvilinea
                code = [code,vect([1,end])];
                dove = [dove,vect];
                lunghe(ii) = length(vect);
                ii = ii + 1;
            end
        end
    end
    if ~isempty(dove)
        idx = dove(1);
        vertici = [];
        inizi = [1, cumsum(lunghe(1:end-1))+1];
        inizi = dove(inizi);
        for ii = 1:length(code)-1
            if any(idx == inizi)
                chi = find(idx == code(1:2:end));
                a = find(idx == dove);
                b = a + lunghe(chi)-1;
                vertici = [vertici,dove(a+1:b)];
            else
                b = a - lunghe(chi)+1;
                vertici = [vertici,dove(a-1:-1:b)];
            end
            idx = vertici(end);
            dove(a:b) = 0;
        end
    else
        vertici = ones(1,3)*idx_bordo_mant(1);
    end
        
        switch tipo
            case 2
                idx_bordo_osta = vertici;
            case 3
                idx_bordo_mant = vertici;
        end
end

