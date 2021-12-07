% LEGGIMI
% 
% 1. Runna questo file assicurandoti di essere dentro la cartella "Funzioni simili"
% 2. Nella cartella principale c'è il main: NUOVO_cloak.
% 3. Runnalo
% 4. Nella stessa cartella c'è gradiente_discreto
% 
% hint risparmia tempo
% 1. Cambiare lato della MESH:                      NUOVO_cloak, riga 9:    DATA.hh_mesh = 0.10;
% 2. Cambiare il numero di iterazioni massime:      NUOVO_cloak, riga 66:   MaxIt = 50;
% 3. Cambiare il numero di frequenze considerate:   Dati_cloak, riga 7:     DATA.n_frq = 2
% 4. Cambiare direzioni e frequenze:                Dati_cloak, riga 8 e 9
% 5. Cambiare il tipo di controllo: NUOVO_cloak, riga 9: DATA.tipo_base \in {'hex', 'corone circolari', 'rettangoli'}
% 6. Cambiare dimensioni dei "quadretti" del controllo
%   _'rettangoli': Dati_cloak, righe 105 e 106
%   _'corone circolari': Dati_cloak, righe 128 e 131
%   _'hex': Dati_cloak, riga 141


fprintf('Assemblaggio 1/6 \n')
mex -R2018a Assembla_A_e_M.c

fprintf('Assemblaggio 2/6 \n')
mex -R2018a Assembla_C.c

fprintf('Assemblaggio 3/6 \n')
mex -R2018a Assembla_D.c

fprintf('Assemblaggio 4/6 \n')
mex -R2018a Assembla_E.c

fprintf('Assemblaggio 5/6 \n')
mex -R2018a Assembla_solo_si.c

fprintf('Assemblaggio 6/6 \n')
mex -R2018a Assembla_solo_mu.c

fprintf('\n Assemblaggio completato con successo\n')

