function DATA = Dati_cloak(DATA)
if ~isfield(DATA,'dim')
    DATA.dim = 2;
end
if ~isfield(DATA,'fem')
    DATA.fem = 'P2';
end
DATA.quad_order =  5;

%% Parametri fisici
DATA.rho_0 = 998;
DATA.B_0 = 2.2e9;
DATA.c_0 = sqrt(DATA.B_0/DATA.rho_0);
DATA.n_frq = 1;
% DATA.omega = [9440  11520  13600]*2*pi;
% DATA.omega = [9440  11520  9440 9440]*2*pi;
% DATA.omega = repmat( linspace(9440,13600,10)*2*pi ,1,3 );
DATA.omega = ([-1 1]+13)*1000 *2*pi;
DATA.omega = [13 14]*2000*pi;
DATA.k_amb = DATA.omega / DATA.c_0;

% DATA.direz = [zeros(DATA.n_frq/3,1) , ones(DATA.n_frq/3,1)*pi/4, ones(DATA.n_frq/3,1)*(-pi/4)];
DATA.direz = [0;  pi/2];

%% Tipo di mesh
% tipi di contorno
% - e: Contorno esterno assorbente
% - r: Contorno interno rigido
% - s: Contorno soft, Dirichlet omogeneo
% - c: Contorno del mantello: nessuna condizione al contorno
% - x: Contorno interno auxiliario: nessuna condizione al contorno
% - z: Contorno con un'impedenza precisa da specificare in DATA.z

% Tipi di dominio
% - a: Fluido ambiente imperturbato
% - c: Dominio su cui agisce il controllo
% - b: Altro fluido imperturbato
switch DATA.tipo_mesh
    case 1
        tipo_contorno(1:8)           = 'e';   % 1 = Gamma esterno
        tipo_contorno(21:24)         = 'r';   % 2 = Gamma interno
        tipo_contorno(13:16)         = 'c';   % 3 = Gamma mantello
        tipo_contorno([9:12, 17:20]) = 'x';   % 4 = Gamma ausiliari: Si comportano allo stesso modo dell'interfaccia fluido-mantello

        tipo_dominio(1:4) = 'a';    % ambiente
        tipo_dominio(5:8) = 'c';    % cloak
    case 2
        tipo_contorno(1:4)           = 'e';     % 1 = Gamma esterno
        tipo_contorno(17:20)         = 'r';     % 2 = Gamma interno rigido
        tipo_contorno(9:12)          = 'c';     % 3 = Gamma mantello
        tipo_contorno([5:8,13:16])   = 'x';    % 4 = Gamma ausiliari: Si comportano allo stesso modo dell'interfaccia fluido-mantello

        tipo_dominio(1:4) = 'a';    % fluido
        tipo_dominio(5:8) = 'c';    % cloak
    case {3,4} % Cloak rigido
        tipo_contorno(1:4)  = 'e';     % e = Gamma esterno (assorbente)
        tipo_contorno(9:12) = 'r';     % r = Gamma interno (rigido)
        tipo_contorno(5:8)  = 'c';     % c = Gamma cloak

        tipo_dominio(1) = 'a';    % fluido
        tipo_dominio(2) = 'c';    % rivestimento
    case 5  % barca rigida
        tipo_contorno(1:4)  = 'e';     % 1 = Gamma esterno
        tipo_contorno(5:9)  = 'r';     % 2 = Gamma interno rigido
        tipo_contorno(10:14)  = 'c';     % 3 = Gamma mantello

        tipo_dominio(1) = 'a';    % fluido
        tipo_dominio(2) = 'c';    % mantello
    case 6  % test controllo
        tipo_contorno([1:3 5 7 8 10:13]) = 'e';
%         tipo_contorno([]) = 2;
        tipo_contorno([4 9]) = 'c';
        
        tipo_dominio([1 3]) = 'a';
        tipo_dominio(2) = 'c';
    case 6.1  % test controllo
        tipo_contorno([1:2:7 8 10]) = 'r';
        tipo_contorno([4 6]) = 'c';
        tipo_contorno(9) = 'e';
        tipo_contorno(2) = 'e';
        
        tipo_dominio(1) = 'a';
        tipo_dominio(2) = 'c';
    case 7  % intensificatore
        tipo_contorno(1:4) = 'e';
        tipo_contorno(13:16) = 'x';
        tipo_contorno(5:12) = 'c';

        tipo_dominio(1) = 'b';
        tipo_dominio([3 4]) = 'a';
        tipo_dominio(2) = 'c';
    case 8  % Cloak silente
        tipo_contorno(1:4)  = 'e';     % Gamma esterno
        tipo_contorno(9:12) = 'x';     % Gamma ausiliario
        tipo_contorno(5:8)  = 'c';     % Gamma mantello

        tipo_dominio(1) = 'a';    % fluido
        tipo_dominio(2) = 'c';    % rivestimento
        tipo_dominio(3) = 'b';    % altro fluido (cioè la zona silente)
    case 9  % Intensificstore_mic
        DATA.z = 2.5e6;%exp(1i*pi/4);
        tipo_contorno([1 10])    = 'e';
        tipo_contorno([2:4 7:9]) = 'r';
        tipo_contorno(5:6)       = 'z';
        tipo_contorno(11:14)     = 'c';
        
        tipo_dominio = 'bca';
    case 10 % Cloak soft
        tipo_contorno(1:4)  = 'e';     % e = Gamma esterno (assorbente)
        tipo_contorno(9:12) = 's';     % r = Gamma interno (soft)
        tipo_contorno(5:8)  = 'c';     % c = Gamma cloak

        tipo_dominio(1) = 'a';    % fluido
        tipo_dominio(2) = 'c';    % rivestimento
    case 11  % barca soft
        tipo_contorno(1:4)  = 'e';     % 1 = Gamma esterno
        tipo_contorno(5:9)  = 's';     % 2 = Gamma interno rigido
        tipo_contorno(10:14)  = 'c';     % 3 = Gamma mantello

        tipo_dominio(1) = 'a';    % fluido
        tipo_dominio(2) = 'c';    % mantello
    case 12 % Meshlente_3D
        model = DATA.model;
        dove = [mphselectcoords(model,'geom1',DATA.raggio_ext*[0; 0; 1],'boundary') , ...
                mphselectcoords(model,'geom1',DATA.raggio_ext*[0; 0; -1],'boundary')];
        tipo_contorno(dove) = 'e';

        dove = [mphselectcoords(model,'geom1',DATA.raggio_osta*[0; 0; 1],'boundary') , ...
                mphselectcoords(model,'geom1',DATA.raggio_osta*[0; 0; -1],'boundary'), ...
                mphselectcoords(model,'geom1',DATA.raggio_mant*[0; 0; 1],'boundary') , ...
                mphselectcoords(model,'geom1',DATA.raggio_mant*[0; 0; -1],'boundary')];
        tipo_contorno(dove) = 'c';

        dove = [mphselectcoords(model,'geom1',DATA.raggio_fuoco*[0; 0; 1],'boundary') , ...
                mphselectcoords(model,'geom1',DATA.raggio_fuoco*[0; 0; -1],'boundary')];
        tipo_contorno(dove) = 'x';

        tipo_dominio([1 3]) = 'a';
        tipo_dominio(4) = 'b';
        tipo_dominio(2) = 'c';
    case 13 % Cloak concavo
        tipo_contorno = 'eeeeccss';
        tipo_dominio = 'ac';
end
DATA.tipo_contorno = tipo_contorno;
DATA.tipo_dominio = tipo_dominio;

%% Onda incidente
switch DATA.dim
    case 2
        DATA.y_inc_fun = @(xy,k_amb,direz) 1 * exp(-1i*k_amb * (cos(direz)*xy{1} + sin(direz)*xy{2}));  % Onda piana progressiva lungo direzione direz
    case 3
        % Onda che prosegue parallela al piano xy
        DATA.y_inc_fun = @(xyz,k_amb,direz) 1 * exp(-1i*k_amb * (cos(direz)*xyz{1} + sin(direz)*xyz{2}));  % Onda piana progressiva lungo direzione direz
        DATA.y_inc_dz_fun = @(xyz,k_amb,direz) 0*xyz{1};
end

% DATA.dy_inc_fun = @(xy,k_amb,direz) -1i*k_amb * cosdir(direz) * DATA.y_inc_fun(xy,k_amb,direz);  % Questo funziona solo in 2D!!
DATA.dy_inc_dx_fun = @(xy,k_amb,direz) -1i*k_amb * cos(direz) * DATA.y_inc_fun(xy,k_amb,direz);
DATA.dy_inc_dy_fun = @(xy,k_amb,direz) -1i*k_amb * sin(direz) * DATA.y_inc_fun(xy,k_amb,direz);
%% Tipo di controllo
switch DATA.tipo_base
    case 'vertici'

    case 'rettangoli'
        DATA.lato_x = 10e-3;
        DATA.lato_y = 10e-3;
    case 'corone circolari'
        DATA.dtheta = 2*pi/4;
%         DATA.dtheta = 2*pi / floor(2*pi / DATA.dtheta);
        DATA.dr = (DATA.raggio_mant - DATA.raggio_osta)/8;
%         DATA.dr = 0.01;
    case 'hex'
        DATA.lato_hex = 35e-3 / 3;
end

end