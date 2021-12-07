function DATA = Dati_WW(DATA)
if ~isfield(DATA,'dim')
    DATA.dim = 2;
end
if ~isfield(DATA,'fem')
    DATA.fem = 'P2';
end
DATA.quad_order = 5;
addpath('C:\Users\Matteo\Documents\GitHub\tensor_toolbox\')
%% Parametri fisici
DATA.n_frq = 1; %0*3;

DATA.omega = 1*2*pi*2;   % [rad/s]
DATA.h = 3;         % [m]
DATA.g = 9.81;      % [m/s2]
DATA.rho = 998;     % [kg/m3]
DATA.T_0 = 0.07275; % [N/m]
fun = @(k) DATA.g *k *tanh(k *DATA.h) - DATA.omega^2;
options = optimoptions("fsolve",'Display','none');
DATA.k_amb = fsolve(fun,DATA.omega^2 / DATA.g, options);    % [rad/m]
DATA.c = DATA.omega ./ DATA.k_amb;  % [m/s]
DATA.direz = -pi/2; % Nel problema 2D ha senso solo se = 0  % [rad]


%% Tipo di mesh
% tipi di contorno
% - e: Contorno esterno assorbente
% - r: Contorno interno rigido
% - c: Contorno del mantello: nessuna condizione al contorno
% - x: Contorno interno auxiliario: nessuna condizione al contorno
% - z: Contorno con un'impedenza precisa da specificare in DATA.z

% Tipi di dominio
% - a: Fluido ambiente imperturbato
% - c: Dominio su cui agisce il controllo
% - b: Altro fluido imperturbato
switch DATA.tipo_mesh
    case 1
        tipo_contorno([1 3])   = 'e';
        tipo_contorno([5:7]) = 'r';
        tipo_contorno([4 8])   = 'z';
%         DATA.z = -DATA.g / DATA.omega^2;
        tipo_dominio = 'a';
    case 2
        tipo_contorno([1 3]) = 'e';
        tipo_contorno([6:8]) = 'r';
        tipo_contorno([5 9]) = 'c';
        tipo_contorno([4 10]) = 'z';
        % il fondale è lasciato ' ' perché è rigido e non voglio metterlo
        % nei costi
        
        tipo_dominio = 'a';
    case 3  % WW massa 3D
        tipo_contorno([1, 3, 9]) = 'r';
        tipo_contorno( 4 ) = 'z';
        tipo_contorno([2, 8]) = 'e';
        tipo_contorno([5]) = 'c';
        tipo_contorno([6, 7]) = 'g';

        tipo_dominio = 'a';
    case 4  % WW massa e rigidezza 3D
        tipo_contorno([1, 3, 9]) = 'r';
        tipo_contorno( 4 ) = 'f';
        tipo_contorno([2, 8]) = 'e';
        tipo_contorno([5]) = 'c';
        tipo_contorno([6, 7]) = 'g';

        tipo_dominio = 'a';

        tipo_contorno_2([4, 6, 9, 12, 13, 14, 21]) = 'r';
        tipo_contorno_2([5, 18]) = 'e';
end

DATA.tipo_contorno = tipo_contorno;
DATA.tipo_dominio = tipo_dominio;
if exist('tipo_contorno_2','var')
    DATA.tipo_contorno_2 = tipo_contorno_2;
end

%% Onda incidente
costante = -1i*DATA.g*1/DATA.omega/cosh(DATA.k_amb*DATA.h);
switch DATA.dim
    case 2
        DATA.y_inc_fun = @(xy,k,direz)  costante* cosh(k*(xy{2}+DATA.h)) .* exp(-1i*k * xy{1});  % Sulla superficie è un'onda piana progressiva lungo direzione direz
        DATA.dy_inc_dx_fun = @(xy,k,direz) -1i*k * DATA.y_inc_fun(xy,k,direz);
        DATA.dy_inc_dy_fun = @(xy,k,direz)  costante* k * sinh(k*(xy{2}+ DATA.h)) .* exp(-1i*k *xy{1});
        DATA.dy_inc_fun = @(xy,k,direz) [DATA.dy_inc_dx_fun(xy,k,direz);
                                         DATA.dy_inc_dy_fun(xy,k,direz)];
    case 3
        DATA.y_inc_fun = @(xyz,k,direz)  costante* cosh(k*(xyz{3}+DATA.h)) .* exp(-1i*k * (cos(direz)*xyz{1} + sin(direz)*xyz{2}));  % Sulla superficie è un'onda piana progressiva lungo direzione direz
        DATA.dy_inc_dx_fun = @(xyz,k,direz) -1i*k * cos(direz) * DATA.y_inc_fun(xyz,k,direz);
        DATA.dy_inc_dy_fun = @(xyz,k,direz) -1i*k * sin(direz) * DATA.y_inc_fun(xyz,k,direz);
        DATA.dy_inc_dz_fun = @(xyz,k,direz)  costante* k * sinh(k*(xyz{3}+ DATA.h)) .* exp(-1i*k * (cos(direz)*xyz{1} + sin(direz)*xyz{2}));
        DATA.dy_inc_fun = @(xyz,k,direz) [DATA.dy_inc_dx_fun(xyz,k,direz);
                                          DATA.dy_inc_dy_fun(xyz,k,direz);
                                          DATA.dy_inc_dz_fun(xyz,k,direz)];
end



% %% Tipo di controllo
% switch DATA.tipo_base
%     case 'vertici'
% 
%     case 'rettangoli'
%         DATA.lato_x = 10e-3;
%         DATA.lato_y = 10e-3;
%     case 'corone circolari'
%         DATA.dtheta = 2*pi;
% %         DATA.dtheta = 2*pi / floor(2*pi / DATA.dtheta);
% %         DATA.dr = (DATA.raggio_mant - DATA.raggio_osta)/20;
%         DATA.dr = 0.01;
%     case 'hex'
%         DATA.lato_hex = 35e-3 / 3;
% end

end