function [press,dpress_dr,r,theta] = sol_ana(vertices,k_amb,r_obs)
n_vertices = size(vertices,2);
press = zeros(n_vertices,1);
dpress_dr = press;
vertices = vertices' * [1;1i];
r = abs(vertices);
theta = angle(vertices);

jacobi = 2*ones(1e3,1);
jacobi(1) = 1;
ka = k_amb*r_obs;

for nu = 0:50
    num = 1/2 * (besselj(nu-1,ka) - besselj(nu+1,ka));
    den = 1/2 * (besselh(nu-1,2,ka) - besselh(nu+1,2,ka));
    B_nu = - jacobi(nu+1) * (-1i)^nu * num/den;
    H = besselh(nu,2,k_amb*r);
    dH =  1/2 * (besselh(nu-1,2,k_amb*r) - besselh(nu+1,2,k_amb*r)) * k_amb;
    
    press     = press     + B_nu*H .* cos(nu*theta);
    dpress_dr = dpress_dr + B_nu*dH.* cos(nu*theta);
end
