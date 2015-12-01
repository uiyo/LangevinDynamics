function var = Diffusivity(x,simPara)
%calculating configuration dependent diffusivity coefficient
%author: Yuguang Yang yyang60@jhu.edu

% parallel diffusivity
h = (x(3)-simPara.radius)/x(3);
if h < 0
    h = 1e-6;
end
f_brenner = (6*h^2+2*h)/(6*h^2+9*h+2);
if f_brenner < 0
    error('Error. f_brenner must be positive, f_brenner: %f.',f_brenner)
end
% perpendicular diffusivity
lambda = simPara.radius/x(3);
f_z = 1 - 1.0/16*lambda^3 + 1.0/8.0*lambda^5 - 31.0/256.0*lambda^6;
if f_z < 0
    error('Error. f_z must be positive, fz: %f.',f_z)
end
var = diag([f_brenner,f_brenner,f_z]);
    