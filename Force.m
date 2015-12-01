function F = Force(x, simPara)
% calculating the configuration dependent force  
%author: Yuguang Yang yyang60@jhu.edu

    F = zeros(3,1);
%% parallel forces
    F(1) = simPara.electr_mobility * simPara.Efield ...
        * cos(simPara.step*simPara.dt*2*pi*simPara.frequency)* abs(x(1))*0.5/simPara.radius;

    F(1) = F(1) -0.2*simPara.electr_mobility * simPara.Efield * x(1)/simPara.radius;
%     % add twDEP+DEP forces
%     [Rep,Imp,E,phi] = fourier(x,y,d);
%     E1 = Rep;
%     E2 = Imp;
%     % here V is size of the particle
%     F = 1/4*V*Vpp^2*(real(alph)*Rep + 2*imag(alph)*Imp);
%     Fdep = 1/4*V*Vpp^2*real(alph)*Rep;
%     Ftw = 1/4*V*Vpp^2*2*imag(alph)*Imp;


%     % Forces for y >= 1.5*d
%     A1 = 0.877354*Vpp;
%     Fx = -pi^3/64/d^3*V*A1^2*imag(alph)*exp(-pi/2/d*y);
%     Fy = -pi^3/64/d^3*V*A1^2*real(alph)*exp(-pi/2/d*y);
%% perpendicular forces
       
% add gravity
    F(3) = simPara.grav;
% add electrostatic replusion from the botton wall
    F(3) = F(3) + simPara.electr_prefactor * exp(-simPara.kappa*(x(3) - simPara.radius));