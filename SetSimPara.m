%set simulation parameter
%author: Yuguang Yang yyang60@jhu.edu


%% simulator parameters
simPara.nstep = 200000; % number of integration step
simPara.step = 0; % step counter
simPara.dt = 1e-4; % integration time step, unit s
simPara.dim = 3; % dimensionality of the simulation
simPara.x0 = [0;0;1470]*1e-9; % initial configuration, unit m 
simPara.saveInterval = 100; % configration save interval


%% simulaiton physical constants
simPara.radius = 1450e-9; % radius of spherical particle, unit m
simPara.kb = 1.38e-23; % Boltzmann constant
simPara.T = 293.0; %temperature unit K
simPara.vis = 1e-3; % viscosity
simPara.rho_p = 1960; % density of particle (silica), SI unit
simPara.rho_f = 1000; % density of medium, SI unit
simPara.dielectric0 = 8.85e-12; % dieclectric constant of free space SI
simPara.dielectricm = 78; % relative dielectric constant of medium
simPara.mass = 4.0/3.0*pi*(simPara.radius)^3*(simPara.rho_p-simPara.rho_f); % mass
simPara.grav = -simPara.mass*9.8; % gravity force
simPara.D = simPara.kb * simPara.T / (6*pi*simPara.vis*simPara.radius); % stokes diffusivity
simPara.mobility = simPara.D / simPara.kb / simPara.T; % mobility equals D/kT

% parameters for electrostatics
simPara.kappa = 1/10.0 * 1e9; % debye length, unit m^-1
simPara.electr_prefactor = 2.29 * (simPara.kb*simPara.T) * simPara.kappa; % prefactor for elecstratic repulsion


% parameters for electrophoresis
simPara.frequency = 0.5; % Hz
simPara.zeta = 50e-3; % zeta potential difference 50mV
simPara.electr_mobility =  simPara.dielectricm * simPara.dielectric0 *6*pi*simPara.radius*simPara.zeta; % electropheresis mobility
simPara.Efield = 20000; % V/m

%% parameter for function handles
simPara.force = @Force; % deterministic force calculation 
simPara.diffusivity = @Diffusivity; % diffusivity calculation

%% parameters for twDEP
addpath('twDEP')
simPara.vol = 4.0/3.0*pi*(simPara.radius)^3;
simPara.d = 20e-6;  % the gap of the electrodes

% Electrode dimensions
de = 10e-6;                                                          % electrode width [m]
dg = 10e-6;                                                          % gap width [m]
d = (de + dg)/2;
ne = 8;                                                                 % number of electrodes
eo = 0.0000000000088542;                                    % Vacuum permittivity
% Medium parameters
em = 78.9*eo;
sm = 0.0001;                                                        % conductivity [S/m]
mu = 0.001002;                                                     % viscosity (Pa*s)
rho_m = 1000;

% Particle parameters
ep = 3*eo;
sp = 0.001;                                                           % conductivity [S/m]
r = simPara.radius;                                                            % radius [m]
rho_p = 1080;

% Field conditions
Vpp = 1;                                                             % peak-to-peak voltage [V]
Vpp = Vpp/sqrt(8);
fr = 1e6;                                                               % field frequency [Hz]

%% CALCULATING THE twDEP FORCE ON THE PARTICLE AT HEIGHTS FAR FROM THE ELECTRODE SURFACE
epc = ep - sp*1i/fr;
emc = em - sm*1i/fr;

fcm = (epc - emc)/(epc + 2*emc);
alph = 3*em*fcm;
Adep = 0.372923;

V = 4/3*pi*r^3;

% Calculate the electric field
part = 20;
xt = 0:d/part:4*d;
% here we fix y
yt = 0.2*d;
[x,y] = meshgrid(xt,yt);

[Rep,Imp,E,phi] = fourierSum(x,y,d);
E1 = Rep;
E2 = Imp;
F_total = 1/4*V*Vpp^2*(real(alph)*Rep + 2*imag(alph)*Imp);
Fdep = 1/4*V*Vpp^2*real(alph)*Rep;
Ftw = 1/4*V*Vpp^2*2*imag(alph)*Imp;
