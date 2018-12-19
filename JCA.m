function mat=JCA(f,mat)


omega = 2*pi*f;

rho_0 = mat.rho;        % [Kg.m-3] density at rest of air at 22C, 1atm
c_0   = mat.c;           % [m.s-1] speed of sound in air at 22C, 1atm
P_0   = 1.0132e+05;      % [N.m-2] atmospheric pressure at 22C, 1atm
a00   = mat.a00;        %tortuosity
p     = mat.p;
L     = mat.L;
eta   = 18.2e-6;%mat.eta;
sigma = mat.sigma ;      % [N.s.m-4] static air flow resistivity of material

% Rockwool Flexi A Plate,
% Density: 30 kg/m3, Air flow resistivity: 9 kPa*s/m2

% Compute variable X and print frequency
% limits of validity for the two models 
mat.rho_s = rho_0*a00./p*(1+ sigma*p./(1i*omega*rho_0*a00)*(1+ (1i*4*a00^2*eta*rho_0)/(p^2*L^2*sigma^2))^(1/2));

g = 1.4;
k = 0.0262;
Cp = 1006;
Pr = Cp*eta/k;
Lt = mat.Lt;

mat.K_s = P_0*g/p/(g-(g-1)/(1+8*eta/(1i*omega*rho_0*Pr*Lt^2)*sqrt(1+Lt^2/16*1i*omega*rho_0*Pr/eta)));
mat.cs = sqrt(mat.K_s/mat.rho_s);