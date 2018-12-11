function  mat = APM(f,mat)

omega = 2*pi*f;

rho_0 = mat.rho;        % [Kg.m-3] density at rest of air at 22C, 1atm
c_0   = mat.c;           % [m.s-1] speed of sound in air at 22C, 1atm
P_0   = 1.0132e+05;      % [N.m-2] atmospheric pressure at 22C, 1atm
sigma = mat.sigma ;      % [N.s.m-4] static air flow resistivity of material

% Rockwool Flexi A Plate,
% Density: 30 kg/m3, Air flow resistivity: 9 kPa*s/m2

% Compute variable X and print frequency
% limits of validity for the two models 

X = f/sigma;
f_min = 0.01*sigma;
f_max = 1.00*sigma;

% Revised expressions of Delany and Bazley model by Miki 
% (NB: gamma = alpha + j beta = j k )

Z_DB70_Mik90 = rho_0*c_0*( 1 + 5.50*(X*1000).^(-0.632) ...
                            - 1i*8.43*(X*1000).^(-0.632) ); 

k_DB70_Mik90 = omega/c_0 .* (-1i) .* ( 11.41*(X*1000).^(-0.618) ...
                                      + 1i* (1 + 7.81*(X*1000).^(-0.618) ) );
                                  
k_DB70_Mik90 = omega/c_0 .* ( 11.41*(X*1000).^(-0.618) ...
                                      + 1i* (1 + 7.81*(X*1000).^(-0.618) ) );
                                  
mat.K_s   = Z_DB70_Mik90.*omega./k_DB70_Mik90;
mat.rho_s = k_DB70_Mik90.*Z_DB70_Mik90./omega;
mat.cs  = sqrt(mat.K_s./mat.rho_s);

assignin('base','mat',mat)
%%%%%
%%%%% Compute sound absorption using the two models 
%%%%% for a sample of thickness d backed by a rigid 
%%%%% and impervious wall under at room temperature
%%%%% and pressure conditions 
%%%%%
% 
% Z = -1i.*Z_DB70./tan(k_DB70*h);
% alpha_DB70 = 1 - ( abs( (Z-rho_0*c_0)./(Z+rho_0*c_0) ) ).^2;
% 
% 
% Z = -1i.*Z_DB70_Mik90./tan(k_DB70_Mik90*h);
% alpha_DB70_Mik90 = 1 - ( abs( (Z-rho_0*c_0)./(Z+rho_0*c_0) ) ).^2;


end

