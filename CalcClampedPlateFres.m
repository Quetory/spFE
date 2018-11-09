function fres = CalcClampedPlateFres(a,b,h,mat)

% a/b ratio
rab = [0.4 2/3 1.0 1.5 2.5];

% for 11
lambda2 = [23.65 27.01 35.99 60.77 147.80];

l2_i = interp1(rab, lambda2,a/b,'pchip');

fres = l2_i/2/pi/a^2*sqrt( mat.E*h^2/12/mat.rho/(1-mat.nu^2));
