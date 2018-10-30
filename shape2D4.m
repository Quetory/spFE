function [N,DN] = shape2D4(XI)

N = zeros(1,4);
DN = zeros(2,4);

s = XI(1);
t = XI(2);

N(1) = 1/4*(1-s)*(1-t);
N(2) = 1/4*(1+s)*(1-t);
N(3) = 1/4*(1+s)*(1+t);
N(4) = 1/4*(1-s)*(1+t); 

DN(1,1) = -1/4*(1-t);
DN(1,2) =  1/4*(1-t);
DN(1,3) =  1/4*(1+t);
DN(1,4) = -1/4*(1+t);

DN(2,1) = -1/4*(1-s);
DN(2,2) = -1/4*(1+s);
DN(2,3) =  1/4*(1+s);
DN(2,4) =  1/4*(1-s);




