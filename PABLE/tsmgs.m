function [R,S] = tsmgs(j,nb,R,S,P,Q,ip)
%
% Purpose
% =======
%
% Two-Sided Modified block Gram-Schmidt process:  Dualize R and S 
% in place against the dual matrices P and Q.
%
% Arguments
% =========
%
% j     (input) integer
%
% nb    (input) integer
%
% R     (input/output) double precision array, dimension( n, nb ) 
%
% S     (input/output) double precision array, dimension( n, nb ) 
%
% P     (input) double precision array, dimension( n, ip(j+1)-1 ) 
%
% Q     (input) double precision array, dimension( n, ip(j+1)-1 ) 
%
% ip    (input) integer array, dimension( j+1 ) 
%
% ================================================================
%
for i=1:j
   ki = ip(i);
   kf = ip(i+1)-1;
   R = R - P(:,ki:kf) * ( Q(:,ki:kf).' * R );
   S = S - Q(:,ki:kf) * ( P(:,ki:kf).' * S );
end
