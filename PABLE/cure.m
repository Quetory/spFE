function [del,P2,Q2,U,sigma,V,omega,tolbd] = ...
               cure(n,nb,del,j,P2,Q2,P,Q,U,sigma,V,ip,omega,tolbd) 
%
% Purpose 
% =======
%
% Increase the current blocksize to cure breakdown and/or adapt
% blocksize to the order of the largest cluster.
%
% Arguments
% =========
%
% n       (input) integer
%
% nb      (input) integer 
%
% del     (input) integer
%
% j       (input) integer
%
% P2      (input/output) double precision array, dimension( n, nb+del ) 
%
% Q2      (input/output) double precision array, dimension( n, nb+del ) 
%
% P       (input) double precision array, dimension( n, ip(j+1)-1 ) 
%
% Q       (input) double precision array, dimension( n, ip(j+1)-1 ) 
%
% ip      (input) integer array, dimension( j+1 )  
%
% U       (output) double precision, dimension( nb+del,nb+del ) 
%
% sigma   (output) double precision, dimension( nb+del )
%
% V       (output) double precision, dimension( nb+del,nb+del ) 
%
% omega   (input/output) double precision array, dimension( j ) 
%
% tolbd   (input/output) double precision 
%
% ===================================================================
%
%
icount = 1;
P2=[P2 zeros(n,del)];
Q2=[Q2 zeros(n,del)];
minofvu = 0;         
%
%  omega(j+1) == old pivot
%
while ( minofvu < tolbd && icount <= 1 )
%
%  choice 1: the augumented vectors are chosen as random vectors. 
%
   P2(1:n,nb+1:nb+del) = randn(n,del);
   Q2(1:n,nb+1:nb+del) = P2(1:n,nb+1:nb+del);
%
%  Dualize P2(:,nb+1:nb+del) and Q2(:,nb+1:nb+del)
%  against the previous Lanczos vectors via TSMGS (with the choice 3,
%  this step should not be necessary) 
%
   [P2(:,nb+1:nb+del),Q2(:,nb+1:nb+del)] = ...
        tsmgs(j,del,P2(:,nb+1:nb+del),Q2(:,nb+1:nb+del),P,Q,ip ); 
%
%  Orthogonalization inside the block
%
   bet = Q2(:, 1:nb).'*Q2(:, nb+1:nb+del);
   Q2(:, nb+1:nb+del) = Q2(:, nb+1:nb+del) - Q2(:, 1:nb)*bet;
   Q2(:, nb+1:nb+del) = orth(Q2(:, nb+1:nb+del));
%
   bet = P2(:, 1:nb).'*P2(:, nb+1:nb+del);
   P2(:, nb+1:nb+del) = P2(:, nb+1:nb+del) - P2(:, 1:nb)*bet;
   P2(:, nb+1:nb+del) = orth(P2(:, nb+1:nb+del));
%
   [Utemp,sigmatemp,Vtemp] = svd(P2.'*Q2);
   sigmatemp = diag(sigmatemp);
   minofvu = min(sigmatemp);
%
   icount=icount+1;
%
end
%
%  If the above procedure fails to increase omega(j+1), 
%  we take omega(j+1) and continue the procedure without increasing
%  blocksize, and set a smaller tolerance value tolbd. 
%
if minofvu < omega(j+1)
%
%  Reset P2 and Q2 to their original values
%
   P2 = P2(1:n,1:nb);
   Q2 = Q2(1:n,1:nb);
   minofvu = omega(j+1);
   del = 0;
%
else
%
%  Update U, Sigma, V, omega. 
%
   U = Utemp;
   sigma = sigmatemp;
   V = Vtemp;
   omega(j+1) = minofvu; 
%
   clear Utemp; clear sigmatemp; clear Vtemp;
%
end
%
fprintf('after %5i iterations, the new minofvu = %9.4e\n',icount-1,minofvu);
%
%   If (T,Omega) is graded, which could happen oftenly, then blocks will 
%   not increase sigma_min ( Omega ); A smaller sigma_min ( Omega ) is 
%   required.
%
tolbd = min( tolbd, omega(j+1) );
