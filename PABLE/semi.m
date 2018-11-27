function [P1,Q1,P2,Q2,P,Q,X1,X2,Y1,Y2,wjj,zjj,isolve,duality,exdual ] = ...
         semi( j,anorm,T,P1,Q1,P2,Q2,P,Q,ip,pqnorm,X1,X2,Y1,Y2,...
               wjj,zjj,duality,exdual ) 
%
% Purpose
% =======
%
% Monitor the loss of duality and correction is invoked to maintain 
% semi-duality if necessary 
%
% notation: P1 = P_j,   Q1 = Q_j --- Aj, Bj (oldB), Cj (oldC) 
%           P2 = P_j+1, Q2 = Q_j+1 --- Bj1, Cj1 
%
% Arguments
% =========
%
%   j      (input) integer
%          Lanczos step 
%
%   anorm  (input) double precision 
%
%   T      (input) double precision array, dimension( nbmax*maxit, * )   
%
%   P1     (input/output) double precision array, dimension( n, nbold ) 
%          where nbold = ip(j+1) - ip(j)
%
%   Q1     (input/output) double precision array, dimension( n, nbold ) 
%
%   P2     (input/output) double precision array, dimension( n, nb ) 
%
%   Q2     (input/output) double precision array, dimension( n, nb ) 
%
%   P      (input/output) double precision array, dimension( n, ip(j+1)-1 ) 
%
%   Q      (input/output) double precision array, dimension( n, ip(j+1)-1 ) 
%
%   ip     (input) integer array, dimension( j+2 ) 
%
%   X1     (input/output) double precision array, dimesion(nbmax*maxit, nbmax )
%
%   X2     (input/output) double precision array, dimesion(nbmax*maxit, nbmax )
%
%   Y1     (input/output) double precision array, dimesion(nbmax*maxit, nbmax )
%
%   Y2     (input/output) double precision array, dimesion(nbmax*maxit, nbmax ) 
%
%   wjj    (input/output) double precision array, dimesion( nb, nbold )  
%
%   zjj    (input/output) double precision array, dimesion( nbold, nb ) 
%
%   duality(input/output) double precision array, dimension( j+1 )
% 
%   exdual (input/output) double precision array, dimension( j+1 )
%
% ===================================================================
%
% ------------------------ monitor the loss of duality
%
%if j == 1, 
%
%   X2 = P1.' * Q2;       Y2 = P2.' * Q1;
%   X1 = zeros(size(X2)); Y1 = zeros(size(Y2));
%   wjj = Y2;             zjj = X2;
%   duality(j+1) = eps;   
%   exdual(j+1) = eps;
%
%else 

   nbold = ip(j+1) - ip(j);  
   oldB = T( ip(j-1):ip(j)-1 , ip(j):ip(j+1)-1);
   oldC = T( ip(j):ip(j+1)-1 , ip(j-1):ip(j)-1);
   Aj = T( ip(j):ip(j+1)-1, ip(j):ip(j+1)-1);  
   Bj1 = T( ip(j):ip(j+1)-1 , ip(j+1):ip(j+2)-1);
   Cj1 = T( ip(j+1):ip(j+2)-1 , ip(j):ip(j+1)-1);
%
%  locerr simulates the rounding errors in computing the 3-term 
%  recurrences in the Lanczos procedure
%
   locerr = eps*(anorm + norm( T(1:ip(j+1)-1,1:ip(j+1)-1), 1 ) )* ...
                   pqnorm(j+1); 
%
% ------------------------ X-term recurrences 
%
   irow = size(X2,1); 
   X2 = [ X2; zeros(nbold)]; 
   X3 = T(1:ip(j+1)-1,1:ip(j+1)-1)*X2 - X2*Aj - [X1;wjj]*oldB; 
   ki = irow+1; kf = irow+nbold;
   wjj = P2.'*Q1; 
   X3( ki:kf,:) = X3(ki:kf,:) + Bj1*wjj;    
   X3 = X3 + randn(size(X3))*locerr;
   X3 = X3*pinv(Cj1);    % pinv can be replaced by inv, if we observe the 
                         % the structure of Cj1;   
   X1 = X2;
   X2 = X3; 
%
% ------------------------ Y-term recurrences 
% 
   icol = size(Y2,2); 
   Y2 = [Y2,zeros(nbold)];
   Y3 = Y2*T(1:ip(j+1)-1,1:ip(j+1)-1) - Aj*Y2 - oldC*[Y1,zjj]; 
   ki = icol+1; kf = icol+nbold;  
   zjj = P1.'*Q2; 
   Y3(:,ki:kf) = Y3(:,ki:kf) + zjj*Cj1;
   Y3 = Y3 + randn(size(Y3))*locerr;
   Y3 = pinv(Bj1)*Y3;    % pinv can be replaced by inv, if we observe the   
                         % the structure of Bj1;  
   Y1 = Y2;
   Y2 = Y3;
%
% ------------------------ duality(j+1) measures the loss of duality 
%                          at step j+1
%
   duality(j+1)= max( norm(X2,inf), norm(Y2,1) )/...
                     ( sum(pqnorm(1:j))*pqnorm(j+1) ); 
%
% ------------------------ compute the exact numerical duality (for 
%                          testing purposes only), will be deleted in
%                          the production code (it's very expensive).  
%
   exdual(j+1)=max( norm( P.' * Q2, inf )/( norm(P,'fro')*norm(Q2,'fro') ), ...
                    norm( P2.'*Q, 1 )/( norm(Q,'fro')*norm(P2,'fro') ) );
%end 
%
% ------------------------ correction is invoked to maintain semi-duality 
%                          if necessary
%
isolve = 0; % flag to determine whether to compute the eigendecomposition 
            % of tridiag. matrix T. In our implmentation, if correction
            % occurs at step j, we compute the eigendecomp. of T_j 
%
if( duality(j+1) > sqrt(eps) )
%
    fprintf('------------ correction at j = %5i\n',j);
%
    isolve = 1;  % set flag to compute the eigenproblem of tridiag. matrix T
%
%   retroactive TSMGS
%
    for i=1:(j-1)
       ki = ip(i);
       kf = ip(i+1)-1;
       P1 = P1 - P(:,ki:kf) * ( Q(:,ki:kf).' * P1);
       P2 = P2 - P(:,ki:kf) * ( Q(:,ki:kf).' * P2);
       Q1 = Q1 - Q(:,ki:kf) * ( P(:,ki:kf).' * Q1);
       Q2 = Q2 - Q(:,ki:kf) * ( P(:,ki:kf).' * Q2);
    end
    ki = ip(j);
    kf = ip(j+1)-1;
    P(:,ki:kf) = P1;
    Q(:,ki:kf) = Q1;
    P2 = P2 - P1 * ( Q1.' * P2); % Maintain
    Q2 = Q2 - Q1 * ( P1.' * Q2); % local duality again.
%
%   Reset level of duality.
%
    newdual = eps * norm(T,1) /min(svd(Bj1));
    X1 = ones(size(X1)) * (newdual);
    X2 = ones(size(X2)) * (newdual);
    Y1 = ones(size(Y1)) * (newdual);
    Y2 = ones(size(Y2)) * (newdual);
    duality(j+1) = newdual;
    duality(j) = newdual; 
%
%   compute the exact duality after correction (it's very expensive) 
%
    exdual(j)   = max( norm( P( :,ip(1):ip(j)-1 ).'*Q1, inf)/...
                     ( norm(P( :,ip(1):ip(j)-1 ),'fro')*norm(Q1,'fro') ), ...
                       norm( P1.'*Q( :,ip(1):ip(j)-1 ), 1 )/...
                     ( norm(Q( :,ip(1):ip(j)-1 ),'fro')*norm(P1,'fro') ) );
    exdual(j+1)=max( norm( P.'*Q2, inf)/( norm(P,'fro')*norm(Q2,'fro') ), ...
                     norm( P2.'*Q, 1 )/( norm(Q,'fro')*norm(P2,'fro') ) );
%
end
