%
% test driver routine of PABLE 
%

t1 = cputime; % start CPU clock 

%
% problem depdent parameters 
% ==========================
%
n = size(A,1);  % the order the matrix A
maxmsz = n;     % maximal memory array size for storing Lanczos vectors (P,Q)

if problem == 1, 
   anorm0 = norm( A, 1 ); % the 1-norm of A;
elseif problem == 2, 
   anorm0 = max( norm(A,1), norm(B,1) );
end 
anorm = anorm0; 
fanorm = anorm0;  
nonnormal = norm(A'*A-A*A','fro')/norm(A*A,'fro'); % departure from normality

%
%  Optional pre-processing: matrix balancing
%	
if matbal == 1 & problem == 1, 
   A0 = A; 
   [D,A1] = balance( full(A) ); 
   A1 = sparse( A1 );  % use the matlab function
%   [D,A1]=spbalance(A);  % use the sparse balance function developed by Yi 
   anorm_b = norm( A1, 1 );                         
   nonnormal_b = norm(A1'*A1 - A1*A1','fro')/norm(A1*A1,'fro')
   anorm = anorm_b;  
   fanorm = anorm_b;  
   A = A1; 
   clear A1;
end

%
%  Optional pre-processing: shift-and-invert spectral transformation 
%
if problem == 1, 
   if precond == 1 | precond == 2,
      [L,U] = lu( A - alpha*speye(n) );
      [c,v] = condest( A - alpha*speye(n) ); 
      fanorm = c/anorm;   % norm of transformed matrix 
      if c >= 1/sqrt(eps),
        fprintf('\nTransformed matrix is ill conditioned, cond= %10.5e\n',c)
%        break;
      end 
   else
      L = 0;
      U = 0;
   end 
elseif problem == 2, 
   if precond == 1,
      [L,U] = lu( A - alpha*B );
      [c,v] = condest( A - alpha*B );  
      fanorm = c/anorm;  
   else
      L = 0;
      U = 0;
   end 
end

printin  % print INPUT information

%
% ==========================================================================
%
% call PABLE 
%
% ==========================================================================
 
[j,szoft,T,ritzvalue,VL,VR,P,Q,ip,neig,ritz,ipritz,resl,resr,omega,...
 duality,exdual,tolconv] = ...
   pable(problem,fulldual,semidual,group,treatbd,n,maxmsz,...
         maxeig,nb,maxit,precond,A,B,alpha,L,U,fanorm); 

%
%  Post-processing I: 
%     transfer the Ritz values of f(A) to the Ritz values of A, 
%
if precond == 1
   ritzvalue = alpha + 1./ritzvalue; 
   ritz = alpha + 1./ritz; 
end 

%
%  Post-processing II: 
%     undo matrix balance transformation if necessary

if fulldual == 1 || semidual == 1
   if matbal == 1 && problem == 1  % if initial matrix is balanced
      Q = D*Q; 
      P = D'\P;
      A = A0; 
   end
end 

%
% for converged Ritz values, compute
%   (1) approximate eigenvectors 
%   (2) ``exact'' residual norms
%   (3) estimated condition numbers    
%

if ( (fulldual == 1 || semidual == 1) && neig > 0)   % Level 2 and up 
%
%     compute the left (W) and right (Z) approximate eigenvectors 
%
      for k = 1:neig 
          Z(:,k) = Q(:,1:szoft)*VR(:,ipritz(k)); 
          W(:,k) = conj(P(:,1:szoft))*VL(:,ipritz(k)); 
      end 
      if problem == 2
         W(:,1:neig) = B'\W(:,1:neig);  
      end  
%
%     compute exact residuals 
%
      for k = 1:neig
          if problem == 1
             if matbal == 1
%                Res_left = W(:,k)'*(D*A\D) - ritz(k)*W(:,k)'; 
%                Res_right = (D*A\D)*Z(:,k) - Z(:,k)*ritz(k);  
                Res_left = W(:,k)'*A - ritz(k)*W(:,k)'; 
                Res_right = A*Z(:,k) - Z(:,k)*ritz(k);  
             else 
                Res_left = W(:,k)'*A - ritz(k)*W(:,k)'; 
                Res_right = A*Z(:,k) - Z(:,k)*ritz(k);  
             end 
          else
             Res_left = W(:,k)'*A - ritz(k)*W(:,k)'*B; 
             Res_right = A*Z(:,k) - B*Z(:,k)*ritz(k);  
          end 
          resl_e(k) = norm(Res_left) /(anorm0*norm(W(:,k)));
          resr_e(k) = norm(Res_right)/(anorm0*norm(Z(:,k)));
      end
%
%     compute condition numbers       
%
      for k = 1:neig
          condnum(k) = norm(W(:,k))*norm(Z(:,k))/ abs(W(:,k)'* Z(:,k));
      end 
%
end

t2 = cputime;   % stop CPU clock 
fprintf('\nTotal CPU time in seconds = %10.5e \n',t2-t1) 

printout  % print/plot output info. 
