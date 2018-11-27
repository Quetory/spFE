function Y = matvec(problem,precond,trans,X,A,alpha,L,U,P,Q,D,B)

%
% MATVEC computes the matrix-vector product:  
% 
% If problem = 1, standard eigenvalue problem A x = lambda x, 
%   this subroutine computes 
%
%                  Y = f(A)*X   if trans = 0
%   or           
%                  Y = X.'*f(A) if trans = 1
%
%   where if precond = 0, f(A) = A;   
%                    = 1, f(A) = inv(A - alpha*I)
%
%   If precond = 1, then the LU decomposition of A-alpha*I has
%   been done, and L and U factors should be passed.  
%
% If problem = 2, generalized eigenvalue problem A x = lambda B x, 
%
%  compute
%            Y = X.'*(inv(A-shift*B)*B)
%  and
%            Y = (inv(A-shift*B)*B)*X
%
%  where f(A,B) = (inv(A-shift*B)*B, assume that A - shift*B = L*U
%  has been computed. 
%
%  ===================================================================
%

if problem == 1   % standard eigenvalue problem, A x = lambda x 
%
   if precond == 0
%
%     compute Y = X.'*A or Y = A*X 
%
      if trans == 1 
         Y = X.'*A; 
      else
         Y = A*X; 
      end 
%
   else 

%
%     compute Y = X.'*inv(A-alpha*I) or 
%             Y = inv(A-alpha*I)*X 
%     assume that A - alpha*I = L*U has been computed  
%
      if trans == 1
         w = X.'/U;   
         Y = w/L;  
      else
         w = L\X;
         Y = U\w; 
      end
%
   end 
%
else % generalized eigenvalue problem
%
%  compute
%            Y = X.'*(inv(A-shift*B)*B)
%  or 
%            Y = (inv(A-shift*B)*B)*X
%
%  assume that A - shift*B = L*U
%
   if trans == 1 
      w = (X.'*Q)/U;
      z = ((w/L)*P)/D;
      Y = z*B;
   else
      w = P*(D\(B*X));
      z = L\w;
      Y = Q*(U\z);
   end
%
end
