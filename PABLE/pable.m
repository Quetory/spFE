function [j,szoft,T,ritzvalue,VL,VR,P,Q,ip,neig,ritz,ipritz,...
          resl,resr,omega,duality,exdual,tolconv] ...
         = pable(problem,fulldual,semidual,group,treatbd,n,maxmsz,...
                 maxeig,nb,maxit,precond,A,B,alpha,L,U,Plu,Qlu,Dlu,fanorm) 
%
%=========================================================================
%
%  Purpose 
%  =======
%
%  This M-file executes the ABLE method, namely, an adaptive block Lanczos 
%  algorithm to extract a few eigenvalues lambda or eigentriplets 
%  (lambd, y, x ) of a non-Hermitian matrix A where the eigentriplet 
%  (lambd, y, x ) approximately satisfies
%
%             y'*A = lambda*y'  and   A*x = x*lambda
%
%  The user is required to furnish a M-file matvec.m to calculate the 
%  products        
%
%             Y = f(A)*X 
%  or 
%             Y = X'*f(A) 
%
%  where f(A) = A if no spectral transformation is selected, or  
%  f(A) = inv(A-alpha*I) if the shift-and-invert spectral transformation
%  is selected. 
%
%  This M-file can also solve the generalized eigenvalue problem
%
%             y'*A = lambda*y'*B    and    A*x = lambda*B*x 
%
%  by using the shift-and-invert the spectral transformation to reduce it 
%  to the standard eigenvalue problem
%
%             z'*C = mu*z'  and  C*x = mu*x 
%
%  where 
%            C  = inv(A-alpha*B)*B, 
%            z' = y'*B,  
%            mu = 1/(lambda-alpha). 
% 
%
%  Input Arguments 
%  ===============
%
%  problem   integer
%            if problem = 1, solve the standard eigenvalue problem 
%            if problem = 2, solve the generalized eigenvalue problem 
%
%  fulldual: logical 
%            Lanczos process with or without rebiorthogonalization
%            at every iteration (default: false). 
%
%  semidual: logical 
%            monitoring the loss of duality and performing the
%            rebiorthogonalization if necessary (default: true).  
%
%  group:    logical 
%            adapting the blocksize to the largest order of converged
%            Ritz values (default: true).  
%
%  treatbd:  logical 
%            treating the breakdown (default: true). 
%
%            The ABLE method can be implemented at different levels by 
%            the choice of the logical variables fulldual, semidual, 
%            group and treatbd. 
%
%            Level     fulldual     semidual     group     treatbd
%            -----------------------------------------------------
%              1          F            F           F          F  
%              2          T            F           F          F  
%              3          F            T           F          F  
%              4          T            T*          T          T  
%
%            F = False, T = True. At level 4, only one of options in fulldual
%            and semidual is requried. In the level 1 implementation, 
%            eigenvalues only are computed. 
%
%            Note: In the current implementation, these logical variables
%            are replaced by integers 0 and 1. O for false, 1 for true. 
%
%  n         integer
%            the order of test matrix(ces)
%
%  maxmsz    integer
%            2*n*maxmsz of memory for storing the computed Lanczos vectors,
%            only needed for the level of implementation >= 2.
%
%  maxeig    integer
%            number of desired eigenvalues 
%
%  nb        integer 
%            number of initial blocksize 
%
%  maxit     integer
%            max number of iteration 
%            
%  precond   integer
%            specify the selected spectral transformation scheme
%              precond = 1, the shift-and-invert 
%            others are described in the file matvec, but has not been tuned.
%
%  A, B      double arrays
%            matrices A (and B) in the above defined eigenvalue problem, 
%            A (and B) is storaged in sparse matrix-format. 
%
%  alpha     double precision 
%            the selected shift 
%
%  L, U      L and U factors of the LU decomposition of (A - alpha*B), where
%            B is an identity matrix for the standard eigenvalue problem
%
%  fanorm    double precision 
%            norm of f(A) matrix or max(norm(A),norm(B))
%
%
%  Output Arguments 
%  ================
%
%  j         integer
%            number of Lanczos steps 
%
%  szoft     integer
%            size of reduced (block) tridiagonal matrix T
% 
%  T         double precision array
%            the reduced (block) tridiagonal matrix 
%
%  ritzvalue double precision array 
%            eigevalues of T (Ritz values)  
%
%  VL, VR    double precision arrays
%            the left and right eigenvectors of T
%
%  P, Q      double precision arrays
%            P and Q store the computed left and right Lanczos vectors 
%            if Level >= 2. At level 1, P and Q are not referenced.
%
%  ip        integer array
%            pointer array, ip(k) is the pointer to the beginning of the 
%            kth block in the block tridiagonal matrix T and Lanczos vectors
%            P and Q if they are saved. Specifically, the blocks of P can be 
%            identified as P( ip(1):ip(2)-1, ip(2):ip(3)-1, .... ).
%
%  neig      integer
%            number of converged Ritz values
%
%  ritz      double precision array
%            converged Ritz values (approximate eigenvalues) 
% 
%  ipritz    integer array
%            pointer array, ritzvalue(ipritz(j)) is a convergent Ritz value
%            and is stored in the position ritz(j). 
%
%  resl,resr double precision arrays
%            left and right alg_residual norms of all Ritz pairs
%
%  omega     double precision array
%            omega(k) = smallest singular values of (P_k^T*Q_k), where P_k
%            and Q_k are the k-th blocks of left and right eigenvectors. 
%            It is the parameter to measure the stability, which is used to 
%            detect the breakdown, or the growth of the norm of Lanczos 
%            vectors.
%
%  duality   double precision array
%            the estimated biorthogonality (duality) of Lanczos vectors     
%            (only returned if level >= 3) 
%
%  exdual    double precision array
%            the exact biorthogonality (duality) of Lanczos vectors     
%            (will be removed in the future, just for testing purpose)
%
%  condnum   double precision array
%            condition number of converged Ritz values, only if level >= 2.  
%
%  tolconv   double precision
%            stopping criterion for convegence test.   
% 
%
%  Internal control parameters, may be altered by the user:
%  ========================================================
%
%  nbmax:    maximum blocksize, say 1 <= nbmax <= 8 (default: 6)
%  maxmsz:   if fulldual or semidual is chosen to be performed, maxmsz
%            specifies the maximum memory size of n-by-maxmsz arrays
%            P and Q to store the computed Lanczos vectors
%            In Level 1 ABLE method, arrays P and Q are not needed.
%  tolconv:  threshold for testing convergence of Ritz values.
%  toldef:   threshold for testing deflation among Lanczos vectors in a block
%  tolbd:    threshold for testing break down
%  tolclust: threshold for testing clustering.
%
%
%  Matlab intrinsic functions used in PABLE:
%  ========================================
%
%  qr, orth, svd, norm, randn, max, min, conj, abs
%
%  External M-files routines called by PABLE:
%  ========================================
%
%  convtst, cure, semi, tsmgs
%
%
%  Additional information on some key Scalars and Arrays used in PABLE
%  ===================================================================
%
%  pqnorm    array, dimension(j) 
%            pqnorm(k) is the F-norm of the kth Lanczos vectors and P_k 
%            and Q_k.  
%
%  P1, Q1, P2, Q2, R, S 
%            double precision arrays, dimension( n,nb )
%            workspace to be used in the basic 3-term recurrences 
%
%  X1, X2, Y1, Y2
%            double precision arrays, dimension( j,nb )
%            workspace used for monitoring the loss of duality. Only 
%            referenced if semidual = true. 
%
%===========================================================================
%
%
% Set up internal control parameters 
% ==================================
%
nbmax=8;                 % limit the maximum blocksize
tolconv = 1e-8;          % tol for testing convergence of Ritz values.
toldef = 1e-10;          % tol for testing deflation among Lanczos vectors 
                         % in a block
tolclust =1e-8;          % tol for testing clustering.
if fulldual == 0 && semidual == 0 % set tol for testing break down 
   tolbd = 10*n*eps;     % Level 1 
else
%   tolbd = 1e-8;        % Level 2 and up.  
   tolbd = 10*n*eps; 
end
  
%
% Set up initial random vectors if necessary 
% ==========================================
 
initvec = 0;             % initial starting vectors   (should be set as input?)
                         % if initvec = 1, the initial vectors P1 and Q1 
                         % have to be passed in the calling sequence of this
                         % file   
if initvec == 0          % generate initial starting vectors 
   rng(1993);            % setup the random number seed 
   P1 = randn(n,nb); 
   Q1 = P1;
   P1 = orth(P1); 
   Q1 = P1; 
   omega(1) = 1; 
else                     % use given initial starting vectors 
   [u,sigma,v] = svd(P1.'*Q1); 
   sigma = diag(sigma);
   omega(1) = min(sigma); 
   if omega(1) == 0 
      fprintf('initlal omega(1) is zero, restart with other initial vectors\n')
      exit
   end 
   sqrtsigma = sqrt(sigma); 
   P1 = P1 * ( conj(u)/diag(sqrtsigma) );
   Q1 = Q1 * ( v/diag(sqrtsigma) );
end 
%
pqnorm(1) = max( norm( P1,'fro'), norm( Q1,'fro' ) );
if fulldual == 1 || semidual == 1  % save Lanczos vectors if desired
   P = P1;
   Q = Q1;
%  usedmemory = nb;      % a counter for the usage of total memory for 
                         % storing the Lanczos vectors  
end 
if semidual == 1         % set up initial value if semi-duality is required
   duality(1) = eps; 
   exdual(1) = eps;
else
    duality = [];
    exdual = [];
end
 
R = (matvec(problem,precond,1,P1,A,alpha,L,U,Plu,Qlu,Dlu,B)).';     
S = matvec(problem,precond,0,Q1,A,alpha,L,U,Plu,Qlu,Dlu,B);   

 
% Initialize data structure for variable blocks.
% =============================================
 
j = 1;
ip(j) = 1;                    % pointer to the beginning of the jth block.
ip(j+1) = ip(j) + nb;
szoft = 0;                    % order of tridiagonal matrix T_j
neig = 0;                     % number of converged RR values 
isolve = 0; 
minstep = 1;                 
jsolve = min( n/nb, maxeig+minstep );

% Note: parameters isolve, minstep and jsolve control the frequence of solving
%       the eigenproblems of the smaller matrix T and test for convergence. 
%       When the Lanczos step j = jsolve, then isolve is step to 1 and the 
%       eigenproblem of T is solved. If we have a good preconditioning 
%       (spectral transformation), minstep should be small. 
%


%  =================== Beginning of main loop ===============================
 

while j <= maxit

%
%  ------------------------ Compute the (j,j) entry (block) of T 
%
   aj = P1.'*S; 
   T( ip(j):ip(j+1)-1, ip(j):ip(j+1)-1) = aj;
   szoft = szoft + nb;
%
%  ------------------------ Continue 3 term recurrence
%
   P2 = R - P1 *aj.'; 
   Q2 = S - Q1 *aj;   
%
%  ------------------------ Compute eigen-decomposition of T and 
%                           Test for convergence if required 
%
   if ( j == jsolve || j == maxit )
      isolve = 1; 
   end 
%
   isolve = 1;              % solve eig of T at every step, otherwise, common
                            % out this.     
%  
   if isolve == 1 && j > 1 
      [neig,ritz,ipritz,ritzvalue,VL,VR,resl,resr] = ...
        convtst(fulldual,semidual,nb,szoft,T,P2,Q2,fanorm,tolconv,tolclust); 
%
      isolve = 0; 
      jsolve = j + minstep;  
%
%      fprintf( 'Solve eig(T) at %3i step, szoft = %3i, conv RR = %3i \n',...
%               j,szoft,neig );
%
   end
%
   if neig >= min(n,maxeig)  
      fprintf( 'desired %3i eigenvalue(s) have been found\n',min(n,maxeig) )  
      break;                % successful exit from the main loop, info = 0  
   elseif j == maxit 
      fprintf( 'maximal number of Lanczos steps have been run out\n' )
      fprintf( 'number of converged eigenvalues is %3i \n', neig )
      break; 
   end
%
%  ------------------------ Extended local duality (dualize the previous
%                           two pairs of Lanczos vectors P_{j-1} and P_j
%                           and Q_{j-1} and Q_j. 
%
   if j > 2 && (fulldual == 1 || semidual == 1)
      ki = ip(j-1);
      kf = ip(j)-1;
      P2 = P2 - P(:,ki:kf) *( Q(:,ki:kf).' * P2 );
      Q2 = Q2 - Q(:,ki:kf) *( P(:,ki:kf).' * Q2 );
   end
   P2 = P2 - P1 * ( Q1.' * P2 ); 
   Q2 = Q2 - Q1 * ( P1.' * Q2 ); 
%
%  ------------------------ Scaling - (rank-revealing) QR decomposition 
%
   [P2, bj] = qr(P2,0);
   bj = bj.';
   [Q2, cj] = qr(Q2,0);
%
%  ------------------------ Detect rank deficient of P2 and Q2 (R and S
%                           in the document)  
%
%  If P2 and/or Q2 is rank deficient (bj and/or cj is singular),  we have
%  the following possible cases and choices:  
%
%     (1) if P2 = 0 (Q2 = 0), left (right) invariant subspace is obtained, 
%         stopping the computation or if desired, restart with the new 
%         initial vectors (dualized against previous Lanczos vectors)   
%     (2) otherwise, choose the columns of P2 (Q2) so that P2 (Q2) is not
%         rank deficient, and biorthogonalize to all previous blocks in 
%         Q (P) (blocksize is unchanged)  
%     (3) or reduce the blocksize 
%  Our current implementation is (3). There are some alternatives to do.  
%
   rankbj = rank(bj,nb*(eps)^(2/3));  
   rankcj = rank(cj,nb*(eps)^(2/3));  
   if ( rankbj < nb || rankcj < nb )
      fprintf('current block vector(s) are rank deficient j = %4i\n',j)
      if ( rankbj == 0 || rankcj == 0 )
          fprintf('P2 = 0 or Q2 = 0, invariant subspace(s) found');
          break;           
      else 
         nb = rankbj;      % reduce the blocksize 
         P2 = P2(:,1:nb);
         Q2 = Q2(:,1:nb);
         bj = bj(:,1:nb);   
         cj = cj(1:nb,:);  
      end
   end
%
%  ------------------------ Compute the scaling factors (matrices)  
%
   [u,sigma,v] = svd(P2.'*Q2);  
   sigma = diag(sigma);
   omega(j+1) = min(sigma);  % parameter for observing the stability
%
   if( omega(j+1) <= 10*n*eps && nb >= nbmax )  
      fprintf(' method fails, incuriable breakdown ...........') 
      break;                 
   end 
%
%  ------------------------ Test for breakdown.
%
   breakdown = 0; 
   delb=0;
   if omega(j+1) < tolbd
      breakdown = 1;   % turn on breakdown flag  
      fprintf('breakdown occurs at %5i step, minofvu = %9.5e\n',j+1,omega(j+1))
      delb = sum( sigma < tolbd );  
   end
%
   if ( breakdown == 1 && treatbd == 0 ) 
      fprintf('breakdown, no treatment requested, ...., method fails ')
      break;         
   end
%
%  ------------------------ Detect the clustering of converged RR values
%
   delc=0;   % parameter for determining the blocksize increment in the
              % presence of larger order of the clustering RR values 
%
   if group == 1 && neig > 1
      for i = 1:neig
          ritz1 = abs( ritz - ritz(i)*ones(neig,1) );
          delc = max( delc, sum( abs(ritz1) < tolclust*abs(ritz(i)) ) );
      end
      delc = max(0, delc - nb);            % the increament block size
      if delc > 0
         fprintf( 'a cluster is found at step = %5i \n',j )
         fprintf( 'proposed increament of the blocksize = %5i\n',delc )
      end
   end
%
%  ------------------------ Cure the breakdown and/or adapt to the largest 
%                           cluster, if requested. 
%
   del = 0;
   if ( group == 1 || ( breakdown == 1  && treatbd == 1 ) ) 
       del = min( max(delc, delb), nbmax - nb );
       if del >= 1
	  fprintf('%3i new-start vecs., new b-size = %3i\n', del, nb+del )
%
	  [del,P2,Q2,u,sigma,v,omega,tolbd] = ...
		  cure(n,nb,del,j,P2,Q2,P,Q,u,sigma,v,ip,omega,tolbd); 
%
          if omega(j+1) < 10*n*eps        % procedure CURE fails 
             fprintf('cure procedure fails, cannot cure breakdown' ) 
             fprintf('j = %4i, omega = %9.4e\n',j,omega(j+1))  
             break;     % exit the main loop, info = .. 
          end
%
%         Update bj and cj
%
          if del > 0
             bj=[bj, zeros( size(bj,1), del )];
             cj=[cj; zeros( del, size(cj,2) )];
          end
       end
   end
%
   nb = nb + del; %  update blocksize nb 
%
%  ------------------------ Compute (j+1)th left and right Lanczos vectors
%
   sqrtsigma = sqrt(sigma);
   bj = bj*u*diag(sqrtsigma) ;       
   cj = diag(sqrtsigma)*v'*cj;
   P2 = P2 * ( conj(u)/diag(sqrtsigma) );
   Q2 = Q2 * ( v/diag(sqrtsigma) );
   pqnorm( j+1 ) = max( norm( P2,'fro'), norm( Q2,'fro' ) );
%
   ip(j+2) = ip(j+1) + nb;
   T( ip(j):ip(j+1)-1, ip(j+1):ip(j+2)-1) = bj;
   T( ip(j+1):ip(j+2)-1, ip(j):ip(j+1)-1) = cj;
%
%  ------------------------ Maintain full or semi-Duality if requested
%
   if fulldual == 1        % Level 2 or 4 implementation   
%
%     call two-sided modified Gram-Schmidt (TSMGS)
%
      [P2,Q2] = tsmgs( j,nb,P2,Q2,P,Q,ip );
%
   elseif semidual == 1    % Level 3 or 4 implementation 
%
%     call procedure SEMI to monitor the loss of duality and 
%     correct the duality if necessary
%
      if j == 1 

         X2 = P1.' * Q2;       Y2 = P2.' * Q1;
         X1 = zeros(size(X2)); Y1 = zeros(size(Y2));
         wjj = Y2;             zjj = X2;
         duality(j+1) = eps;
         exdual(j+1) = eps;

      else 

         [P1,Q1,P2,Q2,P,Q,X1,X2,Y1,Y2,wjj,zjj,isolve,duality,exdual ] = ...
             semi( j,fanorm,T,P1,Q1,P2,Q2,P,Q,ip,pqnorm,X1,X2,Y1,Y2,...
                   wjj,zjj,duality,exdual ); 
      end 
%
   end                                                                   
%
%  ------------------------ Save Lanczos vectors if desired.
%
   if fulldual == 1 || semidual == 1
      P(:, szoft+1:szoft+nb ) = P2;
      Q(:, szoft+1:szoft+nb ) = Q2;
%
%     count the usage of memory and modified the maximal number of 
%     iterations allowed.
%
%     usedmemory = usedmemory + nb; 
%     maxit = min( maxit, j + floor( (maxmsz - usedmemory) / nb ) );
%
   end 
%
%  ------------------------ Prepare for the next iteration 
%
   R = matvec(problem,precond,1,P2,A,alpha,L,U,Plu,Qlu,Dlu,B); 
   R = ( R - cj * P1.' ).'; 
   S = matvec(problem,precond,0,Q2,A,alpha,L,U,Plu,Qlu,Dlu,B); 
   S =   S  - Q1 * bj;      
%
   Q1 = Q2;
   P1 = P2;
%
   j = j+1; 
%
end 
