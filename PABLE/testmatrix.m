% test problems, to be called by pbdy.m
%
if testprob == 1, 
   testmatname = 'real random dense symmetric matrix';
   problem = 1; B = 0;   % standard eigenvalue problem
   randn('seed',1993);
   temp = randn(50);
   A = temp'*temp; 
   aeig = eig(A);
   eig_ext = 1;          % has exact eigenvalues
   precond = 0;
   alpha = 0; 
elseif testprob == 2,
   testmatname = 'complex symmetric matrix';
   problem = 1; B = 0;   % standard eigenvalue problem
   randn('seed',1993);
   temp = randn(50) + i*randn(50);
   A = temp.'*temp;
   aeig = eig(A);
   eig_ext = 1;             % has exact eigenvalues
   precond = 0;
   alpha = 0; 
elseif testprob == 3, 
   testmatname = 'random real sparse matrix'; 
   problem = 1; B = 0;      % standard eigenvalue problem
   randn('seed',1993);
   A=sprandn(50,50,.10);
   aeig = eig(full(A));
   eig_ext = 1;             % has exact eigenvalues
   precond = 0; 
   alpha = 0;
elseif testprob == 4, 
   testmatname = 'Brusselator Model I';  
   problem = 1; B = 0;      % standard eigenvalue problem
   [A,aeig] = brussel(100); % matrix dimension has to be even 
   eig_ext = 1;             % has exact eigenvalues
   precond = 1; 
   alpha = 0;
elseif testprob == 5,
   testmatname='Random generalized eigenvalue problem';
   problem = 2;
   A = randn(50);
   B = randn(50);
   aeig = eig(A,B);
   eig_ext = 1;
   precond = 1;             % choice of preconditioning
   alpha = 0;
end 
