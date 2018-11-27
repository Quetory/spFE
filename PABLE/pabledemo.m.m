

%  Demo I 
%  ======

% Set the following input parameters 

clear; format short e; format compact;  % clean up screen
testprob = 1; % select a test problem 
maxeig = 4;   % number of required eigenvalues
nb = 1;       % initial blocksize
maxit = 100;   % maximum number of Lanczos-iter, altered by maxmsz
matbal = 1;   % whether balance the test matrix (only for standard eigenproblem)
fulldual = 1; % use full-duality if fulldual = 1
semidual = 0; % maintain semi-duality if semidual = 1
group = 0;    % adapt to the order of the largest cluster if group = 1
treatbd = 0;  % cure breakdown if treatbd = 1.

testmatrix;   % include the selected test problem

figure(1); 
pbdy  % call test driver routine of pable    

fprintf('.... Hit Return For The Next Demonstate Example \n')
pause

%  Demo II
%  ========

% Set the following input parameters

clear; format short e; format compact;  % clean up screen
testprob = 2; % select a test problem
maxeig = 4;   % number of required eigenvalues
nb = 1;       % initial blocksize
maxit = 100;   % maximum number of Lanczos-iter, altered by maxmsz
matbal = 1;   % whether balance the test matrix (only for standard eigenproblem)
fulldual = 0; % use full-duality if fulldual = 1
semidual = 1; % maintain semi-duality if semidual = 1
group = 0;    % adapt to the order of the largest cluster if group = 1
treatbd = 0;  % cure breakdown if treatbd = 1.
%
testmatrix;   % include the selected test problem

figure(2); 
pbdy  % call test driver routine of pable

fprintf('.... Hit Return For The Next Demonstate Example \n')
pause

%  Demo III 
%  ========

% Set the following input parameters

clear; format short e; format compact;  % clean up screen
testprob = 3; % select a test problem
maxeig = 4;   % number of required eigenvalues
nb = 1;       % initial blocksize
maxit = 100;   % maximum number of Lanczos-iter, altered by maxmsz
matbal = 1;   % whether balance the test matrix (only for standard eigenproblem)
fulldual = 0; % use full-duality if fulldual = 1
semidual = 1; % maintain semi-duality if semidual = 1
group = 0;    % adapt to the order of the largest cluster if group = 1
treatbd = 0;  % cure breakdown if treatbd = 1.
%
testmatrix;   % include the selected test problem

figure(3);
pbdy  % call test driver routine of pable

fprintf('.... Hit Return For The Next Demonstate Example \n')
pause

%  Demo IV 
%  =======

% Set the following input parameters


clear; format short e; format compact;  % clean up screen
testprob = 4; % select a test problem
maxeig = 4;   % number of required eigenvalues
nb = 1;       % initial blocksize
maxit = 100;   % maximum number of Lanczos-iter, altered by maxmsz
matbal = 1;   % whether balance the test matrix (only for standard eigenproblem)
fulldual = 0; % use full-duality if fulldual = 1
semidual = 1; % maintain semi-duality if semidual = 1
group = 0;    % adapt to the order of the largest cluster if group = 1
treatbd = 0;  % cure breakdown if treatbd = 1.
%
testmatrix;   % include the selected test problem

figure(4);
pbdy  % call test driver routine of pable

fprintf('.... Hit Return For The Next Demonstate Example \n')
pause

%  Demo V
%  ======

% Set the following input parameters

clear; format short e; format compact;  % clean up screen
testprob = 5; % select a test problem
maxeig = 10;   % number of required eigenvalues
nb = 1;       % initial blocksize
maxit = 100;   % maximum number of Lanczos-iter, altered by maxmsz
matbal = 1;   % whether balance the test matrix (only for standard eigenproblem)
fulldual = 0; % use full-duality if fulldual = 1
semidual = 1; % maintain semi-duality if semidual = 1
group = 0;    % adapt to the order of the largest cluster if group = 1
treatbd = 0;  % cure breakdown if treatbd = 1.
%
testmatrix;   % include the selected test problem

figure(5);
pbdy  % call test driver routine of pable

