clear, close all
clc
%%
XY = [ 0,0;...
       2,0;...
       2,3;...
       0,2;...
       0.4,0.4;...
       1.4,0.6;...
       1.5,2;...
       0.3,1.6];

ELEM = [ 1 2, 6, 5 ;...
        2, 3, 7, 6;
        7, 3, 4, 8;
        1, 5, 8, 4;
        5, 6, 7, 8];
NE = size(ELEM,1);
NPE = size(ELEM,2);
NDOF = size(XY,2);


show_mesh(ELEM,XY)
%%
Kelem = zeros(NPE*NDOF,NPE*NDOF,NE);
Melem = zeros(NPE*NDOF,NPE*NDOF,NE);
for e = 1 : NE
    
    [Kelem(:,:,e), B(:,:,e), L(:,:,e),Q(:,:,e)] = ElemStiffness(XY(ELEM(e,:),:), 1,0.3, 1);
    [Melem(:,:,e)] = ElemMass(XY(ELEM(e,:),:), 1.0 );
end

%% Assembly
Elements=NDOF*ELEM(:,kron( 1:NPE, ones(1,NDOF ) ) ) -kron(ones(NE,1),kron(ones(1,4), NDOF-1:-1:0 ) );
     
Y=reshape(repmat(Elements,1,8)',8,8,NE);
X=permute(Y,[2 1 3]);

K=sparse(X(:),Y(:),Kelem(:));
M=sparse(X(:),Y(:),Melem(:));

%remove small numerical round-off errors
Ku = triu(K,1);
K = Ku.'+Ku + diag(diag(K));

any(M-M.')

%% Apply constraints

D  =[0.001, -0.001;...
     0.005, -0.005;...
     0.014, -0.014;...
     0.007, -0.007].';
 
NND = size(D,2);



K11 = K(1:8,1:8);
K12 = K(1:8,9:end);
K21 = K(9:end,1:8);
K22 = K(9:end,9:end);

x = K22\(-K21*D(:));
disp(x)

 