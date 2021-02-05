function [M,K] = assemble_system_matrices(ELEM,XYZ, mat, etype)

[NE, NPE] = size(ELEM);

xyz = XYZ(ELEM(1,:),:);

if strcmpi(etype,'ACOU')

    
    [Melem] = ElemAcouMass(xyz,mat);
    [Kelem] = ElemAcouStiffness(xyz,mat);
    NDOF = 1;
    
elseif strcmpi(etype,'STRUCT')
    rho = mat.rho;
    E = mat.E;
    nu = mat.nu;
    
    eas = 0;
    [Melem] = ElemMass(xyz,rho);
    [Kelem] = ElemStiffness(xyz,E,nu,eas);
    NDOF = 3;
    
end

%% Assemble
Ke = repmat(Kelem,1,1,NE);
Me = repmat(Melem,1,1,NE);

Elements=NDOF*ELEM(:,kron(1:NPE,ones(1,NDOF)))-kron(ones(NE,1),kron(ones(1,NPE),(NDOF-1):-1:0));

Y=reshape(repmat(Elements,1,NDOF*NPE)',NDOF*NPE,NDOF*NPE,NE);
X=permute(Y,[2 1 3]);

K=sparse(X(:),Y(:),Ke(:));
M=sparse(X(:),Y(:),Me(:));


% unsr = 1/2*(real(K).'-real(K));
% unsi = 1/2*(imag(K).'-imag(K));
uns = 1/2*(K.'-K);

K = K+uns;

uns = 1/2*(M.'-M);
% unsr = 1/2*(real(M).'-real(M));
% unsi = 1/2*(imag(M).'-imag(M));

M = M+uns;
