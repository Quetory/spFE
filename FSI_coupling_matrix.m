function [Rfsi] = FSI_coupling_matrix(ELEM, XYZ, NNs, NNf, s)

%
%          Y                    
%          |                    
%      4 ----- 3                
%      |   |   |               
%      |   +------ X          
%      |       |       
%      1 ----- 2            
%
% 4 DoF nodes , UX, UY, UZ, P

XYZs = NaN(size(XYZ));
XYZs(NNs,:) = XYZ(NNs,:);

NDOF = size(XYZs,2);
nd = length(NNs);

en = arrayfun(@(i) find(ELEM - NNs(i)==0),1:nd,'UniformOutput',false);
ind = [];
for ii = 1:length(en)
    [i,~]=ind2sub(size(ELEM),en{ii});  
    ind = [ind ;i];
end

ELEMb = ELEM(unique(ind),:);
% show_mesh(ELEMb,XYZf)
[n, ~, faces] = getsurfacenormals(ELEMb,XYZs);

% Fluid normal is opposite of structural normal
nf = -n.';

[NF,NPF] = size(faces);

if NPF == 4
    IP.XI = [-1  -1 ;...
              1  -1 ;...
              1   1 ;...
             -1   1  ]/sqrt(3);
    IP.WT = ones(1,4);
end


R = zeros(NPF*NDOF,NPF, NF);

for Nn = 1: NF
    for IN = 1:numel(IP.WT)
        XI = IP.XI(IN,:);
        WT = IP.WT(IN);
        if NPF == 4
            [N,DN] = shape2D4(XI);
            J      = [1 1 1;(DN*XYZs(faces(Nn,:),:))];
        end
        dJ   = det(J);
        R(:,:,Nn) = R(:,:,Nn) + kron(N.',eye(3))*(nf(:,Nn)*N)*dJ*WT;
    end
end

R_s_idx=(NDOF*faces(:,kron(1:NPF,ones(1,NDOF)))-kron(ones(NF,1),kron(ones(1,NPF),(NDOF-1):-1:0))).';
R_s_idx = R_s_idx(:) + s.off(1); % correct for single pressure dof, length(XYZa)*2 dofs are not present. 
X = repmat(R_s_idx,1,NPF);

%%
XYZf = NaN(size(XYZ));
XYZf(NNf,:) = XYZ(NNf,:);

en = arrayfun(@(i) find(ELEM - NNf(i)==0),1:nd,'UniformOutput',false);
ind = [];
for ii = 1:length(en)
    [i,~]=ind2sub(size(ELEM),en{ii});  
    ind = [ind ;i];
end

ELEMbf= ELEM(unique(ind),:);
% show_mesh(ELEMb,XYZf)
[~, ~, faces] = getsurfacenormals(ELEMbf,XYZf);

R_f_idx=faces + s.off(2);

Y = repmat(R_f_idx,NDOF*NPF,1);

Rfsi = sparse(X(:),Y(:),R(:),s.tot,s.tot);
% Rfsi = sparse(X(:),Y(:),R(:),s.tot,s.tot);
