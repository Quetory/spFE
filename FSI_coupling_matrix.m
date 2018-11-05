function [R] = FSI_coupling_matrix(ELEMf,XYZfluid,NNf, ELEMs,XYZstruct, NNs)

%
%          Y                    
%          |                    
%      4 ----- 3                
%      |   |   |               
%      |   +------ X          
%      |       |       
%      1 ----- 2            
%



XYZf = NaN(size(XYZfluid));
XYZf(NNf,:) = XYZfluid(NNf,:);

NDOF = size(XYZstruc,2);
nd = length(NNf);

en = arrayfun(@(i) find(ELEMf - NNf(i)==0),1:nd,'UniformOutput',false);
ind = [];
for ii = 1:length(en)
    [i,~]=ind2sub(size(ELEMf),en{ii});  
    ind = [ind ;i];
end

ELEMb = ELEMf(unique(ind),:);
% show_mesh(ELEMb,XYZf)
[n, ~, faces] = getsurfacenormals(ELEMb,XYZf);
% Fluid normal is opposite of structural normal
nf = -n.';

[NF,NPE] = size(faces);

if NPE == 4
    IP.XI = [-1  -1 ;...
              1  -1 ;...
              1   1 ;...
             -1   1  ]/sqrt(3);
    IP.WT = ones(1,4);
end


R = zeros(NPE*NDOF,NPE, NF);

for Nn = 1: NF
    for IN = 1:numel(IP.WT)
        XI = IP.XI(IN,:);
        WT = IP.WT(IN);
        if NPE == 4
            [N,DN] = shape2D4(XI);
            J      = [1 1 1;(DN*XYZ(faces(1,:),:))];
        end
        dJ   = det(J);
        R(:,:,Nn) = R(:,:,Nn) + kron(N.',eye(3))*nf(:,Nn)*N*dJ*WT;
    end
end