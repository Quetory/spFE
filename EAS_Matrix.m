function [M] = EAS_Matrix(XI,T,dJ)

if size(T,1)==3
    
    s = XI(1);
    t = XI(2);
    M = zeros(3,4);
    M(1:3,1:3) = diag([s,t,s]);
    M(3,4) = t;
    
    M = 1/dJ*T*M;
else
    r = XI(1);
    s = XI(2);
    t = XI(3);
    
    % EAS11 Klinkel
%     M = zeros(6,11);
%     M(1,1:2) = [r r*s];
%     M(2,3:4) = [s r*s];
%     M(3,5:8) = [t t*r t*s r*s*t];
%     M(4,9:11)= [r s r*s];

    
%     %EAS24
    M = zeros(6,24);
    M(1:3,1:3)   = diag([r s t]);
    M(4:6,4:9)   = blkdiag([r s], [r t], [s t]);
    M(4:6,10:15) = blkdiag([r*t s*t], [r*s s*t],[r*s r*t]);
    M(1:3,16:21) = blkdiag([r*s r*t],[r*s s*t],[r*t s*t]);
    M(4:6,22:24) = diag([r*s r*t s*t]);
%     M(:,  25:30) = r*s*t*eye(6);
    
    M = 1/dJ*T*M;
end


