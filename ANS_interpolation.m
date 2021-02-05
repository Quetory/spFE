function [Bans, dJ] = ANS_interpolation(DN, XI, XYZ)

    % XYZ0 = [-1  -1  -1;...
        % 1  -1  -1;...
        % 1   1  -1;...
        % -1   1  -1 ;...
        % -1  -1   1;...
        % 1  -1   1;...
        % 1   1   1;...
        % -1   1   1  ];

    [NPE, NDOF] = size(XYZ);

    Bans = zeros(6, NPE * NDOF);

    G = (DN * XYZ).';
    dJ = det(G);

    Bans(1, :) = kron(DN(1, :), G(:, 1).');
    Bans(2, :) = kron(DN(2, :), G(:, 2).');

    Bans(4, :) = kron(DN(2, :), G(:, 1).') + kron(DN(1, :), G(:, 2).');

    %ANS interpolatie punten
    XA = [-1  0 0 ;...
           0 -1 0 ; ...
           1  0 0 ; ...
           0  1 0 ];

    XL = [-1 -1 0 ; ...
           1 -1 0 ; ...
           1  1 0 ; ...
          -1  1 0 ];

    DN3L = zeros(3, 8, 4);
    DN3A = zeros(3, 8, 4);
    GL = zeros(3, 3, 4);
    GA = zeros(3, 3, 4);

    for a = 1:4
        [~, DN3L(:, :, a)] = shape3D8(XL(a, :));
        GL(:, :, a) = (DN3L(:, :, a) * XYZ).';
    end

    for a = 1:4
        [~, DN3A(:, :, a)] = shape3D8(XA(a, :));
        GA(:, :, a) = (DN3A(:, :, a) * XYZ).';
    end

    % E = {e11 e22 e33 2*e12 2*e23 2*e13}
    E13B = 1/2 * (kron(DN3A(1, :, 2), GA(:, 3, 2).') + kron(DN3A(3, :, 2), GA(:, 1, 2).'));
    E13D = 1/2 * (kron(DN3A(1, :, 4), GA(:, 3, 4).') + kron(DN3A(3, :, 4), GA(:, 1, 4).'));
    E23A = 1/2 * (kron(DN3A(2, :, 1), GA(:, 3, 1).') + kron(DN3A(3, :, 1), GA(:, 2, 1).'));
    E23C = 1/2 * (kron(DN3A(2, :, 3), GA(:, 3, 3).') + kron(DN3A(3, :, 3), GA(:, 2, 3).'));

    Bans(5:6, :) = 1/2 * [  (1 - XI(1)) * 2 * E23A + (1 + XI(1)) * 2 * E23C ;...
                            (1 - XI(2)) * 2 * E13B + (1 + XI(2)) * 2 * E13D];
    % E33
    for i = 1:4
        Bans(3, :) = Bans(3, :) + 1/4 * (1 + XL(i, 1) * XI(1)) .* (1 + XL(i, 2) * XI(2)) * kron(DN3L(3, :, i), GA(:, 3, i).');
    end
