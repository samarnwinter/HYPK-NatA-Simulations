% Analytical kinetic model for cotranslational N-terminal acetylation (NTA)
% -------------------------------------------------------------------------
% This MATLAB script is a master equation based model
% describing the kinetic network of NatA–ribosome–HYPK interactions
% during translation. The model calculates the  steady state probabilities  of
% different kinetic states as functions of the translation elongation rate
% (ω), RNC NatA association (k12, k43), dissociation (k21, k34), acetylation (k23), and coding sequence
% length (L).
%----------------------------------------------------------------------------
% Usage:
%   Run directly in MATLAB:
%       HYPK_NatA_acetylation_analytical_model
%
% Author: Inayat Ullah Irshad
% Date: October 2025

%==============================================================================


clear; clc; close all;

enx = [0.1, 1];         %===== %=== Backward rate 
xx  = [4, 0.2];         %===== %=== forward rate
omg = 1:10;             %=== codon translation rate

n = 353;                % number of codon positions / steps
x = 1:n;

% ==============Gaussian mixture used to shape k3 along the chain========

a1 = 0.79;  a2 = 0.20;   a3 = 0.2950;
b1 = 105.0; b2 = 116.2;  b3 = 154.5404;
c1 = 6.414; c2 = 10.6473; c3 = 28.1588;

y = a1 .* exp(-((x - b1) ./ c1).^2) + ...
    a2 .* exp(-((x - b2) ./ c2).^2) + ...
    a3 .* exp(-((x - b3) ./ c3).^2);
y = y ./ max(y);        % normalize to max 1

Avg_k = [0, 0.001373, 0.002281, 0.845, 0.2825, 0.17, 0.285, 0.075, 0.046];

nat = [10 20 40];       % [NatA] 
R_total_in = 70e-9;     % [RNC]


N_total_in = nat(1) * 1e-9;

acetyl = zeros(length(omg), 2);
Nat_conc = zeros(length(omg), 2);
on_rates = zeros(length(omg), 2);

for mm = 1:length(omg)
    for k = 1:2
        % kinetic rates
        k1 = xx(k) * 1e8;                     %=== forward rate
        k2 = enx(k);                          %=== Backward rate
        kd = k2 / k1;                         %=== eq. dissociation constant
        b = R_total_in - N_total_in + kd;

        c = -kd * N_total_in;
        d = b^2 - 4 * c;

        if d > 0
            NatA_s1 = (-b + sqrt(d)) / 2;    %===== free NatA concentration
        else
            NatA_s1 = 0;    % safe fallback
        end

        Nat_conc(mm, k) = NatA_s1 * 1e9;              %===== data_save free NatA concentration
        on_rates(mm, k)  = NatA_s1 * xx(k) * 1e8;     %===== lly, on rates k1

        kk1_scalar = k1 * NatA_s1;                    %===== on rates
        kk2_scalar = k2;                              %===== off rates
        kk3_scalar = max(Avg_k);                      %===== Enzymatic rxn rate at peak

        % construct Enzymatic rate for each condon position i = 1..n
        
        rate_ac = cell(1, n);
        for i = 1:n
            rate_ac{i} = [kk1_scalar, kk2_scalar, kk3_scalar * y(i), kk2_scalar, kk1_scalar];
        end

        omega_val = omg(mm);  % codon elongation rates
        pLen = n;                   %===== CDS len    
        p1 = zeros(1, pLen);        %===== RNC  
        p2 = zeros(1, pLen);        %===== RNC-NatA
        p3 = zeros(1, pLen);        %===== NatA*Ac-RNC
        p4 = zeros(1, pLen);        %===== Ac-RNC

        p1(1) = 1;

        for idx = 2:pLen
            rate2 = rate_ac{2};
            kf1 = rate2(1);
            kf2 = rate2(2);
            kf3 = rate2(4);
            kf4 = rate2(5);

            km1 = rate_ac{idx-1}(3);

            omega = omega_val;

            kp = zeros(4,4);
            kp(1,:) = [-(omega + kf1),      kf2,                 0,      0];
            kp(2,:) = [kf1,               -(omega + kf2 + km1), 0,      0];
            kp(3,:) = [0,                  km1,               -(omega + kf3), kf4];
            kp(4,:) = [0,                  0,                   kf3,   -(omega + kf4)];

            wp = -omega * eye(4);

            mt = kp \ wp;

            prev_vec = [p1(idx-1); p2(idx-1); p3(idx-1); p4(idx-1)];
            pvec = mt * prev_vec;

            p1(idx) = pvec(1);
            p2(idx) = pvec(2);
            p3(idx) = pvec(3);
            p4(idx) = pvec(4);
        end

        p34 = p3 + p4;
        acetyl(mm, k) = p34(300);
    end
end

figure(1);
hold on;
plot(omg, acetyl(:,1), '-o', 'LineWidth', 2.0, 'Color', 'blue', 'MarkerSize', 10, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue'); 
plot(omg, acetyl(:,2), '-s', 'LineWidth', 2.0, 'Color', 'red', 'MarkerSize', 10, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
%xlabel('\rm\omega (s^{-1})', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
xlabel('ω (s^{-1})', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Acetylation', 'Interpreter', 'tex', 'FontSize', 16,  'FontName', 'Arial');
leg = legend({'- HYPK', '+ HYPK'}, 'Location', 'best', 'FontSize', 16, 'Interpreter', 'tex');
set(leg, 'Box', 'off', 'FontName', 'Arial');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.015], 'FontName', 'Arial');
hold off;
