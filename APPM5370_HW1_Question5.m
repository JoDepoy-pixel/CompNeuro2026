%% Hodgkin-Huxley Rheobase and f-I Curve Analysis
clear; clc; close all;

%% --- Part (a): Locate Rheobase by Bisection ---
tmax = 300;
skipfrac = 0.4;
tol = 1e-4;

% Criterion: Repetitive firing is defined as rate > 0 Hz 
% (i.e., >= 2 spikes in the post-transient window t in [120, 300] ms).

gK1 = 36;
I_rheo1 = bisect_rheobase(gK1, 120, tmax, skipfrac, tol);
fprintf('Part (a): Rheobase for gK = 36 mS/cm^2 is %.3f uA/cm^2\n', I_rheo1);

% Simulate just below and just above rheobase
eps_I = 0.05;
[t_below, V_below] = hh_sim(I_rheo1 - eps_I, tmax, gK1, 120, skipfrac);
[t_above, V_above] = hh_sim(I_rheo1 + eps_I, tmax, gK1, 120, skipfrac);

% Plot comparison traces with identical axes
figure('Name', 'Part (a): Traces Near Rheobase');
subplot(2,1,1);
plot(t_below, V_below, 'b', 'LineWidth', 1.2);
title(sprintf('Just Below Rheobase (I_d = %.3f \\muA/cm^2)', I_rheo1 - eps_I));
xlabel('Time (ms)'); ylabel('V (mV)');
ylim([-20 120]); xlim([0 tmax]); grid on;

subplot(2,1,2);
plot(t_above, V_above, 'r', 'LineWidth', 1.2);
title(sprintf('Just Above Rheobase (I_d = %.3f \\muA/cm^2)', I_rheo1 + eps_I));
xlabel('Time (ms)'); ylabel('V (mV)');
ylim([-20 120]); xlim([0 tmax]); grid on;


%% --- Part (b): f-I Curve & Onset Rate (gK = 36) ---
I_vec1 = linspace(I_rheo1, 2 * I_rheo1, 40);
rates1 = zeros(size(I_vec1));

for k = 1:length(I_vec1)
    [~, ~, ~, rates1(k)] = hh_sim(I_vec1(k), tmax, gK1, 120, skipfrac);
end

% Onset rate directly above rheobase
[~, ~, ~, onset_rate1] = hh_sim(I_rheo1 + 1e-3, tmax, gK1, 120, skipfrac);
fprintf('Part (b): Onset firing rate for gK = 36 is %.3f Hz\n', onset_rate1);

figure('Name', 'Part (b): f-I Curve (gK = 36)');
plot(I_vec1, rates1, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Injected Current I_d (\\muA/cm^2)');
ylabel('Firing Rate f (Hz)');
title('f-I Curve for Classical Hodgkin-Huxley (g_K = 36 mS/cm^2)');
grid on;


%% --- Part (c): Change gK to 30 mS/cm2 ---
gK2 = 30;
I_rheo2 = bisect_rheobase(gK2, 120, tmax, skipfrac, tol);
fprintf('Part (c): Rheobase for gK = 30 mS/cm^2 is %.3f uA/cm^2\n', I_rheo2);

I_vec2 = linspace(I_rheo2, 2 * I_rheo2, 40);
rates2 = zeros(size(I_vec2));

for k = 1:length(I_vec2)
    [~, ~, ~, rates2(k)] = hh_sim(I_vec2(k), tmax, gK2, 120, skipfrac);
end

% Plot overlay of both f-I curves
figure('Name', 'Part (c): f-I Comparison');
plot(I_vec1, rates1, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'g_K = 36 mS/cm^2'); hold on;
plot(I_vec2, rates2, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'g_K = 30 mS/cm^2');
xline(I_rheo1, 'b--', sprintf('Rheobase 1: %.3f', I_rheo1), 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(I_rheo2, 'r--', sprintf('Rheobase 2: %.3f', I_rheo2), 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel('Injected Current I_d (\\muA/cm^2)');
ylabel('Firing Rate f (Hz)');
title('f-I Curve Comparison for Modified Potassium Conductance');
legend('Location', 'southeast'); grid on;


%% --- Helper Function: Bisection for Rheobase ---
function I_rheo = bisect_rheobase(gK, gNa, tmax, skipfrac, tol)
I_low = 0;
I_high = 15;

while (I_high - I_low) > tol
    I_mid = (I_low + I_high) / 2;
    [~, ~, ~, rate] = hh_sim(I_mid, tmax, gK, gNa, skipfrac);
    if rate > 0
        I_high = I_mid; % Sustained firing occurs
    else
        I_low = I_mid;  % Subthreshold / transient only
    end
end
I_rheo = I_high; % Return current just above boundary
end