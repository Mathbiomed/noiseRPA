

%% ================== Simulation settings ==========================
clearvars -except A
rng('shuffle');

Tmax       = 1e5;     % maximal simulation time
Fs         = 0.1;       % sampling frequency (samples per time unit)
period     = 1/Fs;
time_scale = 0:period:Tmax;
scale_size = numel(time_scale);

simulnum   = 20;     % runs per parameter set
sets_idx   = 17:17;    % run only first 10 rows
nSets      = numel(sets_idx);

% Stoichiometry (species: [X; Z1; Z2])
S = [ 1, -1,  0,  0,  0,  0,  0,  0, -1;   % X
      0,  0,  1,  0, -1,  0,  0,  0,  0;   % Z1
      0,  0,  0,  1, -1,  0,  0,  0,  0;   % Z2
      0,  0,  0,  0,  0,  1,  0, -1, 0; %Z3
      0,  0,  0,  0,  0,  0,  1, -1, 0;
      ]; 
var_num  = size(S,1);

%% ================== Storage for plots/summaries ==================
FANO_all   = zeros(nSets, scale_size);     % Fano(t) for each set
MEANX_last = zeros(nSets,1);               % <X> at final time
FANO_last  = zeros(nSets,1);               % Fano at final time

%% ================== Main loop over parameter sets ================
for s = 1:nSets
    row = sets_idx(s);

    % -------- map columns -> model parameters --------
    P.theta = 13.967363507200000;
    P.mu    = 1.574226815730000;
    P.eta   = 0.132452119018000;
    P.beta  = 1;
    P.nu    = 10;
    P.n     = -19.385364347599999;
    P.kappa = 0.050000000000000;
    P.d1 = 10;
    P.mu2 = (P.mu/P.theta)^2 + 0.01;
    P.theta2 = 1;
    P.eta2 = 1;
    P.d2 = 10;

    % -------- accumulate mean/var across runs ------
    sum_traj   = zeros(var_num, scale_size);     % Σ X
    sumsq_traj = zeros(var_num, scale_size);     % Σ X.^2

    fprintf('Set %d/%d  (theta=%.4g, mu=%.4g, eta=%.4g, beta=%.4g, nu=%.4g, n=%.4g, kappa=%.4g)\n',...
        s, nSets, P.theta, P.mu, P.eta, P.beta, P.nu, P.n, P.kappa);

    % ----- runs -----
    for run = 1:simulnum
        fprintf('  run %d/%d\n', run, simulnum);

        % one-run trajectory buffer (3 x T)
        Xtraj = zeros(var_num, scale_size);

        % initial state
        t = 0; X = [1;1;1;1;1]; idx = 1;

        while t <= Tmax
            % inhibition term
            if (X(3) == 0) & (P.n < 0)
                inhib = 1;
            elseif (X(3) == 0) & (P.n >= 0)
                inhib = 0;
            else
                inhib = (X(3)^P.n) / (X(3)^P.n + P.kappa^P.n);
            end

            % propensities 
            a = [ ...
                P.d1 * X(2);     % R1: X birth
                P.beta * X(1);    % R2: X decay
                P.mu;             % R3: Z1 birth
                P.theta * X(1);   % R4: Z2 birth
                P.eta * X(2)*X(3) % R5: Z1+Z2 -> ∅
                P.mu2;
                P.theta2*X(1)*(X(1)-1);
                P.eta2*X(4)*X(5);
                P.d2*X(5)*X(1);
            ];

            if t > Tmax/2
                a = [ ...
                P.d1 * X(2);     % R1: X birth
                P.beta *X(1);    % R2: X decay
                P.mu;             % R3: Z1 birth
                P.theta * X(1);   % R4: Z2 birth
                P.eta * X(2)*X(3) % R5: Z1+Z2 -> ∅
                P.mu2;
                P.theta2*X(1)*(X(1)-1);
                P.eta2*X(4)*X(5);
                P.d2 *20*X(5)*X(1);
            ];
            end

            a0 = sum(a);
            if a0 <= 0 || ~isfinite(a0), break; end

            % next reaction time
            tau    = -log(rand)/a0;
            t_next = t + tau;

            % record state up to t_next
            while idx <= scale_size && time_scale(idx) <= t_next
                Xtraj(:,idx) = X;
                idx = idx + 1;
            end
            if t_next > Tmax, break; end

            % choose reaction
            r2 = rand * a0;
            reaction_index = find(cumsum(a) >= r2, 1);
            X = X + S(:,reaction_index);
            t = t_next;
        end

        % fill remaining samples at final state
        while idx <= scale_size
            Xtraj(:,idx) = X;
            idx = idx + 1;
        end

        % accumulate
        sum_traj   = sum_traj   + Xtraj;
        sumsq_traj = sumsq_traj + Xtraj.^2;
    end

    % -------- mean/variance across runs --------
    mean_traj = sum_traj / simulnum;                              % E[X] over runs
    var_pop   = sumsq_traj/simulnum - mean_traj.^2;               % population var
    var_traj  = var_pop * (simulnum/(max(simulnum-1,1)));         % sample var

    % -------- Fano --------
    fano = var_traj(1,:) ./ max(mean_traj(1,:), eps);
    FANO_all(s,:) = fano;
    MEANX_last(s) = mean_traj(1,end);
    FANO_last(s)  = fano(end);
end

writematrix(mean_traj,'Fano_PNAS_Mean_NC.csv')
writematrix(var_traj,'Fano_PNAS_Var_NC.csv')

%% ================== Plots ========================================
mean_traj = readmatrix('Fano_PNAS_Mean_NC.csv');
var_traj = readmatrix('Fano_PNAS_Var_NC.csv');
fano2     = var_traj(1,:) ./ mean_traj(1,:);% Fano factor for X2

figure;
plot(time_scale, fano2, 'LineWidth',1.5);
xlabel('Time'); ylabel('Fano factor');
title('Fano factor of X'); grid on;
yline(((P.mu2/P.theta2)+(P.mu/P.theta)-(P.mu/P.theta)^2)/(P.mu/P.theta), '--', LineWidth=2)
saveas(gcf, 'PNAS_NC_FF.eps','epsc');

