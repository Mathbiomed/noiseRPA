% Gillespie simulation for APIF Class 1 with Degradation Inhibition

clear; close all; clc;
rng('shuffle');

%-------------------------------
% Simulation time grid
%-------------------------------
Tmax       = 1e5;    % maximal simulation time
Fs         = 0.1;     % sampling frequency
period     = 1/Fs;
time_scale = 0:period:Tmax;
scale_size = length(time_scale);
simulnum   = 500;    % number of independent runs

%-------------------------------
% Model parameters
%-------------------------------
P.delta   = 150; %%%%%
P.k       = 3; 
P.gamma_1 = 2;
P.gamma_2 = 7;

P.mu      = 10; %%
P.theta   = 2;  %%
P.k_1     = 2;  %%
P.eta     = 1e4;     % degradation inhibition rate
P.kappa_1 = 1;  %%

%-------------------------------
% Stoichiometry matrix
% species: [X1; X2; Z1; Z2]
% reactions R1...R8
S = [ ...
     0,  -1,   0,   0,   0,   0,   1,  -1;   % X1
     1,   0,  -1,   0,   0,   0,   0,   0;   % X2
     0,   0,   0,   1,   0,  -1,   0,   0;   % Z1
     0,   0,   0,   0,   1,  -1,   0,   0    % Z2
];
var_num  = size(S,1);
reac_num = size(S,2);

% Preallocate storage: runs x species x time
simuldataf = zeros(simulnum, var_num, scale_size);

%-------------------------------
% Run Gillespie-sampled trajectories (PARFOR-safe)
%-------------------------------
parfor run = 1:simulnum
    fprintf("%d\n", run);  % (원하면 유지 가능, 병렬 출력은 섞여 보일 수 있음)
    t   = 0;
    X   = [10; 10; 0; 0];
    idx = 1;

    % 로컬 버퍼에 저장 후 마지막에 한 번에 simuldataf(run,:,:)에 기록
    traj = zeros(var_num, scale_size);

    while t <= Tmax
        % compute propensities
        a = [ ...
            P.k_1 * X(1);                                % R1
            P.gamma_1 * X(1);                            % R2
            P.gamma_2 * X(2);                            % R3
            P.mu;                                        % R4
            P.theta * X(2);                              % R5
            P.eta * X(3) * X(4);                         % R6
            P.k   * X(3);                                % R7
            P.delta * X(1) * X(2) / (X(1)+P.kappa_1);    % R8
        ];

        if t > Tmax/2
            a = [ ...
                P.k_1 * 10 * X(1);                            % R1
                P.gamma_1 * X(1);                        % R2
                P.gamma_2 * 10 * X(2);                   % R3
                P.mu;                                    % R4
                P.theta * X(2);                          % R5
                P.eta * X(3) * X(4);                     % R6
                P.k  * X(3);                        % R7
                P.delta * X(1) * X(2) / (X(1)+P.kappa_1);% R8
            ];
        end

        a0 = sum(a);
        if a0 <= 0
            break;
        end

        % time to next reaction
        tau    = -log(rand)/a0;
        t_next = t + tau;

        % record state until next event
        while idx <= scale_size && time_scale(idx) <= t_next
            traj(:,idx) = X;
            idx = idx + 1;
        end
        if t_next > Tmax
            break;
        end

        % choose reaction
        r2      = rand * a0;
        cumA    = cumsum(a);
        r_index = find(cumA >= r2, 1);

        % update state and time
        X = X + S(:, r_index);
        t = t_next;
    end

    % fill remaining samples at final state
    while idx <= scale_size
        traj(:,idx) = X;
        idx = idx + 1;
    end

    % parfor-sliced 쓰기 (한 번에)
    simuldataf(run,:,:) = traj;
end

%-------------------------------
% Compute statistics
%-------------------------------
mean_traj = squeeze(mean(simuldataf,1));   % species x time
var_traj  = squeeze(var(simuldataf,0,1));  % species x time
fano2     = var_traj(2,:) ./ mean_traj(2,:);% Fano factor for X2

writematrix(fano2,'Fano_NatComm.csv')

%% Figures

fano2 = readmatrix('Fano_NatComm.csv');

%-------------------------------
% Plot Fano factor for X2
%-------------------------------
figure;
plot(time_scale, fano2, 'LineWidth',1.5);
xlabel('Time'); ylabel('Fano factor');
title('Fano factor of X_2'); grid on;

