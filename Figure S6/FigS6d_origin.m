% Gillespie simulation for APIF Class 1 with Degradation Inhibition

clear; close all; clc;
rng('shuffle');

%-------------------------------
% Simulation time grid
%-------------------------------
Tmax       = 1e5;    % maximal simulation time
Fs         = 0.1;        % sampling frequency
period     = 1/Fs;
time_scale = 0:period:Tmax;
scale_size = length(time_scale);
simulnum   = 10000;       % number of independent runs

%-------------------------------
% Model parameters
%-------------------------------

P.kp       = 1; 
P.gamma_r = 2;%%
P.gamma_p = 1; %%
P.kd = 3; % degradation rate
P.gammad = 1; % degradation rate for species Z1
P.gammadr = 1; % degradation rate for species Z2
P.kr = 0; % reaction rate for R1
P.k = 0.25;

P.mu      = 10; %%
P.theta   = 2; %%
P.k_1     = 2; %%
P.eta     = 100;    %% % degradation inhibition rate
P.kkp = 5;

%P.mu      = 0;
%P.theta   = 0;
%P.k_1     = 2;
%P.eta     = 0;     % degradation inhibition rate
%P.kappa_1 = 1;

%-------------------------------
% Stoichiometry matrix
% species: [X1; X2; Z1; Z2]
% reactions R1...R10
S = [ ...
     1,   0,  -1,   0,   0,   0,   0,   0, 0, 0, 1;   % X1
     0,   1,   0,  -1,  -2,   2,   0,   0, 0, 0, 0;   % X2
     0,   0,   0,   0,   1,  -1,  -1,   0, 0, 0, 0;%X3
     0,   0,   0,   0,   0,   0,   0,   1, 0, -1, 0;   % Z1
     0,   0,   0,   0,   0,   0,   0,   0, 1, -1, 0    % Z2
];
var_num  = size(S,1);
reac_num = size(S,2);

% Preallocate storage: runs x species x time
simuldataf = zeros(simulnum, var_num, scale_size);

%-------------------------------
% Run Gillespie-sampled trajectories
%-------------------------------
for run = 1:simulnum
    fprintf("%d\n", run);
    t = 0;
    %X = zeros(var_num,1);   % initial state [0;0;0;0]
    X = [10; 10; 10; 10; 10];
    idx = 1;

    while t <= Tmax
        % compute propensities
        a = [ ...
            P.k * X(4);                          % R1
            P.kp * X(1);                     % R2
            P.gamma_r * X(1);
            P.gamma_p * X(2);% R3
            P.kd * X(2) * (X(2)-1);
            P.gammad * X(3);
            P.gammadr * X(3);
            P.mu;                 % R4
            P.theta * X(3);                       % R5
            P.eta * X(4) * X(5);                  % R6
            P.kkp   * max([0, P.mu-P.theta*X(3)]); 
            % R7
             % R8                    % R10
        ];

        if t > Tmax/2
            a = [ ...
            P.k * 10 * X(4);                          % R1
            P.kp * X(1);                     % R2
            P.gamma_r * X(1);
            P.gamma_p * X(2);% R3
            P.kd * X(2) * (X(2)-1);
            P.gammad * X(3);
            P.gammadr * 10 * X(3);
            P.mu;                 % R4
            P.theta * X(3);                       % R5
            P.eta * X(4) * X(5);                  % R6
            P.kkp   * max([0, P.mu-P.theta*X(3)]); 
            % R7
             % R8                    % R10
        ];
        end

        
        

        

        a0 = sum(a);
        if a0 <= 0
            break;
        end
        % time to next reaction
        tau = -log(rand)/a0;
        t_next = t + tau;

        % record state until next event
        while idx <= scale_size && time_scale(idx) <= t_next
            simuldataf(run,:,idx) = X';
            idx = idx + 1;
        end
        if t_next > Tmax
            break;
        end
        % choose reaction
        r2 = rand * a0;
        cumA = cumsum(a);
        reaction_index = find(cumA >= r2,1);

        % update state and time
        X = X + S(:,reaction_index);
        t = t_next;
    end
    % fill remaining samples at final state
    while idx <= scale_size
        simuldataf(run,:,idx) = X';
        idx = idx + 1;
    end
end

%-------------------------------
% Compute statistics
%-------------------------------
mean_traj = squeeze(mean(simuldataf,1));   % species x time
var_traj  = squeeze(var(simuldataf,0,1));  % species x time
fano2     = var_traj(3,:) ./ mean_traj(3,:);% Fano factor for X2


writematrix(fano2,'Fano_Royal.csv') 

%% Figures 

fano2 = readmatrix('Fano_Royal.csv');

%-------------------------------
% Plot Fano factor for X2
%-------------------------------
figure;
plot(time_scale, fano2, 'LineWidth',1.5);
xlabel('Time'); ylabel('Fano factor');
title('Fano factor of X_2'); grid on;