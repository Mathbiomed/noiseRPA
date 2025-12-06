%% initialization
clear
clc
rng('shuffle')
Tmax = 15000; %maximal simualtion time
Fs = 1;


period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time


simulnum = 10000; % number of simulations

%% Global parameters

% original network
global p1 p2 p3 p4
p1 = 0;
p2 = 3;
p3 = 2;
p4 = 1; 

% MC
global mu1 theta1 eta1 k1 
mu1 = 3;
theta1 = 1;
eta1 = 50;
k1 = 1;

% NC
global mu2 theta2 eta2 d2
mu2 = ((mu1/theta1)^2 + 0.6);
theta2 = 1;
eta2 = 50;
d2 = 1;

% Anti-windup
global g1 vm eta_v d_v1 d_v2 h1
g1 = 1; vm = 10; eta_v = 100; d_v1 = 1; d_v2 = 1; h1 = 1;

global g2 wm eta_w d_w1 d_w2 h2
g2 = 1; wm = 20; eta_w = 100; d_w1 = 1; d_w2 = 1; h2 = 1;

global g3 vn d_v3 d_v4 h3
g3 = 1; vn = 10; d_v3 = 1; d_v4 = 1; h3 = 1;

global g4 wn d_w3 d_w4 h4
g4 = 1; wn = 20; d_w3 = 1; d_w4 = 1; h4 = 1;

%           r1 r2 r3 r4 b1 m1 s1 a1 b1 m2 s2 d2
gamma_nc = [ 1 -1  0  0  0  0  0  1  0  0  0 -1; %x1
             0  0  1 -1  0  0  0  0  0  0  0  0; %x2
             0  0  0  0  1  0 -1  0  0  0  0  0; %z1
             0  0  0  0  0  1 -1  0  0  0  0  0; %z2
             0  0  0  0  0  0  0  0  1  0 -1  0; %z3
             0  0  0  0  0  0  0  0  0  1 -1  0; %z4
            ];

gamma_aw = [0  1 -1 -1  0  0  0  0  0  0  0  0;
            1  0 -1  0 -1  0  0  0  0  0  0  0; 
            0  0  0  0  0  0  0  1 -1 -1  0  0;
            0  0  0  0  0  0  1  0 -1  0 -1  0;
            ];

gamma_aw0 = zeros(6, 24);
gamma_aw0(3,6) = 1; gamma_aw0(4,12) = 1; gamma_aw0(5,18) = 1; gamma_aw0(6,24) = 1;

gamma = [gamma_nc, gamma_aw0; zeros(4, size(gamma_nc, 2)), gamma_aw, zeros(4, 12); zeros(4, size(gamma_nc, 2)), zeros(4, 12), gamma_aw];

idx1 = 9; % index of reaction which will be perturbed


%% full model simulation
% stoichiometric

q = size(gamma);
reac_num = q(2); % number of reactions
var_num = q(1); % number of variables
simuldataf = zeros(simulnum,var_num,scale_size); %simulation record data full
simulparam = zeros(simulnum,scale_size);
prop_simuldataf = zeros(simulnum,reac_num,scale_size);
propensities = {};
%% Simulation
init = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
tic
parfor i = 1:simulnum
    fprintf('Simulation %d\n',i);
    t = 0; %current time
    k = init; %Initial Value
    X = zeros(var_num, scale_size); %X(t) record
    P = zeros(reac_num, scale_size); %
    par = zeros(1,scale_size);
    event_time = [];
    P_real = [];
    j = 1; %iterator
    %% iteration
    rho2 = [];
    while t <= Tmax

        params_nc = [p1; p2; p3; p4; mu1; theta1; eta1; k1; mu2; theta2; eta2; d2];
        params_aw = [g1; vm; eta_v; d_v1; d_v2; h1; g2; wm; eta_w; d_w1; d_w2; h2];
        params_aw = [params_aw; g3; vn; eta_v; d_v3; d_v4; h3; g4; wn; eta_w; d_w3; d_w4; h4];
        params = [params_nc; params_aw];

        kinetics_nc = [1; k(1); k(1); k(2); 1; k(2); k(3)*k(4); k(3); 1; k(2)*(k(2)-1); k(5)*k(6); k(6)*k(1)];
        kinetics_aw = [k(4); 1; k(7)*k(8); k(7); k(8); k(8); k(3); 1; k(9)*k(10); k(9); k(10); k(10)];
        kinetics_aw = [kinetics_aw; k(6); 1; k(11)*k(12); k(11); k(12); k(12); k(5); 1; k(13)*k(14); k(13); k(14); k(14)];
        kinetics = [kinetics_nc; kinetics_aw];

        if t < 5000
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        elseif t < 10000
            params(idx1) = 7.5;
            rho = params .* kinetics;
        else
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        end

        lambda = sum(rho);
        
        r = rand([2 1]); %two random numbers r1, r2
        T = -1/lambda * log(r(1));
        
        if t + T > Tmax %end condition
            while j <= scale_size
                X(:,j) = k;
                P(:,j) = rho;
                j = j + 1;
            end
            break
        end
        
        %choose the reaction
        rho_sum = 0;
        for l = 1:reac_num
            rho_sum = rho_sum + rho(l);
            if r(2) * lambda < rho_sum
                reaction_index = l;
                break
            end
        end
        
        %record the X(t)
        while time_scale(1,j) < t + T
            X(:,j) = k;
            par(j) = params(end);
            j = j + 1;
        end
        
        %update k and t
        k = k + gamma(:, reaction_index);
        t = t + T;
    end
    simuldataf(i,:,:) = X;
    simulparam(i,:) = par;
end

X_ss = squeeze(mean(simuldataf,1));
V_ss = squeeze(var(simuldataf,0,1));

save("FigS2_AW_topo1_mean_10000.mat", 'X_ss')
save("FigS2_AW_topo1_variance_10000.mat", 'V_ss')

%% Fig S2a, b

load("FigS2_noAW_mean.mat", 'X_ss')
load("FigS2_noAW_variance.mat", 'V_ss')

figure()
hold on
plot(time_scale, X_ss(3,:)); %index of X_ss = 3(Z1) or 6(Z4)
title('Target Mean')
xlabel('t')
ylabel('mean')

FF_set = repelem(mu2/theta2 * 1/(mu1/theta1) + 1 - mu1/theta1, length(time_scale));
FF_set((time_scale > 5000) & (time_scale < 10000)) = 0.8;

figure()
hold on
plot(time_scale,V_ss(2,:)./X_ss(2,:));
title('Target Fano Factor')
plot(time_scale, FF_set, 'k--', LineWidth=1)
xlabel('t')
ylim([0.6,2])

%% Fig S2c, d

load("FigS2_AW_mean.mat", 'X_ss')
load("FigS2_AW_variance.mat", 'V_ss')

figure()
hold on
plot(time_scale, X_ss(3,:)); %index of X_ss = 3(Z1) or 6(Z4)
title('Target Mean')
xlabel('t')
ylabel('mean')
yline(mu1/theta1, '--', LineWidth=2)

FF_set = repelem(mu2/theta2 * 1/(mu1/theta1) + 1 - mu1/theta1, length(time_scale));
FF_set((time_scale > 5000) & (time_scale < 10000)) = 0.8;

figure()
hold on
plot(time_scale,V_ss(2,:)./X_ss(2,:));
title('Target Fano Factor')
plot(time_scale, FF_set, 'k--', LineWidth=1)
xlabel('t')
ylim([0.6,2])