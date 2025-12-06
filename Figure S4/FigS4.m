%% initialization
clear
clc
rng('shuffle')
Tmax = 170000; %maximal simualtion time
Fs = 1;


period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time


simulnum = 20000; % number of simulations

%% Global parameters
global p1 p2 p3 p4 mu1 mu2 theta1 theta2 eta1 eta2 k1 k2 idx1 idx2 d1 d2
% p1 = 4.5;
p1 = 0;
p2 = 1;
p3 = 2;
p4 = 3; 
mu1 = 3;
theta1 = 1;
eta1 = 50;
k1 = 1;
mu2 = ((mu1/theta1)^2 + 0.05);
theta2 = 1;
eta2 = 50;
d2 = 1;
idx1 = 3; % index of reaction which will be perturbed

per_list = [20, 50, 100];
%%
for n = 1:3
    per = per_list(n);
    
    %% full model simulation
    % stoichiometric
    %        r1 r2 r3 r4 b1 m1 s1 a1 b1 m2 s2 d2
    gamma = [ 1 -1  0  0  0  0  0  1  0  0  0 -1; %x1
              0  0  1 -1  0  0  0  0  0  0  0  0; %x2
              0  0  0  0  1  0 -1  0  0  0  0  0; %z1
              0  0  0  0  0  1 -1  0  0  0  0  0; %z2
              0  0  0  0  0  0  0  0  1  0 -1  0; %z3
              0  0  0  0  0  0  0  0  0  1 -1  0; %z4
             ];    
    q = size(gamma);
    reac_num = q(2); % number of reactions
    var_num = q(1); % number of variables
    simuldataf = zeros(simulnum,var_num,scale_size); %simulation record data full
    simulparam = zeros(simulnum,scale_size);
    prop_simuldataf = zeros(simulnum,reac_num,scale_size);
    propensities = {};
    %% Simulation
    init = [0; 0; 0; 0; 0; 0];
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
    
            params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; mu2; theta2; eta2; d2];
            kinetics = [1; k(1); k(1); k(2); 1; k(2); k(3)*k(4); k(3); 1; k(2)*(k(2)-1); k(5)*k(6); k(6)*k(1)];
    
            if t < 100000
                params(idx1) = params(idx1);
                rho = params .* kinetics;
            else
                params(idx1) = params(idx1) * per;
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
        %% Uncomment this block to draw Fig 1g
    %    figure()
    %    hold on
    %    plot(time_scale,X(1,:));
    %    title('Copy Number')
    %    xlabel('t')
    %    ylabel('mean')
    %    yline(mu1/theta1, '--', LineWidth=2)
    end
    
    X_ss = squeeze(mean(simuldataf,1));
    V_ss = squeeze(var(simuldataf,0,1));
    par_ss = squeeze(mean(simulparam,1));
    
    save(strcat("FigS4_per_",string(per),"_mean_many.mat"), 'X_ss')
    save(strcat("FigS4_per_",string(per),"_variance_many.mat"), 'V_ss')
end
%% plot
per = 50;
load(strcat("FigS4_per_",string(per),"_mean.mat"), 'X_ss')
load(strcat("FigS4_per_",string(per),"_variance.mat"), 'V_ss')

if per == 10000
    time_scale = 0:1:2000000-70000;
else
    time_scale = 0:1:100000;
end

figure('Units', 'inches','Position', [1 1 10 3])
hold on
plot(time_scale,X_ss(2,70001:end));
title('Mean')
xlabel('Time (s)')
ylabel('X_2 mean')
yline(mu1/theta1, '--', LineWidth=2)
xlim([0,time_scale(end)])
ylim([0,10])
legend('x1')

figure('Units', 'inches','Position', [1 1 10 3])
hold on
plot(time_scale,V_ss(2,70001:end)./X_ss(2,70001:end));
title('X_2 Fano Factor')
yline(1 - mu1/theta1 + mu2/theta2 * (theta1/mu1), '--', LineWidth=2)
xlim([0,time_scale(end)])
ylim([0,2])
ylabel('X_2 Fano factor')
xlabel('Time (s)')
