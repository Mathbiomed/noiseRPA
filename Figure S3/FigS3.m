%% initialization
clear
clc
rng('shuffle')
Tmax = 50000; %maximal simualtion time
Fs = 1;


period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time


simulnum = 50000; % number of simulations
folder = 'story1';
char = 'vc2_diffpert';
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
mu2 = 3;
theta2 = 1;
eta2 = 50;
d2 = 1;
idx1 = 3; % index of reaction which will be perturbed

mu2_list = [10.5, 9.9, 9.3, 8.7];
%%
for n = 1:4
    mu2 = mu2_list(n);
    
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
            rho = params .* kinetics;
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
    sc = squeeze(simuldataf(1,:,:));
    
    save(strcat("FigS3_mu2_",string(theta2),"_mean.mat"), 'X_ss')
    save(strcat("FigS3_mu2_",string(theta2),"_sc.mat"), 'sc')

end
%% mu2 mean plot
mu2_list = [10.5, 9.9, 9.3, 8.7];

figure()
for i = 1:4
    hold on
    mu_2 = mu2_list(i);
    load(strcat("FigS3_mu2_",string(mu_2),"_mean.mat"), 'X_ss')
    plot(time_scale,X_ss(5,:));
end
title('Z3 Mean')
xlabel('Time (s)')
ylabel('mean')
% yline(mu1/theta1, '--', LineWidth=2)
ylim([0 1])
% legend(strcat("mu2 = ", string(mu2_list)))
legend(strcat("FF = ", string([1.5 1.3 1.1 0.9])))

%% mu2 single cell
mu2_list = [10.5, 9.9, 9.3, 8.7];

for i=1:4
    subplot(2,2,i);
    mu_2 = mu2_list(i);
    load(strcat("FigS3_mu2_",string(mu_2),"_sc.mat"), 'sc')
    z3_sc = sc(5,:);
    plot(time_scale, z3_sc)
    % title(strcat("mu2 = ", string(mu_2)))
    title(strcat("Fano factor = ", string(mu_2/3 - 2)))
    xlabel('Time (s)')
    ylabel('Z_3 Copy number')
    xlim([0 50000])
    ylim([0, 10])
end
