%% initialization
clear
clc
rng('shuffle')
Tmax = 30000; %maximal simualtion time
Fs = 1;

period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time

simulnum = 1000; % number of simulations

%%
for what=["original", "mc", "nc"]
if what == "mc"
    filename = "telegraph_MC";
elseif what == "nc"
    filename = "telegraph_MC_NC";
else
    filename = "telegraph";
end
% This part changes the which controller you use
% "original": Original network without controllers
% "mc": Original network with MC
% "nc": Original network with MC and NC

%% Global parameters

% original network
global eps p1 p2 p3 p4 p5 p6 p7
eps = 10^(-3);
p1 = 1 * eps; % gene on 
p2 = 1 * eps; % gene off
p3 = 60; % mRNA production_on
p4 = 6; % mRNA production_off
p5 = 1; % mRNA degradation
p6 = 1; % protein production
p7 = 1; % protein degradation


% MC
global mu1 theta1 eta1 k1 

if (what == "mc") || (what == "nc")
    mu1 = 30;
    theta1 = 1;
    eta1 = 50;
    k1 = 1;
else
    mu1 = 0;
    theta1 = 0;
    eta1 = 0;
    k1 = 0;
end


% NC
global mu2 theta2 eta2 d2

if what == "nc"
    mu2 = 906;
    theta2 = 1;
    eta2 = 50;
    d2 = 1;
else
    mu2 = 0;
    theta2 = 0;
    eta2 = 0;
    d2 = 0;
end

idx1 = 10; % index of reaction which will be perturbed


%% full model simulation
% stoichiometric
%        r1 r2 r3 r4 r5 r6 r7 b1 m1 s1 a1 b1 m2 s2 d2
gamma = [ 1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0; %Gon
         -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0; %Goff
          0  0  1  1 -1  0  0  0  0  0  1  0  0  0 -1; %M
          0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0; %P
          0  0  0  0  0  0  0  1  0 -1  0  0  0  0  0; %z1
          0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0; %z2
          0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0; %z3
          0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0; %z4
         ];    
q = size(gamma);
reac_num = q(2); % number of reactions
var_num = q(1); % number of variables
simuldataf = zeros(simulnum,var_num,scale_size); %simulation record data full
simulparam = zeros(simulnum,scale_size);
prop_simuldataf = zeros(simulnum,reac_num,scale_size);
propensities = {};
%% Simulation
init = [1; 0; 0; 0; 0; 0; 0; 0];
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

        params = [p1; p2; p3; p4; p5; p6; p7; mu1; theta1; eta1; k1; mu2; theta2; eta2; d2];
        kinetics = [k(2); k(1); k(1); k(2); k(3); k(3); k(4); 1; k(4); k(5)*k(6); k(5); 1; k(4)*(k(4)-1); k(7)*k(8); k(8)*k(3)];

        if t < 700
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        elseif t < 1400
            params(idx1) = params(idx1);
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

save(filename, 'simuldataf')

end

%% Drawing figures
set(gcf, 'Renderer', 'painters')
figure('Units', 'inches', 'Position', [1, 1, 7, 7]);

for i = 1:9
    if rem(i,3) == 1
        filename = "telegraph";
    elseif rem(i,3) == 2
        filename = "telegraph_MC";
    else
        filename = "telegraph_MC_NC";
    end
    load(filename, 'simuldataf')

    subplot(3, 3, i);  

    if i<=3
        prod_during_on = zeros(30001, 1);
        prod_during_off = zeros(30001, 1);

        for t = 1:30001
            gene_on_cell = find(simuldataf(:,1,t) == 1);
            gene_off_cell = find(simuldataf(:,1,t) == 0);

            prod_during_on(t) = p5 * mean(simuldataf(gene_on_cell, 3, t));
            prod_during_off(t) = p5 * mean(simuldataf(gene_off_cell, 3, t));
        end

        hold on
        plot(time_scale, prod_during_on, 'r', time_scale, prod_during_off, "b");
        xlabel('Time (s)')
        ylabel('Prod. propensity')
        ylim([0,70])
        legend(["fast", "slow"])

    elseif i<=6
        plot(time_scale, squeeze(simuldataf(1,4,:)))
        xlabel('Time (s)')
        ylabel('Num. of P')
        yline(mu1/theta1, '--', LineWidth=1)
        ylim([0,100])
    else
        histogram(simuldataf(:,4,30001), 0:5:100, 'Normalization','pdf', 'EdgeColor','none', 'FaceAlpha',0.3)
        xlim([0,100])
        ylim([0 0.1])
        xline(mu1/theta1, '--', LineWidth=1)
        xlabel('The number of P')
        ylabel('Density')
    end

end
