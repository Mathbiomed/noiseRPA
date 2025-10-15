%% initialization
clear
clc
rng('shuffle')
Tmax = 10000; %maximal simualtion time
Fs = 0.1;

period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time


simulnum = 20; % number of simulations
folder = 'story1';
char = 'vc2_diffpert';
%% Global parameters
global p1 p2 p3 p4 p5 p6 mu1 mu2 theta1 theta2 eta1 eta2 k1 k2 idx1 idx2 d1 d2 idx3 idx4 idx5 idx6 idx7 idx8 idx9 idx10 idx11 idx12
p1 = 40.5;
p2 = 3;
p3 = 2; 
p4 = 1; 
p5 = 0.2;
p6 = 0.2;
mu1 = 27;
theta1 = 1;
eta1 = 50;
k1 = 1;
d1 = 0;

mu2 = (mu1/theta1)^2 + 100.05;
theta2 = 1;
eta2 = 50;
k2 = 0;
d2 = 1;

idx1 = 10; 
idx2 = 11;
idx3 = 12; 
idx4 = 13;
idx5 = 14;

idx6 = 9;
idx7 = 8;
idx8 = 7;
idx9 = 6;
idx10 = 5;
idx11 = 1;
idx12 = 3;


%% full model simulation
% stoichiometric
%        r1 r2 r3 r4 b1 m1 s1 a1 d1 b1 m2 s2 a2 d2
gamma = [ 1 -1  0  0  0  0  0  1 -1  0  0  0  1 -1; %x1
          0  0  1 -1  0  0  0  0  0  0  0  0  0  0; %x2
          0  0  0  0  1  0 -1  0  0  0  0  0  0  0; %z1
          0  0  0  0  0  1 -1  0  0  0  0  0  0  0; %z2
          0  0  0  0  0  0  0  0  0  1  0 -1  0  0; %z3
          0  0  0  0  0  0  0  0  0  0  1 -1  0  0; %z4
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
for i = 1:simulnum
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

        params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; d1; mu2; theta2; eta2; k2; d2];
        kinetics = [1; k(1); k(1); k(2); 1; k(2); k(3)*k(4); k(3); k(4)*k(1); 1; k(2)*(k(2)-1); k(5)*k(6); k(5); k(6)*k(1)];

        if t< 2000
            params(idx6) = params(idx6) * 0;
            params(idx7) = params(idx7) * 0;
            params(idx8) = params(idx8) * 0;
            params(idx9) = params(idx9) * 0;
            params(idx10) = params(idx10) * 0;
            params(idx1) = params(idx1) * 0;
            params(idx2) = params(idx2) * 0;
            params(idx3) = params(idx3) * 0;
            params(idx4) = params(idx4) * 0;
            params(idx5) = params(idx5) * 0;
            params(idx11) = params(idx11);
            params(idx12) = params(idx12);
            rho = params .* kinetics;
        elseif t < 4000
            params(idx1) = params(idx1) * 0;
            params(idx2) = params(idx2) * 0;
            params(idx3) = params(idx3) * 0;
            params(idx4) = params(idx4) * 0;
            params(idx5) = params(idx5) * 0;
            params(idx6) = params(idx6);
            params(idx7) = params(idx7);
            params(idx8) = params(idx8);
            params(idx9) = params(idx9);
            params(idx10) = params(idx10);
            params(idx11) = params(idx11) * 0;
            params(idx12) = params(idx12);
            rho = params .* kinetics;
        elseif t < 8000 
            params(idx1) = params(idx1);
            params(idx2) = params(idx2);
            params(idx3) = params(idx3);
            params(idx4) = params(idx4);
            params(idx5) = params(idx5);
            params(idx6) = params(idx6);
            params(idx7) = params(idx7);
            params(idx8) = params(idx8);
            params(idx9) = params(idx9);
            params(idx10) = params(idx10);
            params(idx11) = params(idx11) * 0;
            params(idx12) = params(idx12);
            rho = params .* kinetics;
        else
            params(idx1) = params(idx1) - 100;
            params(idx2) = params(idx2);
            params(idx3) = params(idx3);
            params(idx4) = params(idx4);
            params(idx5) = params(idx5);
            params(idx6) = params(idx6);
            params(idx7) = params(idx7);
            params(idx8) = params(idx8);
            params(idx9) = params(idx9);
            params(idx10) = params(idx10);
            params(idx11) = params(idx11) * 0;
            params(idx12) = params(idx12);
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
%    figure()
%    hold on
%    plot(time_scale,X(2,:));
%    title('Copy Number')
%    xlabel('t')
%    ylabel('mean')
%    yline(mu1/theta1, '--', LineWidth=2)
end

X_ss = squeeze(mean(simuldataf,1));
V_ss = squeeze(var(simuldataf,0,1));
par_ss = squeeze(mean(simulparam,1));
%% plot


figure()
hold on
for i = 1:20
    plot(time_scale,squeeze(simuldataf(i,2,:))); hold on
end
title('Mean')
xlabel('t')
ylabel('mean')
yline(mu1/theta1, '--', LineWidth=2)
