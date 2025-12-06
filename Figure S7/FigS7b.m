%% initialization
clear
clc
rng('shuffle')
Tmax = 10000; %maximal simualtion time
Fs = 5;


period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time


simulnum = 50; % number of simulations
folder = 'story1';
char = 'vc2_diffpert';
%% Global parameters
global p1 p2 p3 p4 p5 p6 mu1 mu2 theta1 theta2 eta1 eta2 k1 k2 idx1 idx2 d1 d2
p1 = 2 * 10;
p2 = 1* 10;
%p3 = 4; 
p3 = 0* 10;
p4 = 3* 10; 
p5 = 0.2* 10;
p6 = 0.2* 10;
mu1 = 3* 10;
theta1 = 0.5* 10;
eta1 = 50* 10;
k1 = 1* 10;
d1 = 0* 10;
%mu1 = 0;
%theta1 = 0;
%eta1 = 0;
%k1 = 0;
%d1 = 0;
mu2 = 13* 10;
theta2 = 1* 10;
eta2 = 50* 10;
k2 = 0* 10;
d2 = 1* 10;

%mu2 = 0;
%theta2 = 0;
%eta2 = 0;
%k2 = 0;
%d2 = 0;
idx1 = 1; % index of reaction which will be perturbed
idx2 = 1;


%% full model simulation
% stoichiometric
%        r1 r2 r3 r4 b1 m1 s1 a1 d1 b1 m2 s2 a2 d2
gamma = [ 1 -1  0  0  0  0  0  0  0  0  0  0  0  0; %x1
          0  0  1 -1  0  0  0  1 -1  0  0  0  1 -1; %x2
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
        %Omega = 1
        params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; d1; mu2; 2*theta2; eta2; k2; d2];
        kinetics = [k(2); k(1); 1; k(2); 1; k(1); k(3)*k(4); k(3); k(4)*k(2); 1; (k(1)*k(1)*25*25)/((k(1)*k(1))+25*25); k(5)*k(6); k(5); k(6)*k(2)];
        % if k(5) > 350
        %     params(end-1) = 1;
        %     params(end) = 0;
        % end
        if t < 3000
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        else
            params(idx1) = params(idx1)*8;
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
%    plot(time_scale,X(1,:));
%    title('Copy Number')
%    xlabel('t')
%    ylabel('mean')
%    yline(mu1/theta1, '--', LineWidth=2)
end

X_ss1 = squeeze(mean(simuldataf,1));
V_ss1 = squeeze(var(simuldataf,0,1));
par_ss1 = squeeze(mean(simulparam,1));
%% plot



figure()
hold on
plot(time_scale,V_ss1(1,:)./X_ss1(1,:));
yline(((mu2/theta2)-(mu1/theta1)^2)/(mu1/theta1), '--', LineWidth=2)
title('Fano Factor')
xlabel('t')

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
        %Omega = 1
        params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; d1; mu2; 2*theta2; eta2; k2; d2];
        kinetics = [k(2); k(1); 1; k(2); 1; k(1); k(3)*k(4); k(3); k(4)*k(2); 1; (k(1)*k(1)*50*50)/((k(1)*k(1))+50*50); k(5)*k(6); k(5); k(6)*k(2)];
        % if k(5) > 350
        %     params(end-1) = 1;
        %     params(end) = 0;
        % end
        if t < 3000
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        else
            params(idx1) = params(idx1)*8;
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
%    plot(time_scale,X(1,:));
%    title('Copy Number')
%    xlabel('t')
%    ylabel('mean')
%    yline(mu1/theta1, '--', LineWidth=2)
end

X_ss2 = squeeze(mean(simuldataf,1));
V_ss2 = squeeze(var(simuldataf,0,1));
%% plot



figure()
hold on
plot(time_scale,V_ss2(1,:)./X_ss2(1,:));
yline(((mu2/theta2)-(mu1/theta1)^2)/(mu1/theta1), '--', LineWidth=2)
title('Fano Factor')
xlabel('t')

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
        %Omega = 1
        params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; d1; mu2; 2*theta2; eta2; k2; d2];
        kinetics = [k(2); k(1); 1; k(2); 1; k(1); k(3)*k(4); k(3); k(4)*k(2); 1; (k(1)*k(1)*75*75)/((k(1)*k(1))+75*75); k(5)*k(6); k(5); k(6)*k(2)];
        % if k(5) > 350
        %     params(end-1) = 1;
        %     params(end) = 0;
        % end
        if t < 3000
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        else
            params(idx1) = params(idx1)*8;
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
%    plot(time_scale,X(1,:));
%    title('Copy Number')
%    xlabel('t')
%    ylabel('mean')
%    yline(mu1/theta1, '--', LineWidth=2)
end

X_ss3 = squeeze(mean(simuldataf,1));
V_ss3 = squeeze(var(simuldataf,0,1));
%% plot



figure()
hold on
plot(time_scale,V_ss3(1,:)./X_ss3(1,:));
yline(((mu2/theta2)-(mu1/theta1)^2)/(mu1/theta1), '--', LineWidth=2)
title('Fano Factor')
xlabel('t')

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
        %Omega = 1
        params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; d1; mu2; 2*theta2; eta2; k2; d2];
        kinetics = [k(2); k(1); 1; k(2); 1; k(1); k(3)*k(4); k(3); k(4)*k(2); 1; (k(1)*k(1)*100*100)/((k(1)*k(1))+100*100); k(5)*k(6); k(5); k(6)*k(2)];
        % if k(5) > 350
        %     params(end-1) = 1;
        %     params(end) = 0;
        % end
        if t < 3000
            params(idx1) = params(idx1);
            rho = params .* kinetics;
        else
            params(idx1) = params(idx1)*8;
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
%    plot(time_scale,X(1,:));
%    title('Copy Number')
%    xlabel('t')
%    ylabel('mean')
%    yline(mu1/theta1, '--', LineWidth=2)
end

X_ss4 = squeeze(mean(simuldataf,1));
V_ss4 = squeeze(var(simuldataf,0,1));
%% plot



figure()
hold on
plot(time_scale,V_ss4(1,:)./X_ss4(1,:));
yline(((mu2/theta2)-(mu1/theta1)^2)/(mu1/theta1), '--', LineWidth=2)
title('Fano Factor')
xlabel('t')

figure()
hold on
plot(time_scale,V_ss1(1,:)./X_ss1(1,:)); hold on
plot(time_scale,V_ss2(1,:)./X_ss2(1,:)); hold on
plot(time_scale,V_ss3(1,:)./X_ss3(1,:)); hold on
plot(time_scale,V_ss4(1,:)./X_ss4(1,:)); hold on
yline(((mu2/theta2)-(mu1/theta1)^2)/(mu1/theta1), '--', LineWidth=2)
legend('25', '50', '75', '100')
title('Fano Factor')
xlabel('t')

writematrix(V_ss1(1,:)./X_ss1(1,:),'Fano_1.csv') 
writematrix(V_ss2(1,:)./X_ss2(1,:),'Fano_2.csv') 
writematrix(V_ss3(1,:)./X_ss3(1,:),'Fano_3.csv') 
writematrix(V_ss4(1,:)./X_ss4(1,:),'Fano_4.csv') 