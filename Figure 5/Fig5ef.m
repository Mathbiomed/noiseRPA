%% initialization
clear
clc
rng('shuffle')
Tmax = 10005; %maximal simualtion time
Fs = 1;
period = 1/Fs;
time_scale = 0:period:Tmax; % record time scale
scale_size = length(time_scale); % number of record time
simulnum = 500; % number of simulations

%% Global parameters
global k1 k2 k3 k4 k5 k6 mu1 theta1 eta1 eta2 theta2 mu2 d1 d2 r1 r2 k7
k1 = 7.5;
k3 = 8;
k4 = 1;
k5 = 0.2;
k6 = 10;
k7 = 1;
d1 = 1;
mu1 = 1.5 * 3;
theta1 = 1;
eta1 = 50;
r1 = 1;
%d1 = 0;
%mu1 = 0;
%theta1 = 0;
%eta1 = 0;
%r1 = 0;
r2 = 1;
theta2 = 1;
mu2 = (mu1/theta1)*(mu1/theta1) + 0.05;
eta2 = 50;
d2 = 1;
%r2 = 0;
%theta2 = 0;
%mu2 = theta2 * 0;
%eta2 = 0;
%d2 = 0;
%% full model simulation
% stoichiometric
gamma = [-1 1 0 0 0 0 0 0 0 0 0 0 0  -1 0 0; %Ada
          0 0 0 0 0 1 0 0 0 -1 1 -1 0 0 0 1;%mRNA
          0 0 1 0 -1 0 0 0 0 0 0 0 0 0 0 0; %Z1
          0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0; %Z2
          0 0 0 0 0 0 1 0 -1 0 0 0 0 0 0 0; %Z3
          0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0; %Z4
          0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0; %MMS
          0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0; %me-Ada
];
q = size(gamma);
reac_num = q(2); % number of reactions
var_num = q(1); % number of variables
simuldataf = zeros(simulnum,var_num,scale_size); %simulation record data full0
%% Simulation
for i = 1:simulnum
fprintf('Simulation %d\n',i);
alpha = 10;
theta = 1.5/alpha;
k2 = gamrnd(alpha,theta);
%k2 = normrnd(1.5, 0.5)*(5/1.6);
if k2 <= 0
k2 = 0.01;
end
t = 0; %current time
k = [1; 1; 1; 1; 1; 1; 0; 0]; %Initial Value
X = zeros(var_num, scale_size); %X(t) record
j = 1; %iterator
MMScheck = 0;
%% iteration
while t <= Tmax
%Omega = 1
if t>10000 & MMScheck == 0
k(7) = k(7) + 100;
MMScheck = 1;
end
rho = [k1*k(1); k2*k(2); mu1; theta1*k(1); eta1*k(3)*k(4); d1*k(3); mu2; theta2*k(1)*(k(1)-1); eta2*k(5)*k(6); d2*k(6)*k(2); k3; k4*k(2); k5*0.01*k(7); k6*k(7)*((k(1)*k(1))/(1 + k(1)*k(1))); k1*0.4*k(8); k7*k(8)*(k(8)-1)];
for m = 1:size(gamma, 2)
if rho(m) < 0
rho(m) = 0;
end
end
lambda = sum(rho);
r = rand([2 1]); %two random numbers r1, r2
T = -1/lambda * log(r(1));
if t + T > Tmax %end condition
while j <= scale_size
X(:,j) = k;
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
j = j + 1;
end
%update k and t
k = k + gamma(:, reaction_index);
t = t + T;
end
simuldataf(i,:,:) = X;

end
X_ss = squeeze(mean(simuldataf,1));
V_ss = squeeze(var(simuldataf,0,1));

%%
figure()
histogram(simuldataf(:,1,10005)+simuldataf(:,8,10005), 'Normalization','probability', 'BinWidth', 900);
%%

writematrix(simuldataf(:,1,10005)+simuldataf(:,8,10005),'MMC_original.csv') 
