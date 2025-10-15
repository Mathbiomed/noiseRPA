function simuldataf = stochastic_simulation_1000(mu1,mu2,eta1, eta2,k1,k2, d1, d2)
    rng('shuffle')
    Tmax = 10000; %maximal simualtion time
    Fs = 1;
    
    period = 1/Fs;
    time_scale = 0:period:Tmax; % record time scale
    scale_size = length(time_scale); % number of record time
    
    simulnum = 1000; % number of simulations
    %% Global parameters
    global p1 p2 p3 p4 theta1 theta2
    p1 = 0;
    p2 = 3;
    p3 = 2; 
    p4 = 1; 

    theta1 = 1;
    %eta1 = 50;

    theta2 = 1;
    %eta2 = 50;
    
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
        while t <= Tmax
            %Omega = 1
            params = [p1; p2; p3; p4; mu1; theta1; eta1; k1; d1; mu2; theta2; eta2; k2; d2];
            kinetics = [1; k(1); k(1); k(2); 1; k(2); k(3)*k(4); k(3); k(4)*k(1); 1; k(2)*(k(2)-1); k(5)*k(6); k(5); k(6)*k(1)];
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
    end

end
