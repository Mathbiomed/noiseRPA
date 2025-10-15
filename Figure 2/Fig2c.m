clear
clc
rng('shuffle')

simuldir = 'Fig2c';
mkdir(simuldir)

target_mean = linspace(2, 5, 16);
target_fanofactor = linspace(0.1, 2.5, 25);

control_idx = 2;
speci_num = 2;

x_list = target_mean;
y_list = target_fanofactor;

heat_mat_ff = zeros(length(y_list), length(x_list));

for i = 1:length(target_mean)
    for j = 1:length(target_fanofactor)
        mn = target_mean(i);
        ff = target_fanofactor(j);
        vr = mn * ff;
                                                        %k1 k2 d1 d2
        x = stochastic_simulation_1000(mn, vr+mn^2-mn, 50, 50, 1, 0, 0, 1);

        X_ss = squeeze(mean(x,1));
        V_ss = squeeze(var(x,0,1));
        heat_mat_ff(j,i) = mean(V_ss(control_idx,8000:10000))/ (ff * mean(X_ss(control_idx,8000:10000)));

        fprintf([repmat('-',1,20),'%d %d', repmat('-',1,20),'\n'], i, j )
    end
end

save(strcat("Fig2c_heatmap_ff.mat"), "heat_mat_ff")