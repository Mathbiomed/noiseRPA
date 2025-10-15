%%
control_idx = 2;
speci_num = 2;
x_list = linspace(2, 5, 16);
y_list = linspace(0.1, 2.5, 25);
% y_list = linspace(0.75, 1.5, 26);

target_mean = linspace(2, 5, 16);
sim_ff = zeros(1, length(target_mean));

%% Drawing Fig. 3a

figure('Position', [-100 -100 1000 800])
eta1_list = [5, 50, 500];
eta2_list = [5, 50, 500];

x_list = linspace(2, 5, 16);
y_list = linspace(0.1, 2.5, 25);

target_mean = linspace(2, 5, 16);
sim_ff = zeros(1, length(target_mean));

for eta1_idx = 1:3

    eta1 = eta1_list(eta1_idx);

    for k=1:length(target_mean)
        mn = target_mean(k);
        y = load(strcat("Fig3a_mn/eta1_",string(eta1),"_tar_mn_",string(mn),".mat"), "x");
        V_ss = squeeze(var(y.x,0,1));
        sim_ff(k) = mean(V_ss(1,800:1000)) / mn;
    end

    for eta2_idx = 1:3

        plot_loc_idx = eta1_idx + 3 * (eta2_idx - 1);
        subplot(3,3,plot_loc_idx);
        hold on

        eta1 = eta1_list(eta1_idx);
        eta2 = eta2_list(eta2_idx);

        load(strcat("Fig3a/heatmap_ff_eta1_",string(eta1),"_eta2_",string(eta2),".mat"), "heat_mat_ff")

         % Drawing a heatmap
        heat = log(heat_mat_ff)/log(10);
        h = imagesc(x_list, y_list, heat); 
        c = colorbar;
        c.Ticks = [0, 0.5, 1]; % Define the positions of the ticks
        c.TickLabels = {'10^0', '10^{0.5}', '10^1'}; % Custom labels

        targetColor = [0.3, 0.7, 0.3];
        customColormap = createMonochromaticColormap(targetColor, 256);
        colormap(customColormap);
        hold on;
        set(gca, 'YDir', 'normal')

        red_heat = [];
        for i = 1:length(x_list)
            red_heat = [red_heat; heat(y_list <= sim_ff(i), i)];
        end
        
        clim([0 1]);
        axis tight;
        plot(x_list, ones(length(x_list),1), 'r-', 'LineWidth',2);
        plot(x_list, sim_ff, 'b-', 'LineWidth',2);

        x_full = [fliplr(x_list), x_list]; % Use flipped x values to close the region
        y_top = 3 * ones(size(x_list)); % Bottom y-boundary to close the region (within heatmap)
        y_full = [y_top, sim_ff]; % Combine the y values to close the shaded region
        
        % Step 3: Plot the shaded region using `fill`
        fill(x_full, y_full, [0.5 0.5 0.5], 'FaceAlpha', 1); % Gray shaded region with transparency

        xlabel('Target mean')
        ylabel('Target Fano factor')
        xlim([2,5])
        ylim([y_list(1), y_list(end)])
        title('Fano factor ratio')
        hold off;

        fprintf(strcat("Finished eta1: ",string(eta1)," eta2: ",string(eta2),"\n"))
    end
end

%% Drawing Fig. 3b

figure('Position', [-100 -100 1000 800])
k1_list = [0.1, 1, 10];
d2_list = [0.1, 1, 10];

x_list = linspace(2, 5, 16);
y_list = linspace(0.1, 2.5, 25);

target_mean = linspace(2, 5, 16);
sim_ff = zeros(1, length(target_mean));

for k1_idx = 1:3

    k1 = k1_list(k1_idx);
    shade_lim = 3;
    if k1_idx == 1
        Ylim = [0.5, 2];
    elseif k1_idx == 2
        Ylim = [y_list(1), y_list(end)];
    else
        y_list = linspace(0.3, 7.5, 25);
        shade_lim = 10;
        Ylim = [y_list(1), y_list(end)];
    end

    for k=1:length(target_mean)
        mn = target_mean(k);
        y = load(strcat("Fig3b_mn/k1_",string(k1),"_tar_mn_",string(mn),".mat"), "x");
        V_ss = squeeze(var(y.x,0,1));
        sim_ff(k) = mean(V_ss(1,800:1000)) / mn;
    end

    for d2_idx = 1:3

        plot_loc_idx = k1_idx + 3 * (d2_idx - 1);
        subplot(3,3,plot_loc_idx);
        hold on

        k1 = k1_list(k1_idx);
        d2 = d2_list(d2_idx);

        load(strcat("Fig3b/heatmap_ff_k1_",string(k1),"_d2_",string(d2),".mat"), "heat_mat_ff")

         % Drawing a heatmap
        heat = log(heat_mat_ff)/log(10);
        h = imagesc(x_list, y_list, heat); 
        c = colorbar;
        c.Ticks = [0, 0.5, 1]; % Define the positions of the ticks
        c.TickLabels = {'10^0', '10^{0.5}', '10^1'}; % Custom labels

        targetColor = [0.3, 0.7, 0.3];
        customColormap = createMonochromaticColormap(targetColor, 256);
        colormap(customColormap);
        hold on;
        set(gca, 'YDir', 'normal')
        
        clim([0,1])
        % 
        axis tight;
        plot(x_list, ones(length(x_list),1), 'r-', 'LineWidth',2);
        plot(x_list, sim_ff, 'b-', 'LineWidth',2);

        x_full = [fliplr(x_list), x_list]; % Use flipped x values to close the region
        y_top = shade_lim * ones(size(x_list)); % Bottom y-boundary to close the region (within heatmap)
        y_full = [y_top, sim_ff]; % Combine the y values to close the shaded region
        
        % Step 3: Plot the shaded region using `fill`
        fill(x_full, y_full, [0.5 0.5 0.5], 'FaceAlpha', 1); % Gray shaded region with transparency

        xlabel('Target mean')
        ylabel('Target Fano factor')
        xlim([2,5])
        ylim(Ylim)
        title('Fano factor ratio')
        hold off;

        fprintf(strcat("Finished k1: ",string(k1)," d2: ",string(d2),"\n"))
    end
end

%% Drawing Fig. 3c

control_idx = 2;
speci_num = 2;

x_list = linspace(2, 5, 16);
y_list = linspace(0.1, 2.5, 25);

target_mean = linspace(2, 5, 16);
sim_ff = zeros(1, length(target_mean));

figure('Position', [-100 -100 1000 800])

for k=1:length(target_mean)
    mn = target_mean(k);
    y = load(strcat("Fig3c_mn/k_3_tar_mn_",string(mn),".mat"), "x");
    V_ss = squeeze(var(y.x,0,1));
    sim_ff(k) = mean(V_ss(control_idx,800:1000)) / mn;
end

load(strcat("Fig3c_heatmap_ff.mat"), "heat_mat_ff")

 % Drawing a heatmap
heat = log(heat_mat_ff)/log(10);
h = imagesc(x_list, y_list, heat); 
c = colorbar;
c.Ticks = [0, 0.5, 1]; % Define the positions of the ticks
c.TickLabels = {'10^0', '10^{0.5}', '10^1'}; % Custom labels

targetColor = [0.3, 0.7, 0.3];
customColormap = createMonochromaticColormap(targetColor, 256);
colormap(customColormap);
hold on;
set(gca, 'YDir', 'normal')

clim([0 1]);
axis tight;
plot(x_list, ones(length(x_list),1), 'r-', 'LineWidth',3);
plot(x_list, sim_ff, 'b-', 'LineWidth',3);

% Complete the region to close the shape (this is required for `fill`)
x_full = [fliplr(x_list), x_list]; % Use flipped x values to close the region
y_top = 3 * ones(size(x_list)); % Bottom y-boundary to close the region (within heatmap)
y_full = [y_top, sim_ff]; % Combine the y values to close the shaded region

% Step 3: Plot the shaded region using `fill`
fill(x_full, y_full, [0.5 0.5 0.5], 'FaceAlpha', 1); % Gray shaded region with transparency

xlabel('Target mean')
ylabel('Target Fano factor')
xlim([2,5])
ylim([y_list(1), y_list(end)])
% ylim([y_list(1), 2])
title('Fano factor ratio')
hold off;

%% Drawing Fig. 3d

control_idx = 1;
speci_num = 3;

x_list = linspace(2, 5, 16);
y_list = linspace(0.75, 1.5, 26);

target_mean = linspace(2, 5, 16);
sim_ff = zeros(1, length(target_mean));

figure('Position', [-100 -100 1000 800])

for k=1:length(target_mean)
    mn = target_mean(k);
    y = load(strcat("Fig3d_mn/k_2_tar_mn_",string(mn),".mat"), "x");
    V_ss = squeeze(var(y.x,0,1));
    sim_ff(k) = mean(V_ss(control_idx,800:1000)) / mn;
end

load(strcat("Fig3d_heatmap_ff.mat"), "heat_mat_ff")

 % Drawing a heatmap
heat = log(heat_mat_ff)/log(10);
h = imagesc(x_list, y_list, heat); 
c = colorbar;
c.Ticks = [0, 0.1, 0.2]; % Define the positions of the ticks
c.TickLabels = {'10^0', '10^{0.1}', '10^{0.2}'}; % Custom labels

targetColor = [0.3, 0.7, 0.3];
customColormap = createMonochromaticColormap(targetColor, 256);
colormap(customColormap);
hold on;
set(gca, 'YDir', 'normal')

clim([0 0.2]);

axis tight;
plot(x_list, ones(length(x_list),1), 'r-', 'LineWidth',3);
plot(x_list, sim_ff, 'b-', 'LineWidth',3);

% Complete the region to close the shape (this is required for `fill`)
x_full = [fliplr(x_list), x_list]; % Use flipped x values to close the region
y_top = 3 * ones(size(x_list)); % Bottom y-boundary to close the region (within heatmap)
y_full = [y_top, sim_ff]; % Combine the y values to close the shaded region

% Step 3: Plot the shaded region using `fill`
fill(x_full, y_full, [0.5 0.5 0.5], 'FaceAlpha', 1); % Gray shaded region with transparency

xlabel('Target mean')
ylabel('Target Fano factor')
xlim([2,5])
ylim([y_list(1), y_list(end)])
title('Fano factor ratio')
hold off;
%%
function customColormap = createMonochromaticColormap(color, numColors)
    % Create a colormap that transitions from black to the specified color
    customColormap = [linspace(1, color(1), numColors)', ...
                      linspace(1, color(2), numColors)', ...
                      linspace(1, color(3), numColors)'];
end