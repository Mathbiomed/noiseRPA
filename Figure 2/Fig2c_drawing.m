%% Drawing Fig. 2c

control_idx = 2;
speci_num = 2;

x_list = linspace(2, 5, 16);
y_list = linspace(0.1, 2.5, 25);

target_mean = linspace(2, 5, 16);
sim_ff = zeros(1, length(target_mean));

figure('Position', [-100 -100 1000 800])

for k=1:length(target_mean)
    mn = target_mean(k);
    y = load(strcat("Fig2c_mn/eta1_50_tar_mn_",string(mn),".mat"), "x");
    V_ss = squeeze(var(y.x,0,1));
    sim_ff(k) = mean(V_ss(1,800:1000)) / mn;
end

load(strcat("Fig2c_heatmap_ff.mat"), "heat_mat_ff")

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

%%
function customColormap = createMonochromaticColormap(color, numColors)
    % Create a colormap that transitions from black to the specified color
    customColormap = [linspace(1, color(1), numColors)', ...
                      linspace(1, color(2), numColors)', ...
                      linspace(1, color(3), numColors)'];
end