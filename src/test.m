clear all, close all; clc

data = dir("data\");
data = data(~ismember({data.name}, {'.', '..'})); % Filter /. and /.. subfolders

% Test data
data1 = load([data(1).name]);
samples = 1:data1.N; em_iter = 1:size(data1.A_hat, 3);

%% TEST CONNECTIVITY
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile(1), imshow(data1.A), title('A'), colormap("parula")
c = colorbar; c.Layout.Tile = 'south'; c.Label.String = 'Connection Strength';
for i = 1:67 %em_iter
    nexttile(2), imshow(data1.A_hat(:, :, i)), colormap("parula")
    title(['A hat ' num2str(i)])
    drawnow
end

% E-M estimated connection strength 1 \to 2
em_conn12 = squeeze(data1.A_hat(1, 2, :));
figure, plot(em_iter, em_conn12)
title('Estimated Connection Strength 1 to 2 Evolution')

%% TEST BOLD
reg = 1; % Time-series of the first region
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile(1), plot(samples, data1.y(:, reg)), title('y')
for i = 1:67 %em_iter
    nexttile(2), plot(data1.y_hat(:, reg, i))
    title(['y hat ' num2str(i)])
    drawnow
end

%% TEST HAEMODYNAMIC RESPONSE
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile(1), plot(data1.h)
for i = 1:67 %em_iter
    nexttile(2), plot(data1.h_hat(:, :, i))
    title(['h hat ' num2str(i)])
    drawnow
end
