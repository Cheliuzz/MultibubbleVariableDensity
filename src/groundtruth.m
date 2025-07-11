%% --- 1. Ground Truth Pursuit Bouts (in frames) ---
ground_truth_bouts = [
    72000, 73419;
    73967, 74087;
    74443, 79525;
    80086, 80283;
    80438, 81098;
    81098, 82544;
    82544, 82768;
    82991, 85073;
    85336, 85430;
    87257, 87511;
    88034, 90102;
    90210, 92899;
    92968, 94160;
    94498, 94846;
    95095, 99318;
    99490, 99591;
    100026, 100453;
    100509, 101143;
    101211, 103597;
    103793, 107717;
    107717, 109095;
    109163, 110194;
    110707, 110969;
    111403, 114626;
    114668, 114773;
    115766, 116732;
    117467, 118177;
    118315, 118449;
    119016, 128725;
    128866, 131809;
    131976, 132700;
    132933, 139193;
    139535, 140515;
];

start_frame = 72000;
end_frame = 140515;
n_window_frames = end_frame - start_frame + 1;

% Initialize ethogram
ground_truth_ethogram = false(n_window_frames, 1);
for i = 1:size(ground_truth_bouts, 1)
    s = max(ground_truth_bouts(i, 1) - start_frame + 1, 1);
    e = min(ground_truth_bouts(i, 2) - start_frame + 1, n_window_frames);
    ground_truth_ethogram(s:e) = true;
end

%% --- 2. Compare to .fig Ethograms ---
fig_prefix = 'Figure_';
fig_indices = 401:412;
FPS = 60;

match_scores = zeros(length(fig_indices), 1);
ethogram_data_all = cell(length(fig_indices), 1);

for i = 1:length(fig_indices)
    fig_name = sprintf('%s%d.fig', fig_prefix, fig_indices(i));
    fig = openfig(fig_name, 'invisible');

    ax = findobj(fig, 'Type', 'Axes');
    lines = findobj(ax, 'Type', 'Line');

    if isempty(lines)
        fprintf('Figure_%d: No lines found.\n', fig_indices(i));
        continue;
    end

    % Assume ethogram is the first or longest line
    [~, idx_longest] = max(arrayfun(@(l) length(l.YData), lines));
    ethogram_line = lines(idx_longest);
    ethogram_y = ethogram_line.YData;
    ethogram_x = ethogram_line.XData;

    % Heuristic: check XData units and convert if needed
    if max(ethogram_x) <= 60*60 && max(ethogram_x) > 1000
        % Probably in seconds
        ethogram_x = round(ethogram_x * FPS);
    elseif all(diff(ethogram_x) == 1) && length(ethogram_x) == length(ethogram_y)
        % Just indices
        ethogram_x = start_frame:end_frame;
    end

    % Crop to 20–40 min (72000–140515)
    idx_window = (ethogram_x >= start_frame & ethogram_x <= end_frame);
    if sum(idx_window) == 0
        fprintf('Figure_%d: No ethogram data in 72000–140515 window.\n', fig_indices(i));
        continue;
    end

    cropped_ethogram = ethogram_y(idx_window);

    % Convert to binary (target ID > 0 = pursuit)
    binary_ethogram = cropped_ethogram > 0;

    % Resize ground truth to match
    this_gt = ground_truth_ethogram(1:length(binary_ethogram));

    % Compare and store
    ethogram_data_all{i} = binary_ethogram;
    match_scores(i) = mean(binary_ethogram(:) == this_gt(:));

    close(fig);
end

%% --- 3. Display Best Matches ---
[sorted_scores, sort_idx] = sort(match_scores, 'descend');

fprintf('\n--- Match Scores vs Ground Truth ---\n');
for i = 1:length(sort_idx)
    fprintf('Figure_%d --> Match: %.2f%%\n', fig_indices(sort_idx(i)), sorted_scores(i)*100);
end

%% --- 4. Plot Best Match ---
best_idx = sort_idx(1);
figure;
plot(ground_truth_ethogram, 'k', 'LineWidth', 1.5); hold on;
plot(ethogram_data_all{best_idx}, 'r');
legend('Ground Truth', sprintf('Figure_%d', fig_indices(best_idx)));

title(sprintf('Best Match: Figure_%d (%.2f%%)', fig_indices(best_idx), sorted_scores(1)*100));
ylim([-0.1, 1.1]);
xlabel('Frame (20–40 min)');

%% --- 5. Plot ALL Ethograms in a Grid (Fixed version) ---
FPS = 60;
time_axis = (frame_window - frame_window(1)) / FPS / 60; % in minutes

figure('Name', 'Ethogram Comparisons: Ground Truth vs Computed', 'Position', [100 100 1200 800]);
for i = 1:length(fig_indices)
    subplot(nrows, ncols, i);
    
    this_comp = ethogram_data_all{i};
    
    if isempty(this_comp)
        title(sprintf('Fig %d: No data', fig_indices(i)));
        continue;
    end

    % Resize ground truth to match length of computed ethogram
    this_gt = ground_truth_ethogram(1:length(this_comp));

    % Ensure computed ethogram is binary
    binary_comp = this_comp > 0;

    % Time axis for plotting
    time_axis_local = time_axis(1:length(binary_comp));

    % Plot
    plot(time_axis_local, this_gt, 'k', 'LineWidth', 1.2); hold on;
    plot(time_axis_local, binary_comp, 'r');
    ylim([-0.1, 1.1]);
    xlim([time_axis(1), time_axis(end)]);
    
    title(sprintf('Fig %d: %.2f%%', fig_indices(i), match_scores(i)*100));
    xlabel('Time (min)');
    yticks([0 1]);
    yticklabels({'Off','On'});
end
legend({'Ground Truth', 'Computed'}, 'Location', 'best');

%% --- 6. Save Summary Table with Match Scores ---
match_table = table(fig_indices(:), match_scores(:)*100, ...
    strcat("Figure_", string(fig_indices(:)), ".fig"), ...
    'VariableNames', {'FigureNumber', 'MatchScorePercent', 'FileName'});

% Display top matches in console
disp('--- Top Matches ---');
disp(sortrows(match_table, 'MatchScorePercent', 'descend'));

% Save the table
save('ethogram_match_scores.mat', 'match_table');
writetable(match_table, 'ethogram_match_scores.csv');
fprintf('Saved match scores as .mat and .csv\n');

FPS = 60;
time_axis = (frame_window - frame_window(1)) / FPS / 60; % in minutes

figure('Name', 'Ethogram Comparisons: Ground Truth vs Computed (Blocks)', 'Position', [100 100 1400 900]);
ncols = 3;
nrows = ceil(length(fig_indices) / ncols);

for i = 1:length(fig_indices)
    subplot(nrows, ncols, i);
    
    this_comp = ethogram_data_all{i};
    if isempty(this_comp)
        title(sprintf('Fig %d: No data', fig_indices(i)));
        continue;
    end
    
    % Resize ground truth to match length of computed ethogram
    this_gt = ground_truth_ethogram(1:length(this_comp));
    binary_comp = this_comp > 0;
    
    % Use time axis cropped for current length
    t = time_axis(1:length(binary_comp));
    
    % Plot ground truth as blocks at y=1
    hold on;
    for k = 1:length(t)-1
        if this_gt(k)
            fill([t(k) t(k+1) t(k+1) t(k)], [0.6 0.6 1.4 1.4], 'k', 'EdgeColor', 'none');
        end
    end
    
    % Plot computed ethogram as blocks at y=0
    for k = 1:length(t)-1
        if binary_comp(k)
            fill([t(k) t(k+1) t(k+1) t(k)], [-0.4 -0.4 0.4 0.4], 'r', 'EdgeColor', 'none');
        end
    end
    
    ylim([-0.5 1.5]);
    xlim([t(1) t(end)]);
    yticks([-0.4, 1]);
    yticklabels({'Computed','Ground Truth'});
    xlabel('Time (min)');
    title(sprintf('Fig %d: Match %.2f%%', fig_indices(i), match_scores(i)*100));
    box on;
end