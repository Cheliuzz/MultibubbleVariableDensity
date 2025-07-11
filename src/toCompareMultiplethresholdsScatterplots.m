% To compare multiple threshold combinations:
cd('/Volumes/otopaliklab/flydisco_data/2025-06-25/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_4x1_1hr_CSMH_20250625T102433')
load('registered_trx.mat')
load('movie-track.mat')
angle_struct = load('perframe/anglefrom1to2_nose2ell.mat');
angle_data = angle_struct.data;

%% Parameters
fly_IDs = [
    1, 2, 3, 4, 5, 6, 7;
    %8, 9,10,11,12,13,14;
    16,15,17,18,19,20,21;
    23,25,24,22,26,27,28;
    30,29,31,32,33,34,35;
    37,36,38,39,40,41,42;
    %44,43,45,44,47,48,49;
    52,50,51,53,54,55,56;
    59,58,63,60,61,62,57;
];


FPS = 60;
chamber_to_plot = 1;
angle_threshold = 25; % degrees
female_ids = fly_IDs(chamber_to_plot, 2:end);

for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);

% Threshold ranges
time_thresholds = [5, 8, 10] * FPS;
dist_thresholds = [3, 5, 8, 10];



%% --- Updated Pursuit Analysis Script with Per-Chamber Duration Tracking and Visualizations ---

% Define and initialize result storage
summary_all = struct();
combo = 1;

for tt = 1:length(time_thresholds)
    for dt = 1:length(dist_thresholds)

        time_thr = time_thresholds(tt);
        dist_thr = dist_thresholds(dt);

        target_counts_all = zeros(size(fly_IDs,1), length(female_ids));
        total_bouts_all = zeros(size(fly_IDs,1), 1);
        bout_durations_all = cell(size(fly_IDs, 1), 1);
        mean_bout_durations = nan(size(fly_IDs, 1), 1);

        for chamber = 1:size(fly_IDs,1)
            male_id = fly_IDs(chamber, 1);
            f_ids = fly_IDs(chamber, 2:end);
            angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));
            target_per_frame = nan(1, endframe_all);
            bout_durations = [];

            for f = 1:length(f_ids)
                f_id = f_ids(f);
                dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
                bin = (dist < dist_thr) & (angle_deg < angle_threshold);

                [bouts, lens] = detect_binarybouts(bin);
                for i = 1:length(lens)-1
                    if bouts(i+1,1) - bouts(i,2) < 3
                        bin(bouts(i,2):bouts(i+1,1)) = 1;
                    end
                end

                [bouts, lens] = detect_binarybouts(bin);
                for i = 1:length(lens)
                    if lens(i) < time_thr
                        bin(bouts(i,1):bouts(i,2)) = 0;
                    else
                        target_per_frame(bouts(i,1):bouts(i,2)) = f;
                        bout_durations(end+1) = lens(i) / FPS;
                    end
                end
            end

            for f = 1:length(f_ids)
                target_counts_all(chamber, f) = sum(target_per_frame == f);
            end
            total_bouts_all(chamber) = numel(bout_durations);
            bout_durations_all{chamber} = bout_durations;
            if ~isempty(bout_durations)
                mean_bout_durations(chamber) = mean(bout_durations);
            end
        end

        summary_all(combo).time_thr = time_thr / FPS;
        summary_all(combo).dist_thr = dist_thr;
        summary_all(combo).target_counts = target_counts_all;
        summary_all(combo).total_bouts = total_bouts_all;
        summary_all(combo).mean_bout_duration = mean(cellfun(@mean, bout_durations_all, 'UniformOutput', true), 'omitnan');
        summary_all(combo).mean_bout_duration_per_chamber = mean_bout_durations;
        summary_all(combo).bout_durations_all = bout_durations_all;

        combo = combo + 1;
    end
end

%% --- Create Per-Chamber Bout Duration Table ---
num_combos = length(summary_all);
rows = [];
chamber_ids = (1:size(fly_IDs, 1))';

for i = 1:num_combos
    time_thr = summary_all(i).time_thr;
    dist_thr = summary_all(i).dist_thr;
    durations = summary_all(i).mean_bout_duration_per_chamber;
    for ch = 1:length(durations)
        rows = [rows;
            time_thr, dist_thr, chamber_ids(ch), durations(ch), summary_all(i).total_bouts(ch)];
    end
end

bout_durations_chamber_table = array2table(rows, ...
    'VariableNames', {'Time_Threshold_s', 'Dist_Threshold_mm', 'Chamber', 'Mean_Bout_Duration_s', 'Total_Bouts'});

% Save table
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
writetable(bout_durations_chamber_table, fullfile(results_dir, 'bout_durations_per_chamber.csv'));
save(fullfile(results_dir, 'bout_durations_per_chamber.mat'), 'bout_durations_chamber_table');

%% --- Visualization: Bout Count vs. Mean Duration ---

% 1. Simple 2D scatter (color by time threshold, size by distance threshold)
figure;
scatter(bout_durations_chamber_table.Total_Bouts, ...
        bout_durations_chamber_table.Mean_Bout_Duration_s, ...
        40, bout_durations_chamber_table.Time_Threshold_s + 0.01 * bout_durations_chamber_table.Dist_Threshold_mm, ...
        'filled');
xlabel('Total Bouts per Chamber');
ylabel('Mean Bout Duration (s)');
title('Bout Count vs Duration (Numeric Encoding of Thresholds)');
colorbar;
grid on;

% 2. 3D scatter: Count, Duration, and Distance (colored by time threshold)
figure;
scatter3(bout_durations_chamber_table.Total_Bouts, ...
         bout_durations_chamber_table.Mean_Bout_Duration_s, ...
         bout_durations_chamber_table.Dist_Threshold_mm, ...
         40, bout_durations_chamber_table.Time_Threshold_s, ...
         'filled');
xlabel('Total Bouts');
ylabel('Mean Duration (s)');
zlabel('Distance Threshold (mm)');
title('3D Scatter: Count, Duration, and Distance');
colorbar;
grid on;

% 3. Grouped 2D scatter using categorical labels
combo_labels = strcat("T", string(bout_durations_chamber_table.Time_Threshold_s), ...
                      "-D", string(bout_durations_chamber_table.Dist_Threshold_mm));
combo_categorical = categorical(combo_labels, unique(combo_labels, 'stable'));

figure;
gscatter(bout_durations_chamber_table.Total_Bouts, ...
         bout_durations_chamber_table.Mean_Bout_Duration_s, ...
         combo_categorical);
xlabel('Total Bouts per Chamber');
ylabel('Mean Bout Duration (s)');
title('Bout Count vs Duration by Threshold Combination');
legend('Location', 'bestoutside');
grid on;
%% box plot

% Prepare boxplot data
num_combos = length(summary_all);
combo_labels = strings(1, num_combos);
box_data = [];

for i = 1:num_combos
    box_data = [box_data; summary_all(i).total_bouts(:)];
    combo_labels(i) = sprintf('T%ds-D%dmm', ...
        summary_all(i).time_thr, summary_all(i).dist_thr);
end

% Generate group indices for boxplot
group = repelem(1:num_combos, size(summary_all(1).total_bouts,1))';

% Plot
figure;
boxplot(box_data, group, 'Labels', combo_labels);
xtickangle(45);
ylabel('Pursuit Bouts per Chamber');
xlabel('Threshold Combination');
title('Variation Across Chambers for Each Threshold Combination');


%% figure;
hold on
boxplot(box_data, group, 'Labels', combo_labels);
xtickangle(45);
ylabel('Pursuit Bouts per Chamber');
xlabel('Threshold Combination');
title('Boxplot + Chamber-Level Points');

% Overlay scatter dots (jittered x positions)
for i = 1:num_combos
    y = summary_all(i).total_bouts(:);
    x = i + 0.15*randn(size(y));  % jitter around the box
    scatter(x, y, 25, 'filled', 'MarkerFaceAlpha', 0.6);
end
hold off

%% Save
results_dir = '/Volumes/otopaliklab/Chelo/VariableDensity/Thresholds/07102025/MultipleParams/2025-06-26/20mins/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_6x1_1hr_CSMH_20250626T081139';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);
    fig_name = get(fig, 'Name');
    if isempty(fig_name), fig_name = ['Figure_' num2str(fig.Number)]; end
    fig_name = regexprep(fig_name, '[^\w]', '_');
    figure(fig.Number);
    saveas(fig, fullfile(results_dir, [fig_name '.png']));
    savefig(fig, fullfile(results_dir, [fig_name '.fig']));
end
save(fullfile(results_dir, 'bout_summary_table.mat'), 'bout_summary');
save(fullfile(results_dir, 'summary_all_struct.mat'), 'summary_all');
summary_table = struct2table(summary_all);
writetable(summary_table, fullfile(results_dir, 'summary_all_table.csv'));
save(fullfile(results_dir, 'summary_all_table.mat'), 'summary_table');

disp(['Saved all figures and summary_table to: ' results_dir])


