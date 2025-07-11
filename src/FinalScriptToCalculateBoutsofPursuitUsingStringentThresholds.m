% Script to detect bouts of pursuit using facing angle, distance and time
% as thresholds. To be revised by Adriane
cd('/Volumes/otopaliklab/flydisco_data/2025-06-25/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_4x1_1hr_CSMH_20250625T102433')
load('registered_trx.mat')
load('movie-track.mat')
angle_struct = load('perframe/anglefrom1to2_nose2ell.mat');
angle_data = angle_struct.data;
dnose_struct = load('perframe/dnose2ell_angle_min20to20.mat');
dnose_data = dnose_struct.data;

%% Parameters
%you have to modify this for each experiment folder
fly_IDs = [
    1, 2, 3, 4, 5;
    6, 8, 7, 9, 10;
    11, 13, 14, 15, 12;
    16, 20, 17, 18, 19;
    21, 23, 22, 24, 25;
    27, 28, 26, 29, 30;
    31, 33, 32, 34, 35;
    36, 37, 38, 39, 40;
    42, 41, 44, 43, 45;
]; %modify accordingly



FPS = 60;
time_threshold = 10 * FPS; %10 sec calculated from frames
dist_threshold = 5; % mm
angle_threshold = 25; %deg 
join_threshold = 3; % frames
dnose_threshold = 5; % mm

% Compute endframe across all flies
for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);

%% Interfly Distance Plot (Raw Distances)
figure(1); clf
set(gcf, 'Name', 'Figure_DistancePlots');

for chamber = 1:size(fly_IDs, 1)
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:end);
   % n_targets = size(fly_IDs, 2) - 1;% check if it works without this

    subplot(size(fly_IDs,1), 1, chamber);
    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        plot(trx(male_id).timestamps(1:endframe_all), dist, 'LineWidth', 1); hold on
    end
    legend(arrayfun(@(x) sprintf('Target %d', x), 1:length(female_ids), 'UniformOutput', false))
    box off
    ylabel('Distances (mm)', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
end
%% Raster Plot with Angle Constraint
figure(2); clf
set(gcf, 'Name', 'Figure_RasterPlot');
for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    subplot(size(fly_IDs,1), 1, chamber)

    angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end

        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)
            if lens(i) < time_threshold
                bin(bouts(i,1):bouts(i,2)) = 0;
            end
        end

        inds = find(bin == 1);
        plot(trx(male_id).timestamps(inds), f * ones(size(inds)), '|'); hold on
    end
    ylim([0.5 length(female_ids)+0.5])
    ylabel('Target No.', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
  
end
%% Ethogram
figure(3); clf
set(gcf, 'Name', 'Figure_Ethogram');
for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end

        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)
            if lens(i) < time_threshold
                bin(bouts(i,1):bouts(i,2)) = 0;
            end
        end

        inds = find(bin == 1);
        plot(trx(male_id).timestamps(inds), chamber * ones(size(inds)), '|', 'MarkerSize', 50); hold on
    end
end
ylim([0 size(fly_IDs,1) + 1])
ylabel('Chamber No.', 'FontSize', 12)
xlabel('Time (s)', 'FontSize', 12)


%% Stacked Bar Graphs: Time spent with each target (distance + angle threshold)

figure(4); clf
clear bar_stats target_stats
set(gcf, 'Name', 'Figure_StackedBars');

n_chambers = size(fly_IDs, 1);
n_targets = size(fly_IDs, 2) - 1;
angle_threshold = 25; % degrees
dist_threshold = 5;   % mm

% Add at the top of the bar graph section:
fly_id_list = [trx.id];

for chamber = 1:n_chambers
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:end);
    n_inds = zeros(1, n_targets);  % one per female

    for k = 1:n_targets
        f_id = female_ids(k); 
        dist = zeros(1, endframe_all);

        % Compute interfly distances
        for i = 1:endframe_all
            dist(i) = pdist([
                trx(male_id).x_mm(i), trx(male_id).y_mm(i);
                trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
            ]);
        end

        % Convert female ID to row index in angle matrix
        f_index = find(fly_id_list == f_id);

        % Get angle from male to this female
        angle = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

        % Apply combined distance and angle condition
        inds = find((dist < dist_threshold) & (angle < angle_threshold));
        n_inds(k) = length(inds);
    end

    total_inds = length(trx(male_id).x_mm);
    none_inds = total_inds - sum(n_inds);

    % Percent time over total
    bar_stats(chamber,:) = ([n_inds, none_inds] ./ total_inds) * 100;

    % Percent pursuit time only
    if sum(n_inds) > 0
        target_stats(chamber,:) = n_inds ./ sum(n_inds) * 100;
    else
        target_stats(chamber,:) = zeros(1, n_targets);
    end
end
% Sort target preferences for bar plotting
target_stats = sort(target_stats, 2, 'descend');
bar_stats(:,1:n_targets) = sort(bar_stats(:,1:n_targets), 2, 'descend');

% Plot 1: Percent total time
subplot(1,2,1)
b1 = bar(bar_stats, 'stacked');
box off
ylabel('% Total Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12)

b1(end).FaceColor = [0.94 0.94 0.94];% Force the last bar segment (Disengaged) to be light gray

%legend_labels = [arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), 'Disengaged'];
%legend(b1, legend_labels, 'FontSize', 12, 'Location', 'northeast')

% Plot 2: Percent pursuit time only
subplot(1,2,2)
bar(target_stats, 'stacked')
box off
ylabel('% Pursuit Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12)

%legend(arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), ...
 %   'FontSize', 12, 'Location', 'northeast')


 %% Compute Pursuit Bouts and Durations (With Summary)
pursuit_summary = struct();
fly_id_list = [trx.id];

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

    n_targets = length(female_ids);
    dist_all = zeros(n_targets, endframe_all);
    bin_all = false(n_targets, endframe_all);

    % Step 1: Compute bin for each female
    for f = 1:n_targets
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([
            trx(male_id).x_mm(i), trx(male_id).y_mm(i);
            trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
        ]), 1:endframe_all);

        dist_all(f,:) = dist;
        bin_all(f,:) = (dist < dist_threshold) & (angle_deg < angle_threshold);
    end

    % Step 2: For each frame, select the closest female among candidates
    target_per_frame = nan(1, endframe_all); % NaN = no pursuit
    for t = 1:endframe_all
        candidates = find(bin_all(:,t));
        if ~isempty(candidates)
            [~, min_idx] = min(dist_all(candidates, t));
            target_per_frame(t) = candidates(min_idx); % 1, 2, ..., n_targets
        end
    end

    % Step 3: Convert labeled frame array into per-target bouts
    total_bouts = 0;
    bout_durations = [];
    target_counts = zeros(1, n_targets);

    for f = 1:n_targets
        bin = (target_per_frame == f); % binary array: 1 = pursuing female f

        % Join close bouts
        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end

        % Filter short bouts
        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)
            if lens(i) < time_threshold
                bin(bouts(i,1):bouts(i,2)) = 0;
            end
        end

        % Recompute final bouts and stats
        [bouts, lens] = detect_binarybouts(bin);
        bout_durations = [bout_durations, (lens / FPS)'];
        total_bouts = total_bouts + length(lens);
        target_counts(f) = length(lens);
    end

    % Store results
    pursuit_summary(chamber).chamber = chamber;
    pursuit_summary(chamber).male_id = male_id;
    pursuit_summary(chamber).total_bouts = total_bouts;
    pursuit_summary(chamber).mean_bout_duration = mean(bout_durations);
    pursuit_summary(chamber).std_bout_duration = std(bout_durations);
    pursuit_summary(chamber).bout_durations = bout_durations;  
    pursuit_summary(chamber).target_counts = target_counts;
    pursuit_summary(chamber).target_per_frame = target_per_frame;
end

% Convert to table
summary_table = struct2table(pursuit_summary);
disp(summary_table);

 %% Stacked Bar Graphs: Time spent pursuing each female (closest-only logic)
figure(5); clf
clear bar_stats target_stats
set(gcf, 'Name', 'Figure_StackedBar_ClosestOnly');

n_chambers = size(fly_IDs, 1);
n_targets = size(fly_IDs, 2) - 1;

for chamber = 1:n_chambers
    summary = pursuit_summary(chamber);
    target_labels = summary.target_per_frame; % 1 to N or NaN
    total_frames = length(target_labels);

    n_inds = zeros(1, n_targets);
    for k = 1:n_targets
        n_inds(k) = sum(target_labels == k);
    end

    none_inds = sum(isnan(target_labels));

    % Percent of total time
    bar_stats(chamber,:) = ([n_inds, none_inds] ./ total_frames) * 100;

    % Percent of pursuit time only (exclude disengaged)
    if sum(n_inds) > 0
        target_stats(chamber,:) = n_inds ./ sum(n_inds) * 100;
    else
        target_stats(chamber,:) = zeros(1, n_targets);
    end
end

% Sort target preferences for prettier bars
target_stats = sort(target_stats, 2, 'descend');
bar_stats(:,1:n_targets) = sort(bar_stats(:,1:n_targets), 2, 'descend');

% Plot 1: % total time
subplot(1,2,1)
bar(bar_stats, 'stacked')
box off
ylabel('% Total Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12)
title('Time with Each Target (Closest Only)', 'FontSize', 13)

% Plot 2: % of pursuit time
subplot(1,2,2)
bar(target_stats, 'stacked')
box off
ylabel('% Pursuit Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12)
title('Pursuit Target Preference (Closest Only)', 'FontSize', 13)

legend_labels = [arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), 'Disengaged'];
legend(legend_labels, 'FontSize', 10, 'Location', 'northeastoutside');



%% Total Number of Pursuit Bouts per Chamber
figure(6); clf
set(gcf, 'Name', 'Figure_TotalBouts');
total_bouts = arrayfun(@(s) s.total_bouts, pursuit_summary);
bar(total_bouts);
xlabel('Chamber');
ylabel('Total Pursuit Bouts');
title('Total Pursuit Bouts per Chamber');
box off

%% Bout Duration Plot per Chamber: Bar + Scatter + SEM
figure(7); clf
set(gcf, 'Name', 'Figure_BoutDurations');

n_chambers = length(pursuit_summary);
mean_durations = zeros(1, n_chambers);
sem_durations  = zeros(1, n_chambers); % Standard Error of the Mean

all_durations = cell(1, n_chambers); % Store each chamber's durations

for chamber = 1:n_chambers
    bout_durations = pursuit_summary(chamber).bout_durations;
    all_durations{chamber} = bout_durations;
    mean_durations(chamber) = mean(bout_durations);
    sem_durations(chamber) = std(bout_durations) / sqrt(length(bout_durations));
end

% Plot bars for mean
bar(mean_durations, 'FaceColor', [0.6 0.6 0.9]); hold on

% Plot scatter of individual bouts
for chamber = 1:n_chambers
    xvals = repmat(chamber, 1, length(all_durations{chamber})) + 0.1*randn(1,length(all_durations{chamber}));
    scatter(xvals, all_durations{chamber}, 8, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
end

% Error bars = SEM
errorbar(1:n_chambers, mean_durations, sem_durations, '.k', 'LineWidth', 1.2);
% Filter out empty cells
nonempty_durations = all_durations(~cellfun(@isempty, all_durations));

% Compute max safely
if ~isempty(nonempty_durations)
    max_val = max(cellfun(@max, nonempty_durations));
    ylim([0, max_val * 1.1]);
else
    ylim([0 1]); % fallback in case all cells are empty
end
xlabel('Chamber');
ylabel('Bout Duration (s)');
title('Mean Pursuit Bout Duration per Chamber');
box off

%% Determine Primary Pursuit Target (Closest) During Overlaps
figure(8); clf
set(gcf, 'Name', 'Figure_PrimaryTarget_Overlaps');
closest_target_frames = cell(size(fly_IDs,1), 1); % one cell per chamber

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));
    
    pursuit_matrix = false(length(female_ids), endframe_all);
    distances_all = zeros(length(female_ids), endframe_all); % store all distances

    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([
            trx(male_id).x_mm(i), trx(male_id).y_mm(i);
            trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
        ]), 1:endframe_all);
        distances_all(f,:) = dist;

        bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

        % Join and filter
        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end
        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)
            if lens(i) < time_threshold
                bin(bouts(i,1):bouts(i,2)) = 0;
            end
        end

        pursuit_matrix(f,:) = bin;
    end

    % Now identify primary female (closest one in simultaneous pursuit)
    primary_target = nan(1, endframe_all);  % NaN = no pursuit; otherwise, index of closest female
    for t = 1:endframe_all
        active_females = find(pursuit_matrix(:,t)); % indices of females being pursued at t
        if ~isempty(active_females)
            [~, min_idx] = min(distances_all(active_females, t)); % closest among active
            primary_target(t) = active_females(min_idx);  % store index of female (1=first, 2=second, etc.)
        end
    end

    closest_target_frames{chamber} = primary_target;

    % Plot
    subplot(size(fly_IDs,1),1,chamber)
    hold on
    overlapping = sum(pursuit_matrix,1) > 1;
    plot(trx(male_id).timestamps(overlapping), chamber*ones(1,sum(overlapping)), '|r', 'MarkerSize', 4); 
    % Identify valid (non-NaN) frames where male is pursuing someone
valid_idx = ~isnan(primary_target);

% Get timestamps for those frames
timestamps_valid = trx(male_id).timestamps(valid_idx);

% Compute y-positions and color values for scatter plot
y_positions = chamber + 0.2 * (primary_target(valid_idx) - 1) / length(female_ids);
color_values = primary_target(valid_idx);

% Plot
scatter(timestamps_valid, y_positions, 2, color_values, 'filled');

    title(sprintf('Chamber %d - Closest Target (color) & Overlaps (red)', chamber));
    ylabel('Chamber');
    xlabel('Time (s)');
end
%% Simplified Ethogram (Single Chamber, 20–40 min) - Chamber-Level Style

male_chamber = 1; % <-- Change this as needed
male_id = fly_IDs(male_chamber, 1);
female_ids = fly_IDs(male_chamber, 2:end);

% Time window (seconds)
start_time_sec = 20 * 60;
end_time_sec   = 40 * 60;

% Convert to frame indices
start_frame = find(trx(male_id).timestamps >= start_time_sec, 1, 'first');
end_frame   = find(trx(male_id).timestamps <= end_time_sec, 1, 'last');

% Precompute angle data
angle_deg = abs(rad2deg(angle_data{male_id}(start_frame:end_frame)));

% Set up figure
figure(9); clf
set(gcf, 'Name', 'Figure_20mins_Ethogram');
hold on

for f = 1:length(female_ids)
    f_id = female_ids(f);

    % Compute distance
    dist = arrayfun(@(i) pdist([
        trx(male_id).x_mm(i), trx(male_id).y_mm(i);
        trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
    ]), start_frame:end_frame);

    % Apply distance + angle threshold
    bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

    % Join nearby bouts
    [bouts, lens] = detect_binarybouts(bin);
    for i = 1:length(lens)-1
        if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
            bin(bouts(i,2):bouts(i+1,1)) = 1;
        end
    end

    % Filter out short bouts
    [bouts, lens] = detect_binarybouts(bin);
    for i = 1:length(lens)
        if lens(i) < time_threshold
            bin(bouts(i,1):bouts(i,2)) = 0;
        end
    end

    % Plot
    inds = find(bin == 1);
    timestamps = trx(male_id).timestamps(start_frame:end_frame);
    plot(timestamps(inds), f * ones(size(inds)), '|', 'MarkerSize', 50);
end

% Format
ylim([0.5 length(female_ids) + 0.5]);
yticks(1:length(female_ids));
yticklabels(arrayfun(@(i) sprintf('F%d', i), 1:length(female_ids), 'UniformOutput', false));
xlabel('Time (s)');
ylabel('Target Female');
title(sprintf('Simplified Ethogram: Male %d (Chamber %d)', male_id, male_chamber));
box off


%% Compare Pursuit Bouts Across 3 × 20-Minute Bins
figure(10); clf
set(gcf, 'Name', 'ComparePursuitBouts');
n_bins = 3;
bin_duration = 20 * 60; % seconds
bin_edges = 0 : bin_duration : 60*60; % [0, 1200, 2400, 3600]

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    timestamps = trx(male_id).timestamps;
    
    % Calculate number of valid frames shared across all flies
    n_frames = min([
        length(trx(male_id).x_mm), ...
        min(arrayfun(@(f_id) length(trx(f_id).x_mm), female_ids)), ...
        length(timestamps)
    ]);

    angle_deg = abs(rad2deg(angle_data{male_id}));
    angle_deg = angle_deg(1:n_frames);
    n_targets = length(female_ids);
    
    dist_all = zeros(n_targets, n_frames);
    bin_all = false(n_targets, n_frames);

    % Compute distance and logical mask for each female
    for f = 1:n_targets
        f_id = female_ids(f);

        dist_all(f,1:n_frames) = arrayfun(@(i) pdist([
            trx(male_id).x_mm(i), trx(male_id).y_mm(i);
            trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
        ]), 1:n_frames);

        bin_all(f,:) = (dist_all(f,:) < dist_threshold) & (angle_deg < angle_threshold);
    end
    
    % Determine closest female per frame that meets criteria
    target_per_frame = nan(1, n_frames);
    for t = 1:n_frames
        candidates = find(bin_all(:,t));
        if ~isempty(candidates)
            [~, idx] = min(dist_all(candidates, t));
            target_per_frame(t) = candidates(idx);
        end
    end

    % Plot results for this chamber
    subplot(size(fly_IDs,1),1,chamber); hold on
    colors = lines(n_bins);

    for bin = 1:n_bins
        bin_start = bin_edges(bin);
        bin_end   = bin_edges(bin+1);
        
        frame_idx = find(timestamps(1:n_frames) >= bin_start & timestamps(1:n_frames) < bin_end);
        tpf_bin = target_per_frame(frame_idx);

        total_bouts_bin = 0;
        total_duration_bin = 0;

        for f = 1:n_targets
            bin_logic = (tpf_bin == f);
            [bouts, lens] = detect_binarybouts(bin_logic);

            % Join close bouts
            for i = 1:length(lens)-1
                if bouts(i+1,1) - bouts(i,2) < join_threshold
                    bin_logic(bouts(i,2):bouts(i+1,1)) = 1;
                end
            end

            [bouts, lens] = detect_binarybouts(bin_logic);
            lens(lens < time_threshold) = [];

            total_bouts_bin = total_bouts_bin + length(lens);
            total_duration_bin = total_duration_bin + sum(lens); % In frames
        end

        % Plot: bout count (tall bar) and duration (shorter dashed)
        yyaxis left
        bar(bin, total_bouts_bin, 'FaceColor', colors(bin,:), 'EdgeColor', 'none');

        yyaxis right
        bar(bin, total_duration_bin / FPS, 0.4, 'FaceColor', [0.5 0.5 0.5], ...
            'EdgeColor', 'k', 'LineStyle', '--');
    end

    % Labels and formatting
    title(sprintf('Chamber %d - Bouts per 20-min Bin', chamber));

    yyaxis left
    ylabel('Bout Count');

    yyaxis right
    ylabel('Duration (sec)');

    xticks(1:n_bins);
    xticklabels({'0-20', '20-40', '40-60'});
    xlabel('Time Bin (min)');
    ylim([0 inf]);
    box off
end
%% Save Results for Further Analysis
results_dir = '/Volumes/otopaliklab/Chelo/VariableDensity/Thresholds/07102025/MultipleParams/2025-06-26/20mins/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_6x1_1hr_CSMH_20250626T081139';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);
    fig_name = get(fig, 'Name');
    if isempty(fig_name)
        fig_name = ['Figure_' num2str(fig.Number)];
    end
    fig_name = regexprep(fig_name, '[^\w]', '_');
    figure(fig.Number);
    saveas(fig, fullfile(results_dir, [fig_name '.png']));
    savefig(fig, fullfile(results_dir, [fig_name '.fig']));
end

save(fullfile(results_dir, 'pursuit_summary_table.mat'), 'summary_table');
save(fullfile(results_dir, 'pursuit_summary_struct.mat'), 'pursuit_summary');
writetable(summary_table, fullfile(results_dir, 'pursuit_summary_table.csv'));

disp(['Saved all figures and summary_table to: ' results_dir])
