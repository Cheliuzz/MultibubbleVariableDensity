%Script for Amy to Detect Pursuit Bouts and Durations for Males in Varying Female Density, to be revised by Adriane

cd('/Volumes/otopaliklab/flydisco_data/2025-06-25/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_4x1_1hr_CSMH_20250625T102433')%change accordingly
load('registered_trx.mat')
load('movie-track.mat')
angle_struct = load('perframe/anglefrom1to2_nose2ell.mat');
angle_data = angle_struct.data;
%change and try old videos

%% Parameters
%remember first is always the male
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
time_threshold = 10 * FPS;
dist_threshold = 5;
angle_threshold = 25; % degrees
join_threshold = 3;

% Compute endframe across all flies
for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);
%% Interfly Distance Plot
figure(1); clf
for chamber = 1:size(fly_IDs, 1)
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:end);

    subplot(size(fly_IDs,1), 1, chamber);
    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        plot(trx(male_id).timestamps(1:endframe_all), dist, 'LineWidth', 1); hold on
    end
    %legend(arrayfun(@(x) sprintf('Target %d', x), 1:length(female_ids), 'UniformOutput', false))
    box off
    ylabel('Distances (mm)', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
    title('Male-Female Distances (mm)');
end
%% Raster Plot with Angle Constraint
figure(2); clf
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

figure(5); clf
clear bar_stats target_stats

n_chambers = size(fly_IDs, 1);
n_targets = size(fly_IDs, 2) - 1;
angle_threshold = 20; % degrees
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
bar(bar_stats, 'stacked')
box off
ylabel('% Total Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12)

% Auto legend: Target 1, 2, ..., N, Disengaged
legend_labels = [arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), 'Disengaged'];
legend(legend_labels, 'FontSize', 12, 'Location', 'northeast')

% Plot 2: Percent pursuit time only
subplot(1,2,2)
bar(target_stats, 'stacked')
box off
ylabel('% Pursuit Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12)

%legend(arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), ...
 %   'FontSize', 12, 'Location', 'northeast')

%% Compute Pursuit Bouts and Durations (With Angle Filter)
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
figure(50); clf
clear bar_stats target_stats

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

total_bouts = arrayfun(@(s) s.total_bouts, pursuit_summary);
bar(total_bouts);
xlabel('Chamber');
ylabel('Total Pursuit Bouts');
title('Total Pursuit Bouts per Chamber');
box off

%% Bout Duration Plot per Chamber: Bar + Scatter + SEM
figure(60); clf

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

xlabel('Chamber');
ylabel('Bout Duration (s)');
title('Mean Pursuit Bout Duration per Chamber');
box off

%% Stacked Bar Plot: Target Preference (Bout Count per Female)

%Stack data (ensure consistent size)
%max_targets = max(cellfun(@(s) length(s), {pursuit_summary.target_counts}));
%stacked_counts = zeros(length(pursuit_summary), max_targets);

%for i = 1:length(pursuit_summary)
  %  counts = pursuit_summary(i).target_counts;
 %   stacked_counts(i, 1:length(counts)) = counts;
%end

%figure;
%bar(stacked_counts, 'stacked');
%xlabel('Chamber');
%ylabel('Pursuit Bouts');
%title('Target Preference per Chamber');
%legend(arrayfun(@(x) sprintf('Target %d', x), 1:max_targets, 'UniformOutput', false));
%box off

%% Plot Overlapping Pursuit Periods
figure;
for chamber = 1:size(fly_IDs,1)
    subplot(size(fly_IDs,1),1,chamber)
    male_id = fly_IDs(chamber,1);
    angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));
    female_ids = fly_IDs(chamber,2:end);
    pursuit_matrix = false(length(female_ids), endframe_all);

    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([
            trx(male_id).x_mm(i), trx(male_id).y_mm(i);
            trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
        ]), 1:endframe_all);

        bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

         %Join and filter
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

    overlapping = sum(pursuit_matrix,1) > 1;
    plot(trx(male_id).timestamps(overlapping), chamber*ones(1,sum(overlapping)), '|r', 'MarkerSize', 4); hold on
    title(sprintf('Chamber %d - Overlapping Pursuit Frames', chamber));
    ylabel('Chamber');
    xlabel('Time (s)');
end
%% Determine Primary Pursuit Target (Closest) During Overlaps
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

%% Compare filtered vs. unfiltered bins for a sample male & female
%sample_chamber = 1;
%male_id = fly_IDs(sample_chamber, 1);
%female_id = fly_IDs(sample_chamber, 2);

% Distances and angles
%dist = arrayfun(@(i) pdist([
 %   trx(male_id).x_mm(i), trx(male_id).y_mm(i);
  %  trx(female_id).x_mm(i), trx(female_id).y_mm(i)
%]), 1:endframe_all);
%angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

% Original bin
%raw_bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

% Filtered bin (using same logic as above)
%bin = raw_bin;

%[bouts, lens] = detect_binarybouts(bin);
%for i = 1:length(lens)-1
 %   if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
  %      bin(bouts(i,2):bouts(i+1,1)) = 1;
   % end
%end

%[bouts, lens] = detect_binarybouts(bin);
%for i = 1:length(lens)
 %   if lens(i) < time_threshold
  %      bin(bouts(i,1):bouts(i,2)) = 0;
   % end
%end

% Plot comparison
%figure;
%subplot(2,1,1);
%plot(raw_bin, 'b'); ylim([-0.1, 1.1])
%title('Raw Bouts (Distance + Angle only)')
%xlabel('Frame'); ylabel('Raw State')

%subplot(2,1,2);
%plot(bin, 'r'); ylim([-0.1, 1.1])
%title('Filtered Bouts (After Join & Time Threshold)')
%xlabel('Frame'); ylabel('Filtered State')

%% Save Results for Further Analysis
results_dir = '/Volumes/otopaliklab/Chelo/VariableDensity/Thresholds/07032025/2025-06-26/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_6x1_1hr_CSMH_20250626T092344';
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
writetable(summary_table, fullfile(results_dir, 'pursuit_summary_table.csv'));

disp(['Saved all figures and summary_table to: ' results_dir])


