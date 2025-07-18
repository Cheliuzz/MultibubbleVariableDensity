%Script for Amy to Detect Pursuit Bouts and Durations for Males in Varying Female Density, to be revised by Adriane

cd('/Volumes/otopaliklab/flydisco_data/2025-06-23/MultibbubleVariableDensity_multibubble__WhiteOnly1hour_CSMH_2x1_1hr_CSMH_20250623T080740') %change accordingly
load('registered_trx.mat')
load('movie-track.mat')


%% Parameters, change as necessary:

fly_IDs = [
    1,  2,  4;
    4,  5,  6;
    7, 8, 9;
    10, 11, 12;
    13, 14, 15;
    16, 17, 18;
    19, 20, 21;
    22, 23, 24;
    25, 26, 27;
  
];
FPS = 60;  
% Compute endframe across all flies
for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);

%% Plot time series of inter-fly Distances across video:
figure(1); clf
for chamber = 1:size(fly_IDs, 1)
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:3);

    for f = 1:2
        f_id = female_ids(f);

        % Use the original structure
        figure(1)
        clear dist
        for i = 1:endframe_all
           dist(i) = pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); ...
                            trx(f_id).x_mm(i),   trx(f_id).y_mm(i)]); % inter-fly distance (mm)
        end

        subplot(size(fly_IDs, 1), 1, chamber);
        plot(trx(male_id).timestamps(1:endframe_all), dist, 'LineWidth', 2); hold on

        % Optional extra figure just for chamber 2
        if chamber == 2
            figure(111)
            plot(trx(male_id).timestamps(1:endframe_all), dist, 'LineWidth', 2); hold on
            ax = gca;
            ax.FontSize = 14;
            xlabel('Time (s)', 'FontSize', 16)
            ylabel('Male-Target Distances (mm)', 'FontSize', 16)
        end
    end

    legend('Target 1', 'Target 2', 'FontSize', 12)
    box off
    ylabel('Male-Target Distances (mm)', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
end

%% Raster plot of pursuit bouts: Inter-fly Distances < threshold distance

% Define thresholds
time_threshold = 10 * FPS;    % 10*FPS = 5 seconds
dist_threshold = 5;           % mm
join_threshold = 3;           % seconds

figure(2); clf

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:3);

    % For each female
    for f = 1:2
        f_id = female_ids(f);

        % Compute distances
        bin = zeros(endframe_all, 1);
        for i = 1:endframe_all
            dist(i) = pdist([
                trx(male_id).x_mm(i), trx(male_id).y_mm(i);
                trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
            ]);
        end

        % Binary vector: 1 if below threshold
        inds = find(dist < dist_threshold);
        bin(inds) = 1;

        % First round of bout detection
        [bouts, lens] = detect_binarybouts(bin);

        % Join bouts separated by < join_threshold seconds
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end

        % Re-detect after joining
        [bouts, lens] = detect_binarybouts(bin);

        % Remove bouts shorter than time threshold
        for i = 1:length(lens)
            if lens(i) < time_threshold
                bin(bouts(i,1):bouts(i,2)) = 0;
            end
        end

        % Final bout detection
        [bouts, lens] = detect_binarybouts(bin);

        % Raster plot for this female
        subplot(size(fly_IDs, 1), 1, chamber)
        inds = find(bin == 1);
        plot(trx(male_id).timestamps(inds), ...
             f * ones(size(inds)), '|');  % f = 1 or 2
        hold on
    end

    box off
    ylim([0.5 2.5])
    ylabel('Target No.', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
end

%% Ethogram: Pursuit bouts per chamber

% Thresholds
time_threshold = 5 * FPS;     % 5 seconds
dist_threshold = 5;           % mm
join_threshold = 3;           % seconds

figure(33); clf

n_chambers = size(fly_IDs, 1);

for chamber = 1:n_chambers
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:3);

    % For each female
    for f = 1:2
        f_id = female_ids(f);

        % Initialize
        bin = zeros(length(trx(male_id).x_mm), 1);
        dist = zeros(endframe_all, 1);

        % Compute inter-fly distance
        for i = 1:endframe_all
            dist(i) = pdist([
                trx(male_id).x_mm(i), trx(male_id).y_mm(i);
                trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
            ]);
        end

        % Binary pursuit detection
        bin(dist < dist_threshold) = 1;

        [bouts, lens] = detect_binarybouts(bin);

        % Join short gaps
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end

        % Remove short bouts
        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)
            if lens(i) < time_threshold
                bin(bouts(i,1):bouts(i,2)) = 0;
            end
        end

        [bouts, lens] = detect_binarybouts(bin);  % Final detection

        % Plot ethogram
        figure(3); hold on
        inds = find(bin == 1);
        plot(trx(male_id).timestamps(inds), chamber * ones(size(inds)), '|', 'MarkerSize', 50);

        % Optional: plot example chamber 2 separately
        if chamber == 2
            figure(111); hold on
            plot(trx(male_id).timestamps(inds), -1 * ones(size(inds)), '|', 'MarkerSize', 10);
        end
    end
    ylim([0 n_chambers + 1])
ylabel('Chamber No.', 'FontSize', 12)
xlabel('Time (s)', 'FontSize', 12)


    figure(3)
    box off
end

figure(3)
ylim([0 n_chambers + 1])
ylabel('Chamber No.', 'FontSize', 12)
xlabel('Time (s)', 'FontSize', 12)

%% Stacked Bar Graphs: % Time spent with each target

figure(5); clf
clear bar_stats target_stats

n_chambers = size(fly_IDs, 1);

for chamber = 1:n_chambers
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:3);

    n_inds = zeros(1,2);  % to hold counts per female

    for k = 1:2
        f_id = female_ids(k);
        dist = zeros(1, endframe_all);

        for i = 1:endframe_all
            dist(i) = pdist([
                trx(male_id).x_mm(i), trx(male_id).y_mm(i);
                trx(f_id).x_mm(i),    trx(f_id).y_mm(i)
            ]);
        end

        inds = find(dist < dist_threshold);
        n_inds(k) = length(inds);  % frames pursuing this target
    end

    total_inds = length(trx(male_id).x_mm);
    none_inds = total_inds - n_inds(1) - n_inds(2);

    % Percent time over total
    bar_stats(chamber,:) = ([n_inds(1), n_inds(2), none_inds] ./ total_inds) * 100;

    % Percent pursuit time only (excluding disengaged)
    if (n_inds(1) + n_inds(2)) > 0
        target_stats(chamber,:) = [n_inds(1), n_inds(2)] ./ sum(n_inds) * 100;
    else
        target_stats(chamber,:) = [0, 0];  % edge case
    end
end

% Sort target preferences for bar plotting
target_stats = sort(target_stats, 2, 'descend');
bar_stats(:,1:2) = sort(bar_stats(:,1:2), 2, 'descend');

% Plot 1: Percent total time
subplot(1,2,1)
bar(bar_stats, 'stacked')
box off
ylabel('% Total Time', 'FontSize', 12)
ylim([0 100])
xlabel('Expreiment No. (60 Minutes)', 'FontSize', 12) %change for the 1hr long
legend('Target 1', 'Target 2', 'Disengaged', 'FontSize', 12, 'Location', 'northeast')

% Plot 2: Percent pursuit time only
subplot(1,2,2)
bar(target_stats, 'stacked')
box off
ylabel('% Pursuit Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12) %change for 1hr experiments
legend('Target 1', 'Target 2', 'FontSize', 12, 'Location', 'northeast')


%% Compute Pursuit Bouts and Durations

dist_threshold = 5;
join_threshold = 3;
time_threshold = 10 * FPS;

pursuit_summary = struct();

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    total_bouts = 0;
    bout_durations = [];
    target_counts = zeros(1, length(female_ids));

    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        bin = dist < dist_threshold;

        [bouts, lens] = detect_binarybouts(bin);
        for i = 1:length(lens)-1
            if bouts(i+1,1) - bouts(i,2) < join_threshold * FPS
                bin(bouts(i,2):bouts(i+1,1)) = 1;
            end
        end

        [bouts, lens] = detect_binarybouts(bin);
               valid_bouts = lens >= time_threshold;
        lens = lens(valid_bouts);

        bout_durations = [bout_durations, (lens / FPS)'];
        total_bouts = total_bouts + length(lens);
        target_counts(f) = target_counts(f) + length(lens);
    end

    pursuit_summary(chamber).chamber = chamber;
    pursuit_summary(chamber).male_id = male_id;
    pursuit_summary(chamber).total_bouts = total_bouts;
    pursuit_summary(chamber).mean_bout_duration = mean(bout_durations);
    pursuit_summary(chamber).std_bout_duration = std(bout_durations);
    pursuit_summary(chamber).target_counts = target_counts;
end

%Convert to Table and Save
summary_table = struct2table(pursuit_summary);
disp(summary_table);
%% Total Number of Pursuit Bouts per Chamber
figure;
total_bouts = arrayfun(@(s) s.total_bouts, pursuit_summary);
bar(total_bouts);
xlabel('Chamber');
ylabel('Total Pursuit Bouts');
title('Total Pursuit Bouts per Chamber');
box off

%% Mean bouts duration per chamber with STD
figure;
mean_durations = arrayfun(@(s) s.mean_bout_duration, pursuit_summary);
std_durations  = arrayfun(@(s) s.std_bout_duration, pursuit_summary);

bar(mean_durations);
hold on
errorbar(1:length(mean_durations), mean_durations, std_durations, '.k', 'LineWidth', 1.2);
xlabel('Chamber');
ylabel('Mean Bout Duration (s)');
title('Mean Pursuit Bout Duration per Chamber');
box off

%% Stacked Bar Plot: Target Preference (Bout Count per Female)

% Stack data (ensure consistent size)
max_targets = max(cellfun(@(s) length(s), {pursuit_summary.target_counts}));
stacked_counts = zeros(length(pursuit_summary), max_targets);

for i = 1:length(pursuit_summary)
    counts = pursuit_summary(i).target_counts;
    stacked_counts(i, 1:length(counts)) = counts;
end

figure;
bar(stacked_counts, 'stacked');
xlabel('Chamber');
ylabel('Pursuit Bouts');
title('Target Preference per Chamber');
legend(arrayfun(@(x) sprintf('Target %d', x), 1:max_targets, 'UniformOutput', false));
box off
%% Save Results for Further Analysis

% Define the directory to save the results
results_dir = '/Volumes/otopaliklab/Variable_Density_data/Variable Density 2 to 1/2025-05-17/fourbubble_multibubble__WhiteOnly1hour_cs2to1_CSMH_20250517T124000'; % Change date and name of folder experiment accordingly

% Create directory if it doesn't exist
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Save all open figures
figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);

    % Use figure number if Name is empty
    fig_name = get(fig, 'Name');
    if isempty(fig_name)
        fig_name = ['Figure_' num2str(fig.Number)];
    end

    % Sanitize figure name
    fig_name = regexprep(fig_name, '[^\w]', '_');

    % Bring figure to focus
    figure(fig.Number);

    % Save as .png and .fig
    saveas(fig, fullfile(results_dir, [fig_name '.png']));
    savefig(fig, fullfile(results_dir, [fig_name '.fig']));
end

% Save summary table
save(fullfile(results_dir, 'pursuit_summary_table.mat'), 'summary_table');
writetable(summary_table, fullfile(results_dir, 'pursuit_summary_table.csv'));

disp(['Saved all figures and summary_table to: ' results_dir])