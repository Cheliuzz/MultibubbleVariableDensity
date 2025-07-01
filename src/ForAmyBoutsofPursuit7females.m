% Script to Detect Pursuit Bouts and Durations for Males in Varying Female Density
% Experiment: 7 females per male per chamber

cd('/Volumes/otopaliklab/flydisco_data/2025-06-26/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_6x1_1hr_CSMH_20250626T081139'); %modify
load('registered_trx.mat')
load('movie-track.mat')

%% Parameters
fly_IDs = [
    1, 2, 3, 4, 5, 6, 7, 8;
    9,10,11,12,13,14,15,16;
    17,18,19,20,21,22,23,24;
    25,26,27,28,29,30,31,32;
    33,34,35,36,37,38,39,40;
    41,42,43,44,45,46,47,48;
    49,50,51,52,53,54,55,56;
    57,58,59,60,61,62,63,64;
    65,66,67,68,69,70,71,72
    
    % modify if needed
];

FPS = 60;
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
    legend(arrayfun(@(x) sprintf('Target %d', x), 1:length(female_ids), 'UniformOutput', false))
    box off
    ylabel('Male-Target Distances (mm)', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
end

%% Raster Plot
figure(2); clf
time_threshold = 10 * FPS;
dist_threshold = 5;
join_threshold = 3;

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);
    subplot(size(fly_IDs,1), 1, chamber)

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
time_threshold = 5 * FPS;

for chamber = 1:size(fly_IDs,1)
    male_id = fly_IDs(chamber,1);
    female_ids = fly_IDs(chamber,2:end);

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

%% Stacked Bar Graphs
figure(5); clf
clear bar_stats target_stats

n_chambers = size(fly_IDs, 1);
n_targets = size(fly_IDs, 2) - 1;

for chamber = 1:n_chambers
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:end);

    n_inds = zeros(1, n_targets);

    for k = 1:n_targets
        f_id = female_ids(k);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        inds = find(dist < dist_threshold);
        n_inds(k) = length(inds);
    end

    total_inds = length(trx(male_id).x_mm);
    none_inds = total_inds - sum(n_inds);

    bar_stats(chamber,:) = ([n_inds, none_inds] ./ total_inds) * 100;

    if sum(n_inds) > 0
        target_stats(chamber,:) = n_inds ./ sum(n_inds) * 100;
    else
        target_stats(chamber,:) = zeros(1, n_targets);
    end
end

% Sort for plotting
target_stats = sort(target_stats, 2, 'descend');
bar_stats(:,1:n_targets) = sort(bar_stats(:,1:n_targets), 2, 'descend');

subplot(1,2,1)
bar(bar_stats, 'stacked')
box off
ylabel('% Total Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No.', 'FontSize', 12)
legend([arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), 'Disengaged'], 'FontSize', 12)

subplot(1,2,2)
bar(target_stats, 'stacked')
box off
ylabel('% Pursuit Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No.', 'FontSize', 12)
legend(arrayfun(@(x) sprintf('Target %d', x), 1:n_targets, 'UniformOutput', false), 'FontSize', 12)

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
results_dir = '/Volumes/otopaliklab/Amy/'; % Change date and name of folder experiment accordingly

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


