% Script for Amy to Detect Pursuit Bouts and Durations for Males in Varying Female Density, 3 females per male in each chamber

cd('/')%change accordingly
load('registered_trx.mat')
load('movie-track.mat')

%% Parameters
%remember first is always the male
fly_IDs = [
    1, 2, 3, 4;
    5, 6, 7, 8;
    9, 10, 11, 12;
    13, 14, 15, 16;
    17, 18, 19, 20;
    21, 22, 23, 24;
    25, 26, 27, 28;
    29, 30, 31, 32;
    33, 34, 35, 36
];
FPS = 60;
for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);

%% Plot time series
figure(1); clf
for chamber = 1:size(fly_IDs, 1)
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:end);

    subplot(size(fly_IDs,1), 1, chamber);
    for f = 1:length(female_ids)
        f_id = female_ids(f);
        dist = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
        plot(trx(male_id).timestamps(1:endframe_all), dist, 'LineWidth', 2); hold on
    end
    legend(arrayfun(@(x) sprintf('Target %d', x), 1:length(female_ids), 'UniformOutput', false))
    box off
    ylabel('Male-Target Distances (mm)', 'FontSize', 12)
    xlabel('Time (s)', 'FontSize', 12)
end

%% Raster plot
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


%% Stacked Bar Graphs: % Time spent with each target

figure(5); clf
clear bar_stats target_stats

n_chambers = size(fly_IDs, 1);

for chamber = 1:n_chambers
    male_id = fly_IDs(chamber, 1);
    female_ids = fly_IDs(chamber, 2:4);

    n_inds = zeros(1,3);  % to hold counts per female

    for k = 1:3
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
    none_inds = total_inds - sum(n_inds);

    % Percent time over total
    bar_stats(chamber,:) = ([n_inds(1), none_inds] ./ total_inds) * 100;
   
    % Percent pursuit time only
    if sum(n_inds) > 0
    target_stats(chamber,:) = n_inds ./ sum(n_inds) * 100;
else
    target_stats(chamber,:) = [0, 0, 0];
    end

end

% Sort target preferences for bar plotting
target_stats = sort(target_stats, 2, 'descend');
bar_stats(:,1:3) = sort(bar_stats(:,1:3), 2, 'descend');

% Plot 1: Percent total time
subplot(1,2,1)
bar(bar_stats, 'stacked')
box off
ylabel('% Total Time', 'FontSize', 12)
ylim([0 100])
xlabel('Expreiment No. (60 Minutes)', 'FontSize', 12) %change for the 1hr long
legend('Target 1', 'Target 2', 'Target 3', 'Disengaged', 'FontSize', 12, 'Location', 'northeast')

% Plot 2: Percent pursuit time only
subplot(1,2,2)
bar(target_stats, 'stacked')
box off
ylabel('% Pursuit Time', 'FontSize', 12)
ylim([0 100])
xlabel('Experiment No. (60 Minutes)', 'FontSize', 12) %change for 1hr experiments
legend('Target 1', 'Target 2', 'Target 3', 'FontSize', 12, 'Location', 'northeast')

%% Save Results for Further Analysis

% Define the directory to save the results
results_dir = '/Volumes/otopaliklab/Amy/2025-06-23/MultibbubleVariableDensity_multibubble__WhiteOnly1hour_CSMH_2x1_1hr_CSMH_20250623T080740'; % Change date and name of folder experiment accordingly

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

    % Sanitize figure name: replace spaces and punctuation
    fig_name = regexprep(fig_name, '[^\w]', '_');

    % Bring figure to focus (optional, but avoids issues with invisible figures)
    figure(fig.Number);

    % Save as .png and .fig
    saveas(fig, fullfile(results_dir, [fig_name '.png']));
    savefig(fig, fullfile(results_dir, [fig_name '.fig']));
end

disp(['Saved all figures to: ' results_dir])