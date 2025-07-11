% Script: Compare Pursuit Thresholds (No Cropping, 20-40 min Ethograms)
cd('/Volumes/otopaliklab/flydisco_data/2025-06-25/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_4x1_1hr_CSMH_20250625T102433')
load('registered_trx.mat')
load('movie-track.mat')
angle_struct = load('perframe/anglefrom1to2_nose2ell.mat');
angle_data = angle_struct.data;

%% Parameters
FPS = 60;
analysis_start_sec = 20 * 60;
analysis_end_sec   = 40 * 60;
angle_threshold = 25; % degrees
join_threshold = 3; % frames

fly_IDs = [
    3, 2, 1, 4, 5;
    6, 8, 7, 9, 10;
    11, 13, 14, 15, 12;
    16, 20, 17, 16, 19;
    21, 23, 22, 24, 25;
    27, 28, 26, 29, 30;
    31, 33, 32, 34, 35;
    36, 37, 38, 39, 40;
    42, 41, 44, 43, 45
];

% Thresholds to compare
time_thresholds = [5, 8, 10] * FPS;
dist_thresholds = [3, 5, 8, 10];
chamber_to_plot = 1;

summary_all = struct();
combo = 1;

%% Main Loop
for tt = 1:length(time_thresholds)
    for dt = 1:length(dist_thresholds)

        time_thr = time_thresholds(tt);
        dist_thr = dist_thresholds(dt);

        female_ids = fly_IDs(chamber_to_plot, 2:end);
        target_counts_all = zeros(size(fly_IDs,1), length(female_ids));
        total_bouts_all = zeros(size(fly_IDs,1), 1);
        bout_durations_all = cell(size(fly_IDs, 1), 1);
        mean_bout_durations = nan(size(fly_IDs, 1), 1);

        for chamber = 1:size(fly_IDs,1)
            male_id = fly_IDs(chamber, 1);
            f_ids = fly_IDs(chamber, 2:end);

            timestamps = trx(male_id).timestamps;
            start_idx = find(timestamps >= analysis_start_sec, 1, 'first');
            end_idx   = find(timestamps <= analysis_end_sec, 1, 'last');

            if isempty(start_idx) || isempty(end_idx) || end_idx <= start_idx
                continue;
            end

            time_window = start_idx:end_idx;
            angle_deg = abs(rad2deg(angle_data{male_id}(time_window)));
            nFrames = length(time_window);

            target_per_frame = nan(1, nFrames);
            bout_durations = [];

            for f = 1:length(f_ids)
                f_id = f_ids(f);
                dist = arrayfun(@(i) pdist([trx(male_id).x_mm(time_window(i)), trx(male_id).y_mm(time_window(i));
                                            trx(f_id).x_mm(time_window(i)), trx(f_id).y_mm(time_window(i))]), 1:nFrames);
                bin = (dist < dist_thr) & (angle_deg < angle_threshold);

                if any(bin)
                    [bouts, lens] = detect_binarybouts(bin);
                    for i = 1:length(lens)-1
                        if bouts(i+1,1) - bouts(i,2) < join_threshold
                            bin(bouts(i,2):bouts(i+1,1)) = 1;
                        end
                    end
                    [bouts, lens] = detect_binarybouts(bin);
                else
                    bouts = [];
                    lens = [];
                end

                for i = 1:length(lens)
                    if lens(i) >= time_thr
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

        %% Plot Ethogram for Chamber 1
        if chamber_to_plot == 1
            male_id = fly_IDs(chamber_to_plot, 1);
            female_ids = fly_IDs(chamber_to_plot, 2:end);
            timestamps = trx(male_id).timestamps;
            start_idx = find(timestamps >= analysis_start_sec, 1, 'first');
            end_idx   = find(timestamps <= analysis_end_sec, 1, 'last');
            time_window = start_idx:end_idx;
            angle_deg = abs(rad2deg(angle_data{male_id}(time_window)));
            
            figure(400 + combo); clf; hold on
            for f = 1:length(female_ids)
                f_id = female_ids(f);
                dist = arrayfun(@(i) pdist([trx(male_id).x_mm(time_window(i)), trx(male_id).y_mm(time_window(i));
                                            trx(f_id).x_mm(time_window(i)), trx(f_id).y_mm(time_window(i))]), 1:length(time_window));
                bin = (dist < dist_thr) & (angle_deg < angle_threshold);

                if any(bin)
                    [bouts, lens] = detect_binarybouts(bin);
                    for i = 1:length(lens)-1
                        if bouts(i+1,1) - bouts(i,2) < join_threshold
                            bin(bouts(i,2):bouts(i+1,1)) = 1;
                        end
                    end
                    [bouts, lens] = detect_binarybouts(bin);

                    for i = 1:length(lens)
                        if lens(i) < time_thr
                            bin(bouts(i,1):bouts(i,2)) = 0;
                        end
                    end
                end

                inds = find(bin == 1);
                plot(timestamps(time_window(inds)), f * ones(size(inds)), '|', 'MarkerSize', 50);
            end

            ylim([0.5, length(female_ids)+0.5]);
            yticks(1:length(female_ids));
            yticklabels(arrayfun(@(i) sprintf('F%d', i), 1:length(female_ids), 'UniformOutput', false));
            xlabel('Time (s)');
            ylabel('Target Female');
            title(sprintf('Chamber %d | T \geq %.1fs | D < %dmm', chamber_to_plot, time_thr/FPS, dist_thr));
            box off

            results_dir = 'ethogram_output';
            if ~exist(results_dir, 'dir'), mkdir(results_dir); end
            saveas(gcf, fullfile(results_dir, sprintf('ethogram_ch%d_T%ds_D%02dmm.png', chamber_to_plot, time_thr/FPS, dist_thr)));
            savefig(gcf, fullfile(results_dir, sprintf('ethogram_ch%d_T%ds_D%02dmm.fig', chamber_to_plot, time_thr/FPS, dist_thr)));
        end

        combo = combo + 1;
    end
end

fprintf('Done analyzing all threshold combinations.\n');