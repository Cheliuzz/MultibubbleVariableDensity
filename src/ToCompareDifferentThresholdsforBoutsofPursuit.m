% To compare multiple threshold combinations and compute bout statistics including transitions
cd('/Volumes/otopaliklab/flydisco_data/2025-06-26/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_6x1_1hr_CSMH_20250626T081139')
load('registered_trx.mat')
load('movie-track.mat')
angle_struct = load('perframe/anglefrom1to2_nose2ell.mat');
angle_data = angle_struct.data;

%% Parameters
fly_IDs = [... % same fly_IDs matrix as before
    1, 2, 3, 4, 5, 6, 7;
    16,15,17,18,19,20,21;
    23,25,24,22,26,27,28;
    30,29,31,32,33,34,35;
    37,36,38,39,40,41,42;
    52,50,51,53,54,55,56;
    59,58,63,60,61,62,57;];

FPS = 60;
angle_threshold = 25; % degrees
chamber_to_plot = 1;

for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);

%% Threshold ranges
time_thresholds = [5, 8, 10] * FPS;
dist_thresholds = [3, 5, 8, 10];
summary_all = struct();
combo = 1;

for tt = 1:length(time_thresholds)
    for dt = 1:length(dist_thresholds)

        time_thr = time_thresholds(tt);
        dist_thr = dist_thresholds(dt);

        target_counts_all = zeros(size(fly_IDs,1), size(fly_IDs,2)-1);
        total_bouts_all = zeros(size(fly_IDs,1), 1);
        transitions_all = zeros(size(fly_IDs,1), 1);
        bout_durations_chamber = cell(1, size(fly_IDs,1));

        for chamber = 1:size(fly_IDs,1)
            male_id = fly_IDs(chamber, 1);
            f_ids = fly_IDs(chamber, 2:end);
            angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

            dists = nan(length(f_ids), endframe_all);
            valid_mask = false(length(f_ids), endframe_all);

            for f = 1:length(f_ids)
                f_id = f_ids(f);
                dists(f,:) = arrayfun(@(i) pdist([trx(male_id).x_mm(i), trx(male_id).y_mm(i); trx(f_id).x_mm(i), trx(f_id).y_mm(i)]), 1:endframe_all);
                valid_mask(f,:) = (dists(f,:) < dist_thr) & (angle_deg < angle_threshold);
            end

            [~, closest_f] = min(dists, [], 1);
            closest_valid = false(1, endframe_all);
            for i = 1:endframe_all
                if valid_mask(closest_f(i), i)
                    closest_valid(i) = true;
                end
            end

            target_per_frame = nan(1, endframe_all);
            [bouts, lens] = detect_binarybouts(closest_valid);
            for i = 1:length(lens)-1
                if bouts(i+1,1) - bouts(i,2) < 3
                    closest_valid(bouts(i,2):bouts(i+1,1)) = 1;
                end
            end
            [bouts, lens] = detect_binarybouts(closest_valid);

            durations_this_chamber = [];
            engaged_females = [];
            for i = 1:length(lens)
                if lens(i) >= time_thr
                    durations_this_chamber(end+1) = lens(i) / FPS;
                    current_target = mode(closest_f(bouts(i,1):bouts(i,2)));
                    engaged_females(end+1) = current_target;
                    target_per_frame(bouts(i,1):bouts(i,2)) = closest_f(bouts(i,1):bouts(i,2));
                end
            end

            transitions_all(chamber) = sum(diff(engaged_females) ~= 0);
            for f = 1:length(f_ids)
                target_counts_all(chamber, f) = sum(target_per_frame == f);
            end
            total_bouts_all(chamber) = numel(durations_this_chamber);
            bout_durations_chamber{chamber} = durations_this_chamber;
        end

        all_durations = [bout_durations_chamber{:}];

        summary_all(combo).time_thr = time_thr / FPS;
        summary_all(combo).dist_thr = dist_thr;
        summary_all(combo).target_counts = target_counts_all;
        summary_all(combo).total_bouts = total_bouts_all;
        summary_all(combo).mean_bout_duration = mean(all_durations);
        summary_all(combo).all_durations = bout_durations_chamber;
        summary_all(combo).transitions = transitions_all;

        combo = combo + 1;
    end
end

% Build per-chamber summary table
rows = sum(cellfun(@(x) size(x,2), {summary_all.total_bouts}));
bout_durations_chamber_table = table('Size', [rows, 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'Time_Threshold_s','Dist_Threshold_mm','Total_Bouts','Mean_Bout_Duration_s','Transitions'});

row_counter = 1;
for i = 1:length(summary_all)
    t_thr = summary_all(i).time_thr;
    d_thr = summary_all(i).dist_thr;
    n_chambers = length(summary_all(i).total_bouts);
    for c = 1:n_chambers
        bout_durations_chamber_table.Time_Threshold_s(row_counter) = t_thr;
        bout_durations_chamber_table.Dist_Threshold_mm(row_counter) = d_thr;
        bout_durations_chamber_table.Total_Bouts(row_counter) = summary_all(i).total_bouts(c);
        bout_durations_chamber_table.Transitions(row_counter) = summary_all(i).transitions(c);
        chamber_durations = summary_all(i).all_durations{c};
        if ~isempty(chamber_durations)
            bout_durations_chamber_table.Mean_Bout_Duration_s(row_counter) = mean(chamber_durations);
        else
            bout_durations_chamber_table.Mean_Bout_Duration_s(row_counter) = NaN;
        end
        row_counter = row_counter + 1;
    end
end