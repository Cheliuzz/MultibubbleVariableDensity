%% Script to ID Pursuit Bouts and Durations with Angle Filter
% Experiment: 6 females per male per chamber

cd('/Volumes/otopaliklab/flydisco_data/2025-06-26/MultibubbleVariableDensity_multibubble__whiteOnly1hour062025_CSMH_6x1_1hr_CSMH_20250626T081139');
load('registered_trx.mat')
load('movie-track.mat')
angle_struct = load('perframe/anglefrom1to2_nose2ell.mat');
angle_data = angle_struct.data;

%% Parameters
fly_IDs = [
    1, 2, 3, 4, 5, 6, 7;
    14,8, 9,10,11,12,13;
    16,15,17,18,19,20,21;
    23,22,24,25,26,27,28;
    29,30,31,32,33,34,35;
    37,36,38,39,40,41,42;
    44,43,45,46,47,48,49;
    52,51,50,53,54,55,56;
    59,58,57,60,61,62,63;
];

FPS = 60;
time_threshold = 10 * FPS;
dist_threshold = 5;
angle_threshold = 30; % degrees
join_threshold = 3;

for i = 1:length(trx)
    nframes(i) = trx(i).endframe;
end
endframe_all = min(nframes);



%% Plot per-frame feature: anglefrom1to2_nose2ell (radians)
% For each male fly, plot angle to each of their targets over time


figure(101); clf
n_chambers = size(fly_IDs, 1);
fly_id_list = [trx.id];
time = trx(1).timestamps(1:endframe_all);  % Shared time base

for chamber = 1:n_chambers
    male_id = fly_IDs(chamber, 1);
    male_index = find(fly_id_list == male_id);

    if isempty(male_index)
        warning('Male ID %d not found in trx.', male_id);
        continue;
    end

    angle = angle_data{male_index}(1:endframe_all);

    subplot(n_chambers, 1, chamber);
    plot(time, rad2deg(angle), 'k', 'LineWidth', 1.2); hold on
    title(sprintf('Chamber %d | Male %d', chamber, male_id));
    ylabel('Angle to Closest Fly (deg)');
    ylim([-180, 180]);
    xlim([time(1), time(end)]);
end

xlabel('Time (s)');
sgtitle('Angle from Male to Closest Fly (nose2ell)');