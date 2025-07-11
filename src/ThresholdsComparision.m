% Example: Visualize angle at one frame
frame = 10000;  % or pick any i
male_id = 16;    % choose a male
f_id = 17;       % and one female to visualize

% Get positions
xm = trx(male_id).x_mm(frame);
ym = trx(male_id).y_mm(frame);
xf = trx(f_id).x_mm(frame);
yf = trx(f_id).y_mm(frame);

% Get orientation (theta) of male in radians
theta = trx(male_id).theta(frame);

% Plot the male fly
figure(100); clf
plot(xm, ym, 'ko', 'MarkerFaceColor', 'k'); hold on

% Plot male orientation vector
arrow_length = 5;  % mm
quiver(xm, ym, cos(theta)*arrow_length, sin(theta)*arrow_length, 0, 'k', 'LineWidth', 2)

% Plot the angle cone
angle_threshold = 25;  % degrees
cone_angles = linspace(-angle_threshold, angle_threshold, 100);  % degrees
cone_radians = deg2rad(cone_angles);
cone_x = cos(cone_radians + theta) * arrow_length + xm;
cone_y = sin(cone_radians + theta) * arrow_length + ym;
fill([xm cone_x], [ym cone_y], [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

% Plot female position
plot(xf, yf, 'ro', 'MarkerFaceColor', 'r')

% Draw a line from male to female
plot([xm xf], [ym yf], 'r--')

axis equal
xlabel('x (mm)')
ylabel('y (mm)')
title(sprintf('Visual cone: %d° threshold', angle_threshold))
legend({'Male', 'Orientation', 'Cone', 'Female'}, 'Location', 'Best')


%% Compare filtered vs. unfiltered bins for a sample male & female
sample_chamber = 5;
male_id = fly_IDs(sample_chamber, 1);
female_id = fly_IDs(sample_chamber, 2);

% Distances and angles
dist = arrayfun(@(i) pdist([
    trx(male_id).x_mm(i), trx(male_id).y_mm(i);
    trx(female_id).x_mm(i), trx(female_id).y_mm(i)
]), 1:endframe_all);
angle_deg = abs(rad2deg(angle_data{male_id}(1:endframe_all)));

% Original bin
raw_bin = (dist < dist_threshold) & (angle_deg < angle_threshold);

% Filtered bin (using same logic as above)
bin = raw_bin;

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

% Plot comparison
figure;
subplot(2,1,1);
plot(raw_bin, 'b'); ylim([-0.1, 1.1])
title('Raw Bouts')
xlabel('Frame'); ylabel('Raw State')

subplot(2,1,2);
plot(bin, 'r'); ylim([-0.1, 1.1])
title('Filtered Bouts')
xlabel('Frame'); ylabel('Filtered State')