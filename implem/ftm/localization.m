addpath('../');
addpath('../../linux-80211n-csitool-supplementary/matlab/');

% summary_file = "../../ftm/distance/22-07-19-sti/summary-sti.csv";
% summary_file = "../../ftm/distance/22-07-19-crous/summary-crous.csv";
% summary_file = "../../ftm/distance/24-07-19-salle-de-convivialite/summary-salle-de-convivialite.csv";

% summary_file = "../../ftm/distance/22-07-19-sti/summary-sti-corrected.csv";
% summary_file = "../../ftm/distance/22-07-19-crous/summary-crous-corrected.csv";
% summary_file = "../../ftm/distance/24-07-19-salle-de-convivialite/summary-salle-de-convivialite-corrected.csv";

% summary_file = "../../ftm/distance/22-07-19-sti/summary-sti-corrected-2.csv";
% summary_file = "../../ftm/distance/22-07-19-crous/summary-crous-corrected-2.csv";
% summary_file = "../../ftm/distance/24-07-19-salle-de-convivialite/summary-salle-de-convivialite-corrected-2.csv";

% summary_file = "../../ftm/distance/22-07-19-sti/summary-sti-corrected-3.csv";
% summary_file = "../../ftm/distance/22-07-19-crous/summary-crous-corrected-3.csv";
% summary_file = "../../ftm/distance/24-07-19-salle-de-convivialite/summary-salle-de-convivialite-corrected-3.csv";

% summary_file = "../../ftm/distance/22-07-19-sti/summary-sti-corrected-4.csv";
% summary_file = "../../ftm/distance/22-07-19-crous/summary-crous-corrected-4.csv";
summary_file = "../../ftm/distance/24-07-19-salle-de-convivialite/summary-salle-de-convivialite-corrected-4.csv";


% Read the data, ignoring the headers
data = csvread(summary_file, 1);

% AP coordinates
aps_coordinates = data(1:3,2:3);
center = mean(aps_coordinates);

for k=4:length(data)
	ground_truth_pos = data(k,2:3);
	distances_to_aps = data(k,4:6)-0.7; % Offset correction

	% Just for testing. Give the exact distances
	% distances_to_aps(1) = norm(ground_truth_pos-aps_coordinates(1,:));
	% distances_to_aps(2) = norm(ground_truth_pos-aps_coordinates(2,:));
	% distances_to_aps(3) = norm(ground_truth_pos-aps_coordinates(3,:));

	[x, y] = localize(aps_coordinates, distances_to_aps, center);
	computed_pos = [x, y];

	err = norm(ground_truth_pos - computed_pos);
	data(k,7:9) = [computed_pos err];
end


% Localization error statistics
disp("########## Localization error statistics ##########")
loc_err = data(4:end,end);
min_loc_err = min(loc_err)
median_loc_err = median(loc_err)
mean_loc_err = mean(loc_err)
max_loc_err = max(loc_err)


% Distance estimation statistics
disp("########## Distance estimation statistics ##########")
dist_err = zeros(length(data)-3, 3);
for k=4:length(data)
	ground_truth_pos = data(k,2:3);
	distances_to_aps = data(k,4:6)-0.7; % Offset correction

	gt_distances_to_aps(1) = norm(ground_truth_pos-aps_coordinates(1,:));
	gt_distances_to_aps(2) = norm(ground_truth_pos-aps_coordinates(2,:));
	gt_distances_to_aps(3) = norm(ground_truth_pos-aps_coordinates(3,:));

	dist_diffs = abs(distances_to_aps-gt_distances_to_aps);
	dist_err(k-3,:) = dist_diffs;
	data(k,10:12) = gt_distances_to_aps;
	data(k,13:15) = dist_diffs;
end

min_dist_err = min(dist_err)
median_dist_err = median(dist_err)
mean_dist_err = mean(dist_err)
max_dist_err = max(dist_err)

% Distance estimation statistics
disp("########## Overall distance estimation statistics ##########")
dist_err = reshape(dist_err, 3*length(dist_err), 1);
min_dist_err = min(dist_err)
median_dist_err = median(dist_err)
mean_dist_err = mean(dist_err)
max_dist_err = max(dist_err)


disp("########## All the data ##########")
headers = ['pos,x_ground_truth,y_ground_truth,ap1_meas_dist,ap2_meas_dist,ap3_meas_dist,' ...
	'computed_pos_x,computed_pos_y,loc_err,ground_truth_dist_ap1,ground_truth_dist_ap2,ground_truth_dist_ap3', ...
	'dist_est_err_ap1,dist_est_err_ap2,dist_est_err_ap3'];
disp(headers)
disp(data)


% Plot the locations
if false
	close all
	scatter(data(4:end,2), data(4:end,3), 'filled', 'blue')
	hold on
	scatter(data(4:end,7), data(4:end,8), 'filled', 'red')
	scatter(aps_coordinates(:,1), aps_coordinates(:,2), 200, 'filled', '^k')
	xlabel("x (m)")
	ylabel("y (m)")
	legend("Ground truth", "Computed positions", "APs")

	dx = 0.1;
	dy = 0.1;
	c = cellstr(num2str(data(4:end,1)));
	text(data(4:end,2)+dx, data(4:end,3)+dy, c, 'Fontsize', 10);
	text(data(4:end,7)+dx, data(4:end,8)+dy, c, 'Fontsize', 10);
end



function [x, y] = localize(aps_coordinates, distances_to_aps, x0)
	objfun = @(x) sum_of_squared_diffs(aps_coordinates, distances_to_aps, x);
	
	[pos,fval] = fminunc(objfun, x0);

	x = pos(1);
	y = pos(2);
end

function obj = sum_of_squared_diffs(aps_coordinates, distances_to_aps, x)
	obj = 0;
	for ap_idx=1:length(aps_coordinates)
		dif = norm(x-aps_coordinates(ap_idx,:)) - distances_to_aps(ap_idx);
		obj = obj + dif * dif;
	end
end