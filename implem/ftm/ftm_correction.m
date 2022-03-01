addpath('../');
addpath('../../linux-80211n-csitool-supplementary/matlab/');

% csi_directory = "../../ftm/csi/22-07-19-sti/";
% summary_file = "../../ftm/distance/22-07-19-sti/summary-sti.csv";

% csi_directory = "../../ftm/csi/22-07-19-crous/";
% summary_file = "../../ftm/distance/22-07-19-crous/summary-crous.csv";

csi_directory = "../../ftm/csi/24-07-19-salle-de-convivialite/";
summary_file = "../../ftm/distance/24-07-19-salle-de-convivialite/summary-salle-de-convivialite.csv";

% Read the data, ignoring the headers
data = csvread(summary_file, 1);

% AP coordinates
aps_coordinates = data(1:3,2:3)
center = mean(aps_coordinates);

% The resulting matrix
corrected_dist_to_aps = zeros(length(data), 3);

% Correct each distance estimation
for k=4:length(data)
	pos = data(k,1);
	ground_truth_pos = data(k,2:3);
	distances_to_aps = data(k,4:6)-0.7; % Offset correction

	gt_distances_to_aps(1) = norm(ground_truth_pos-aps_coordinates(1,:));
	gt_distances_to_aps(2) = norm(ground_truth_pos-aps_coordinates(2,:));
	gt_distances_to_aps(3) = norm(ground_truth_pos-aps_coordinates(3,:));

	csi_file = strcat(csi_directory, 'ap1/', 'pos', num2str(pos), '.dat');
	[corrected_dist, std_dev] = correct_ftm(distances_to_aps(1), csi_file, gt_distances_to_aps(1));
	corrected_dist_to_aps(k,1) = corrected_dist;
	% pause

	csi_file = strcat(csi_directory, 'ap2/', 'pos', num2str(pos), '.dat');
	[corrected_dist, std_dev] = correct_ftm(distances_to_aps(2), csi_file, gt_distances_to_aps(2));
	corrected_dist_to_aps(k,2) = corrected_dist;
	% pause

	csi_file = strcat(csi_directory, 'ap3/', 'pos', num2str(pos), '.dat');
	[corrected_dist, std_dev] = correct_ftm(distances_to_aps(3), csi_file, gt_distances_to_aps(3));
	corrected_dist_to_aps(k,3) = corrected_dist;

	% pause
	close all;
end

% Replace the non corrected distances by the corrected ones
non_corrected_dist_to_aps = data(:,4:6);
data(:,4:6) = corrected_dist_to_aps;

% Display all
disp("########## All the data ##########")
headers = ['pos,x_ground_truth,y_ground_truth,ap1_meas_dist,ap2_meas_dist,ap3_meas_dist'];
disp(headers)
disp(data)