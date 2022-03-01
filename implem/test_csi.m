% Close all opened windows
close all;

% Indoor office
antenna_distance = 0.029;
frequency = 5190 * 10^6;
sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz

% Read csi from linux-80211n-csitool output file
addpath('../linux-80211n-csitool-supplementary/matlab/');
% csi_file = 'test-data/csi-files/calibration/csi-calibration.3.dat';
csi_file = 'test-data/csi-files/25-06-19/test.dat';

csi_trace = read_bf_file(csi_file);

start_ind=1;
num_packets=200;
num_packets=min(num_packets, length(csi_trace)-start_ind+1);
for ii=start_ind:start_ind+num_packets-1
	csi_entry = csi_trace{ii};
	csi = get_scaled_csi(csi_entry);
	csi = squeeze(csi(1,:,:));

	plot(rad2deg(angle(csi))', '--x');

	pause
end