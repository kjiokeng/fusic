% Close all opened windows
% close all;

% Indoor office
antenna_distance = 0.029;
frequency = 5190 * 10^6;
sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz
bw = 40 * 10^6;

% Read csi from linux-80211n-csitool output file
addpath('../linux-80211n-csitool-supplementary/matlab/');
% csi_file = 'test-data/csi-files/5ghz/27-02-19/csi-_2m-2m.dat';
% csi_file = 'test-data/csi-files/5ghz/27-02-19/csi-0m-8m.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/_2m-2m.2.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/1m-0d.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/4m-los_45d-nlos.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/22-05-19/a0m45.3.dat';
% csi_file = 'test-data/csi-files/29-05-19/csi-4m-0d.dat';
% csi_file = 'test-data/csi-files/29-05-19/csi-_2m-4m.dat';
% csi_file = 'test-data/csi-files/29-05-19/csi-7m-0d-ref-45d.dat';
% csi_file = 'test-data/csi-files/03-06-19/office/4m-0d-ref-4m.dat';
% csi_file = 'test-data/csi-files/04-06-19/5m-0d-ref-+_3.5m.dat';
% csi_file = 'test-data/csi-files/04-06-19/3m-0d-ref-3m.dat';
% csi_file = 'test-data/csi-files/04-06-19/7.50m-0d-ref-+_6.75m.dat';
% csi_file = 'test-data/csi-files/calibration/csi-calibration.3.dat';
% csi_file = 'test-data/csi-files/office/test.2.dat';
% csi_file = 'test-data/csi-files/office/65v.dat';

% csi_file = '../ftm/csi/5m-0d-ref-6.5m.3.dat';
csi_file = '../ftm/csi/14-06-19/5m-0d-ref-7m-noobs.dat';
% csi_file = '../ftm/csi/14-06-19/5m-0d-ref-7m-obs.2.dat';
% csi_file = '../ftm/csi/14-06-19/5m-0d-ref-7m-pcasobs.dat';


c = 3 * 10^8;

N_ifft_bins = 2^16;
pdp = zeros(1, N_ifft_bins);
n_peaks = 2;
peaks_pw = zeros(1, n_peaks);
peaks_time = zeros(1, n_peaks);


csi_trace = read_bf_file(csi_file);
start_ind=20;
num_packets=length(csi_trace);
num_packets=100;
for ii=start_ind:start_ind+num_packets-1
	idx = ii - start_ind + 1;
	csi_entry = csi_trace{ii};

	csi = get_scaled_csi(csi_entry);
	csi = squeeze(csi(1,:,:));

	% Generate csi matrix
	% aoas = [45, -20, -80];
	% tofs = [4e-8, 8e-8, 12e-8, 16e-8, 24e-8, 32e-8, 36e-8, 40e-8];
	% n_antennas = 3;	
	% n_subcarriers = 30;
	% csi = generate_csi(aoas, tofs, n_antennas, antenna_distance, n_subcarriers, frequency, sub_freq_delta);

	[packet_pdp, idx_values, packet_peaks_pw, packet_peaks_time] = power_delay_profile(csi, sub_freq_delta, N_ifft_bins);
	pdp(idx, :) = packet_pdp(1,:);
	peaks_pw(idx, 1:n_peaks) = packet_peaks_pw(1, 1:n_peaks);
	peaks_time(idx, 1:n_peaks) = packet_peaks_time(1, 1:n_peaks);


	plot(idx_values, packet_pdp(1,:), 'LineWidth', 1.8);
	% % stem(idx_values, packet_pdp(1,:));
	% % hist(pdp(1,:), 30);
	% legend("Antenna 1", "Antenna 2", "Antenna 3");
	xlabel("Time (s)")
	ylabel("Power")
	title("More obstacles: 10m")

	% return
	
end

% Convert power to dBm
peaks_pw = 10 * log10(peaks_pw * 10^3);

plot(repmat(5,size(peaks_time,1),1), '-.k', 'LineWidth', 1.8)
hold on
plot(peaks_time * c, 'LineWidth', 1.8);
mean(peaks_time(:,2)-peaks_time(:,1)) * c
hold off
legend("Ground truth", "Direct path", "Reflected path")
title("Resolved distance of propagation paths")
xlabel("Packet index")
ylabel("Distance (in m)")
set(gca, 'linewidth', 2)

% x = idx_values;
% y = start_ind:start_ind+num_packets-1;
% z = pdp;
% surf(idx_values, start_ind:start_ind+num_packets-1, pdp);
% colormap jet;
% pcolor(pdp);
% colormap(gray(200))