% Close all opened windows
% close all;

% Indoor office
antenna_distance = 0.029;
frequency = 5190 * 10^6;
sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz
bw = 40 * 10^6;

% Read csi from linux-80211n-csitool output file
addpath('../');
addpath('../../linux-80211n-csitool-supplementary/matlab/');
csi_file = '../../ftm/2nd/ftm-csi/14-06-19/5m-0d-ref-7m-obs.2.dat';

c = 3 * 10^8;

N_ifft_bins = 2^16;
pdp = zeros(1, N_ifft_bins);
n_peaks = 2;
peaks_pw = zeros(1, n_peaks);
peaks_time = zeros(1, n_peaks);


% dir_name = '10-07-19/';
dir_name = '08-07-19-terrass/';
% directory = strcat('test-data/csi-files/', dir_name);
directory = strcat('../ftm/csi/', dir_name);
csi_files = dir(directory);
csi_files = csi_files(3:end); % We ignore '.' and '..' files
num_files = length(csi_files);
mkdir(strcat('results/', dir_name));

used_antenna = 1;
n_sig = 2;

for k=1:num_files
	close all;
	figure;
	subplot(2, 1, 1)
	hold on;

	csi_file = csi_files(k).name;
	% csi_file = '../ftm/csi/08-07-19-terrass/10m-0d-ref-_8.3m.dat'
	csi_file = '../ftm/csi/08-07-19-terrass/10m-0d-ref-_8.3m-2-pers.dat'
	% csi_trace = read_bf_file(strcat(directory, csi_file));
	csi_trace = read_bf_file(csi_file);
	
	if length(csi_trace) < 10
		continue
	end

	start_ind=min(100, length(csi_trace)-10);
	num_packets=min(100, length(csi_trace)-start_ind);
	num_packets_to_show=2;
	noise_levels = zeros(num_packets, 1);
	for ii=start_ind:start_ind+num_packets-1
		idx = ii - start_ind + 1;
		fprintf("File nb: %d/%d, %s, packet nb: %d/%d\n", k, num_files, csi_file, idx, num_packets);

		csi_entry = csi_trace{ii};
		csi = get_scaled_csi(csi_entry);
		csi = squeeze(csi(1,:,:));

		% Generate csi matrix
		% aoas = [45, -20, -80];
		% tofs = [4e-8, 8e-8, 12e-8, 16e-8, 24e-8, 32e-8, 36e-8, 40e-8];
		% n_antennas = 3;	
		% n_subcarriers = 30;
		% csi = generate_csi(aoas, tofs, n_antennas, antenna_distance, n_subcarriers, frequency, sub_freq_delta);

		% [packet_pdp, idx_values, packet_peaks_pw, packet_peaks_time, used_antenna] = power_delay_profile(csi, sub_freq_delta, N_ifft_bins);
		[packet_pdp, packet_peaks_time, packet_peaks_pw, idx_values] = music_tofs(csi, frequency, sub_freq_delta, n_sig);
		% spotfi_music_spectrum(csi, antenna_distance, frequency, sub_freq_delta);
		if length(packet_peaks_pw) < n_peaks
			continue
		end

		% Keep only peaks that are at least 10% of the max peak
		% noise_level = min(packet_pdp);
		% max_peak = max(packet_peaks_pw);
		% peaks_pw_zero = packet_peaks_pw - noise_level;
		% valid_peaks_locs = packet_peaks_pw > 0.1 * max_peak;
		% packet_peaks_pw = packet_peaks_pw(valid_peaks_locs);
		% packet_peaks_time = packet_peaks_time(valid_peaks_locs)


		noise_level = min(packet_pdp);
		pdp(idx, :) = packet_pdp(used_antenna,:);
		peaks_pw(idx, 1:n_peaks) = packet_peaks_pw(1:n_peaks);
		peaks_time(idx, 1:n_peaks) = packet_peaks_time(1:n_peaks);
		noise_levels(idx, 1) = min(packet_pdp);

		hold off
		% plot(idx_values * c, packet_pdp(1,:), 'LineWidth', 1.8);
		% plot(idx_values * c, 10*log10(packet_pdp/length(pdp)), 'LineWidth', 1.8);
		plot(idx_values * c, packet_pdp-noise_level, 'LineWidth', 1.8);
		rssi = 10*log10(sum(packet_pdp)/length(pdp))
		pause

		if idx <= num_packets_to_show
			% plot(idx_values * c, packet_pdp(1,:), 'LineWidth', 1.8);
			% plot(idx_values, 10 * log10(packet_pdp(1,:)), 'LineWidth', 1.8);
			plot(idx_values * c, packet_pdp, 'LineWidth', 1.8);
		end
	end
	title("Power Delay Profile")
	xlabel("Distance (m)")
	ylabel("Power")
	legend("Packet 1", "Packet 2")

	% Convert power to dBm
	% peaks_pw = 10 * log10(peaks_pw * 10^3);
	subplot(2, 1, 2)
	% plot(repmat(5,size(peaks_time,1),1), '-.k', 'LineWidth', 1.5)
	% hold on
	plot(peaks_time * c, 'LineWidth', 1.8);
	delta_diff = (peaks_time(1,2)-peaks_time(1,1)) * c
	relative_strenght = (peaks_pw(:,1) - noise_levels) ./ (peaks_pw(:,2) - noise_levels);
	pause
	[relative_strenght, 10*log10(relative_strenght)];
	%peaks_pw_db = 10 * log10(peaks_pw)
	hold off
	% legend("Ground truth", "Direct path", "Reflected path", 'Location','southeast')
	legend("Direct path", "Reflected path", 'Location','southeast')
	title("Resolved distance")
	xlabel("Packet index")
	ylabel("Distance (in m)")


	th1 = 1.6;
	th2 = 0.4;
	t = 0 / c;
	corrected_dist = zeros(length(relative_strenght), 1);
	
	% 6m-ref-_8.3m
	tftm = 13 / c;
	delta_t = 11.65 / c;

	% 4m-ref-_8.3m
	% tftm = 12 / c;
	% delta_t = 13.07 / c;
	for l=1:length(relative_strenght)
		delta_t = (peaks_time(l,2) - peaks_time(l,1));
		R = relative_strenght(l);
		if R > th1
			disp("case R > th1")
			t = tftm;
		elseif R < th2
			disp("case R < th2")
			t = tftm - delta_t;
		else
			disp("case default")
			% t = tftm - delta_t * 0.8 * exp(10 * log10(R));
			% t = tftm - delta_t * R^2;
			t = tftm - delta_t * 0.8 * R ;
		end

		corrected_dist(l) = t * c;
	end

	corrected_dist
	figure
	hold on
	plot(6 * ones(1,length(corrected_dist)), 'LineWidth', 1.8);
	plot(tftm * c * ones(1,length(corrected_dist)), 'LineWidth', 1.8);
	plot(corrected_dist, 'LineWidth', 1.8);
	legend("Ground truth", "Measured", "Corrected")
	xlabel("Packet index")
	ylabel("Distance (m)")

	pause

	str = input('Enter title: ','s');
	title({csi_file, str})

	if strcmp(str, '')
		continue;
	end

	image_name = strcat('results/', dir_name, csi_file, '.png');
	saveas(gcf, image_name);
	close all;
end

% return


% x = idx_values;
% y = start_ind:start_ind+num_packets-1;
% z = pdp;
% surf(idx_values, start_ind:start_ind+num_packets-1, pdp);
% colormap jet;
% pcolor(pdp);
% colormap(gray(200))