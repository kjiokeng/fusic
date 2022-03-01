function [corrected_distance, std_dev, corrected_dist] =  correct_ftm(ftm_meas, csi_file, ground_truth)
	% Read csi from linux-80211n-csitool output file
	addpath('../');
	addpath('../../linux-80211n-csitool-supplementary/matlab/');

	if nargin < 3
		ground_truth = -1;
	end

	% Indoor office
	antenna_distance = 0.029;
	frequency = 5190 * 10^6;
	sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz
	bw = 40 * 10^6;
	c = 3 * 10^8;

	N_ifft_bins = 2^16;
	pdp = zeros(1, N_ifft_bins);
	n_peaks = 2;
	peaks_pw = zeros(1, n_peaks);
	peaks_time = zeros(1, n_peaks);

	used_antenna = 1;
	n_sig = 2;


	csi_trace = read_bf_file(csi_file);

	start_ind=min(1, length(csi_trace));
	num_packets=min(1000, length(csi_trace)-start_ind);
	n_treated_packets=0;
	num_packets_to_show=2;
	noise_levels = zeros(num_packets, 1);
	for ii=start_ind:start_ind+num_packets-1
		idx = ii - start_ind + 1;
		% fprintf("File: %s, packet nb: %d/%d\n", csi_file, idx, num_packets);

		try
			csi_entry = csi_trace{ii};
			csi = get_scaled_csi(csi_entry);
			csi = squeeze(csi(1,:,:));
		catch
			continue
		end

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
		rssi = 10*log10(sum(packet_pdp)/length(pdp));
		pause

		if idx <= num_packets_to_show
			% plot(idx_values * c, packet_pdp(1,:), 'LineWidth', 1.8);
			% plot(idx_values, 10 * log10(packet_pdp(1,:)), 'LineWidth', 1.8);
			figure
			plot(idx_values * c, packet_pdp, 'LineWidth', 1.8);
		end
		n_treated_packets = n_treated_packets + 1;
	end

	% No valid packet, stop
	if n_treated_packets < 1
		corrected_distance = ftm_meas * 3/4;
		std_dev = 0;
		return
	end
	% title("Power Delay Profile")
	% xlabel("Distance (m)")
	% ylabel("Power")
	% legend("Packet 1", "Packet 2")

	% % Convert power to dBm
	% % peaks_pw = 10 * log10(peaks_pw * 10^3);
	% subplot(2, 1, 2)
	% % plot(repmat(5,size(peaks_time,1),1), '-.k', 'LineWidth', 1.5)
	% % hold on
	% plot(peaks_time * c, 'LineWidth', 1.8);
	% delta_diff = (peaks_time(1,2)-peaks_time(1,1)) * c
	% relative_strenght = (peaks_pw(:,1) - noise_levels) ./ (peaks_pw(:,2) - noise_levels);
	% hold off
	% % legend("Ground truth", "Direct path", "Reflected path", 'Location','southeast')
	% legend("Direct path", "Reflected path", 'Location','southeast')
	% title("Resolved distance")
	% xlabel("Packet index")
	% ylabel("Distance (in m)")


	% th1 = 0.42;
	th1 = 0.7;

	ftm_meas = ftm_meas - 0.7; % Offset correction
	tftm = ftm_meas / c;
	t = 0;
	corrected_dist = zeros(length(peaks_pw), 1);
	r = zeros(length(peaks_pw), 1);

	for l=1:length(peaks_pw)
		% delta_t = (peaks_time(l,2) - peaks_time(l,1));
		pw = peaks_pw(l,:);
		pw = peaks_pw(l,:) - noise_levels(l);
		% pw = pw .* pw;
		sum_pw = sum(pw);
		% R = peaks_pw(l,1) * peaks_pw(l,1) / sum(pw .* pw);
		% R = peaks_pw(l,1) / sum(pw);
		R = peaks_pw(l,1) / sum(pw(1:end));
		r(l) = R;
		% R = peaks_pw(l,1) / (sum_pw - peaks_pw(l,1));

		% R = relative_strenght(l);
		if R > th1
			% disp("case R > th1")
			t = tftm;
		else
			% disp("case default")
			% t = tftm - delta_t * 0.8 * exp(10 * log10(R));
			% t = tftm - delta_t * R^2;
			% t = tftm - delta_t * 0.8 * R;

			delta_t = peaks_time(l,:) - peaks_time(l,1);
			% delta_t(2) = 4.42/c;

			mean_excess_delay = sum(pw .* delta_t) / sum_pw;
			t = tftm - mean_excess_delay;
			% t = max(t, tftm * 2/3);
		end

		corrected_dist(l) = t * c;
	end

	corrected_dist = corrected_dist(corrected_dist>0.5);

	plot(r)
	r = r(~isnan(r));
	m = mean(r)
	m_dB= 10*log10(m)

	[m, m_dB]

	% corrected_dist
	if false
		plot(r)
		figure
		hold on
		plot(ftm_meas * ones(1,length(corrected_dist)), 'LineWidth', 1.8);
		plot(corrected_dist, 'LineWidth', 1.8);
		if ground_truth ~= -1
			plot(ground_truth * ones(1,length(corrected_dist)), 'LineWidth', 1.8);
		end
		legend("Measured", "Corrected", "Ground truth")
		xlabel("Packet index")
		ylabel("Distance (m)")
	end

	% corrected_dist(100) = 8
	% ssss = mean(corrected_dist)
	% ssss = std(corrected_dist)

	corrected_dist = corrected_dist(~isnan(corrected_dist));

	% k = 5;
	% vals = zeros(length(corrected_dist)/k, 1);
	% for l=1:length(vals)
	% 	vals(l) = mean(corrected_dist(k*(l-1)+1:k*l));
	% end

	% Result
	std_dev = std(corrected_dist);
	corrected_distance = mean(corrected_dist);
	corrected_distance = r;
	% corrected_dist = corrected_dist;
end