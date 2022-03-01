function [m, v] = evaluate_spectrum(filepath, expected_angles)
	
	% Anechoic chamber
	antenna_distance = 0.028;
	frequency = 5180 * 10^6;
	sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz

	% Read csi from linux-80211n-csitool output file
	addpath('../linux-80211n-csitool-supplementary/matlab/');
	csi_trace = read_bf_file(filepath);

	for ii=30:length(csi_trace)-30
		csi_entry = csi_trace{ii};
		valid_csi = is_valid_csi(csi_entry);
		if ~valid_csi
			continue;
		end

		csi = get_scaled_csi(csi_entry);
		csi = squeeze(csi);
		csi = csi(csi_entry.perm,:);

		accurate_tof = -3e-8;
		corrected_paths = find_and_correct_paths(csi, accurate_tof, ...
			antenna_distance, frequency, sub_freq_delta)

		disp('############# Estimated distances : (in meters) #############');
		[corrected_paths(:,1), 3e8 * corrected_paths(:,2)]
	end
end