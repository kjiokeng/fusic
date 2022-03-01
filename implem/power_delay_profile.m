function [pdp, idx_values, peaks_vals, peaks_idx, used_antenna] = power_delay_profile(csi, sub_freq_delta, N)
	n_subcarriers = size(csi, 2);
	time_res = 1 / (N * sub_freq_delta);
	used_antenna = 2;

	% hamming_window = hamming(n_subcarriers)';
	% csi = csi .* repmat(hamming_window, size(csi, 1), 1);
	
	pdp = double(abs(ifft(csi, N, 2)).^2);
	idx_values = (1:N) * time_res;

	[peaks_vals, peaks_loc] = findpeaks(pdp(used_antenna,:));

	% Conserve only peaks that have at least 10% of the strongest peak power
	max_peak = max(peaks_vals);
	valid_peaks_idx = peaks_vals > 0.2 * max_peak;
	% valid_peaks_idx = peaks_vals > 6.3e-9;
	peaks_vals = peaks_vals(valid_peaks_idx);
	peaks_loc = peaks_loc(valid_peaks_idx);

	peaks_idx = idx_values(peaks_loc);
end