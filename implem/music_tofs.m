function [music_spectrum, tofs, power_db, tof_range] = music_tofs(csi_matrix, freq, freq_delta, n_sig)

	% Default value for number of signals
	if nargin < 4
		n_antennas = size(csi_matrix, 1);
        n_sig = n_antennas - 1;
    end

	% Number of subcarriers
	n_subcarriers = size(csi_matrix, 2);

	% Covariance matrix
	csi = csi_matrix.'; % Transpose without conjugate
	covariance_matrix = csi * ctranspose(csi);

	% Eigen decomposition of the covariance matrix
	[V, D] = eig(covariance_matrix);

	% Sorting the eigen vectors
	[d,ind] = sort(diag(D));
	Vs = V(:,ind);
	En = Vs(:,1:size(Vs, 2) - n_sig);
	En_squared = En * ctranspose(En);

	% error("----------------------")

	% Compute the music spectrum
	tof_range = (0:1.0:167) * 10^-9;
	distance_range = 0:0.05:100;
	speed_of_light = 3e8;
	tof_range = distance_range / speed_of_light;
	music_spectrum = zeros(length(tof_range), 1);
	for tof_ind = 1:length(tof_range)
		sv = steering_vector(tof_range(tof_ind), freq_delta, n_subcarriers);
		P = sv' * En_squared * sv;
		music_spectrum(tof_ind) = abs(1 / P);
	end

	% Convert to decibels
    % music_spectrum = 10 * log10(music_spectrum);

	% Finding tofs where peaks occur
	isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        pkg load signal;
    end
	[power_db, loc] = findpeaks(music_spectrum);
	tofs = tof_range(loc);

	if true
		% Conserve only peaks that have at least 20% of the strongest peak power
		max_peak = max(power_db);
		noise_level = min(power_db);
		valid_peaks_idx = (power_db - noise_level) > 0.05 * (max_peak - noise_level);
		% valid_peaks_idx = power_db > 6.3e-9;
		power_db = power_db(valid_peaks_idx);
		loc = loc(valid_peaks_idx);

		tofs = tof_range(loc);
	end

	if false
		figure('Name', 'ToF MUSIC spectrum');
		plot(tof_range, music_spectrum);
		xlabel('Time of Flight (ns)');
        ylabel('Power (db)');
        title('ToF Estimation from MUSIC Algorithm');
        grid on
	end
end

function a = steering_vector(tof, freq_delta, n_subcarriers)
	mult_factor = exp(-2i * pi * tof * freq_delta);

    a = ones(n_subcarriers, 1);
    for i=2:n_subcarriers
    	a(i) = a(i-1) * mult_factor;
    end
end

function power_db = to_db(power)
	power_db = 10*log10(power);
end