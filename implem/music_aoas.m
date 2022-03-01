function [music_spectrum, aoas, power_db] = music_aoas(csi, freq, antenna_distance, n_sig)

	% cov_mat = 0;
	% csi_trace = csi;
	% for ii=length(csi_trace):-1:1
	% 	csi_entry = csi_trace{ii};
	% 	csi = get_scaled_csi(csi_entry);
	% 	csi = squeeze(csi);

	% 	covariance_matrix = csi * csi';
	% 	cov_mat = cov_mat + covariance_matrix;
	% end

	% covariance_matrix = cov_mat / length(csi_trace);



	% Number of antennas
	n_antennas = size(csi, 1);
	
	% Default value for number of signals
	if nargin < 4
        n_sig = n_antennas - 1;
    end

	% Covariance matrix
	covariance_matrix = csi * ctranspose(csi);

	% Eigen decomposition of the covariance matrix
	[V, D] = eig(covariance_matrix);

	% Sorting the eigen vectors
	[d,ind] = sort(diag(D));
	Vs = V(:,ind);
	En = Vs(:,1:n_antennas - n_sig);
	En_squared = En * ctranspose(En);

	%error("----------------------")

	% Compute the music spectrum
	theta_range = [-90:1:90];
	music_spectrum = zeros(length(theta_range), 1);
	for theta_ind = 1:length(theta_range)
		sv = steering_vector(theta_range(theta_ind), freq, antenna_distance, n_antennas);
		P = sv' * En_squared * sv;
		music_spectrum(theta_ind) = abs(1 / P);
	end

	% Convert to decibels
    % music_spectrum = 10 * log10(music_spectrum);

	% Finding aoas where peaks occur
	isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        pkg load signal;
    end
	[power_db, loc] = findpeaks(music_spectrum);
	aoas = theta_range(loc);

	if true
		figure('Name', 'AoA MUSIC spectrum');
		plot(theta_range, music_spectrum);
		xlabel('Angle of Arival (°)');
        ylabel('Power (db)');
        title('AoA Estimation from MUSIC Algorithm');
        grid on
	end
end

function a = steering_vector(theta, freq, antenna_distance, n_antennas)
	% Speed of light (in m/s)
    c = 3e8;

    mult_factor = exp(-2i * pi * antenna_distance * sind(theta) * freq / c);

    a = ones(n_antennas, 1);
    for i=2:n_antennas
    	a(i) = a(i-1) * mult_factor;
    end
end

function power_db = to_db(power)
	power_db = 10*log10(power);
end