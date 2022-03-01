function [spectrum, theta_range, tof_range] = music_spectrum(csi_matrix, n_sig, theta_range, tof_range)

	% Defaults values
	if nargin < 2
		n_sig = 3;
	end
	if nargin < 3
		theta_range = [-90:45:90] * pi / 180;
	end
	if nargin < 4
		max_range = 51; % Max range in meters
		spatial_resolution = 5; % Spatial resolution in meters 
		speed_of_light = 299792458; % Speed of light in m/s
		tof_range = [0:spatial_resolution:max_range] / speed_of_light;
	end


	% Set physical layer parameters (frequency, subfrequency spacing, and antenna spacing
    antenna_distance = 0.1;
    % frequency = 5 * 10^9;
    frequency = 5.785 * 10^9;
    % frequency = 5.32 * 10^9;
    sub_freq_delta = (40 * 10^6) / 30;


    % Sanitize ToFs with Algorithm 1
    csi_matrix = spotfi_algorithm_1(csi_matrix, sub_freq_delta);
    % Acquire smoothed CSI matrix
    csi_matrix = smooth_csi(csi_matrix);


	% Number of subcarriers
	n_antennas = size(csi_matrix, 1);
	n_subcarriers = size(csi_matrix, 2);
	
	% Compute the covariance matrix XX'
	covariance_matrix = csi_matrix * csi_matrix';
	
	% Eigen decomposition of the covariance matrix
	[V, D] = eig(covariance_matrix);

	% Sorting the eigen vectors
	[d,ind] = sort(diag(D));
	Vs = V(:,ind);
	En = Vs(:,1:n_antennas - n_sig);
	En_squared = En * En';

	% Evaluate MUSIC spectrum over the specified grid
	n_theta = length(theta_range);
	n_tof = length(tof_range);
	spectrum = zeros(n_theta, n_tof);
	for theta_ind = 1:n_theta
		for tof_ind = 1:n_tof
			steering_vector = compute_steering_vector(theta_range(theta_ind), tof_range(tof_ind), ...
                    frequency, sub_freq_delta, antenna_distance);
			spectrum(theta_ind, tof_ind) = 1 / (steering_vector' * En_squared * steering_vector);
		end
	end
end


%% Computes the steering vector for SpotFi. 
% Each steering vector covers 2 antennas on 15 subcarriers each.
% theta           -- the angle of arrival (AoA) in degrees
% tau             -- the time of flight (ToF)
% freq            -- the central frequency of the signal
% sub_freq_delta  -- the frequency difference between subcarrier
% ant_dist        -- the distance between each antenna
% Return:
% steering_vector -- the steering vector evaluated at theta and tau
%
% NOTE: All distance measurements are in meters
function steering_vector = compute_steering_vector(theta, tau, freq, sub_freq_delta, ant_dist)
    steering_vector = zeros(30, 1);
    k = 1;
    base_element = 1;
    for ii = 1:2
        for jj = 1:15
            steering_vector(k, 1) = base_element * omega_tof_phase(tau, sub_freq_delta)^(jj - 1);
            k = k + 1;
        end
        base_element = base_element * phi_aoa_phase(theta, freq, ant_dist);
    end
end

%% Compute the phase shifts across subcarriers as a function of ToF
% tau             -- the time of flight (ToF)
% frequency_delta -- the frequency difference between adjacent subcarriers
% Return:
% time_phase      -- complex exponential representing the phase shift from time of flight
function time_phase = omega_tof_phase(tau, sub_freq_delta)
    time_phase = exp(-1i * 2 * pi * sub_freq_delta * tau);
end

%% Compute the phase shifts across the antennas as a function of AoA
% theta       -- the angle of arrival (AoA) in degrees
% frequency   -- the frequency of the signal being used
% d           -- the spacing between antenna elements
% Return:
% angle_phase -- complex exponential representing the phase shift from angle of arrival
function angle_phase = phi_aoa_phase(theta, frequency, d)
    % Speed of light (in m/s)
    c = 3.0 * 10^8;
    % Convert to radians
    theta = theta / 180 * pi;
    angle_phase = exp(-1i * 2 * pi * d * sin(theta) * (frequency / c));
end

%% Creates the smoothed CSI matrix by rearranging the various csi values in the default CSI matrix.
% csi          -- the regular CSI matrix to use for creating the smoothed CSI matrix
% Return:
% smoothed_csi -- smoothed CSI matrix following the construction put forth in the SpotFi paper.
%                   Each column in the matrix includes data from 2 antennas and 15 subcarriers each.
%                   Has dimension 30x32. 
function smoothed_csi = smooth_csi(csi)
    smoothed_csi = zeros(size(csi, 2), size(csi, 2));
    % Antenna 1 (values go in the upper left quadrant)
    m = 1;
    for ii = 1:1:15
        n = 1;
        for j = ii:1:(ii + 15)
            smoothed_csi(m, n) = csi(1, j); % 1 + sqrt(-1) * j;
            n = n + 1;
        end
        m = m + 1;
    end
    
    % Antenna 2
    % Antenna 2 has its values in the top right and bottom left
    % quadrants, the first for loop handles the bottom left, the second for
    % loop handles the top right
    
    % Bottom left of smoothed csi matrix
    for ii = 1:1:15
        n = 1;
        for j = ii:1:(ii + 15)
            smoothed_csi(m, n) = csi(2, j); % 2 + sqrt(-1) * j;
            n = n + 1;
        end
        m = m + 1;
    end
    
    % Top right of smoothed csi matrix
    m = 1;
    for ii = 1:1:15
        n = 17;
        for j = ii:1:(ii + 15)
            smoothed_csi(m, n) = csi(2, j); %2 + sqrt(-1) * j;
            n = n + 1;
        end
        m = m + 1;
    end
    
    % Antenna 3 (values go in the lower right quadrant)
    for ii = 1:1:15
        n = 17;
        for j = ii:1:(ii + 15)
            smoothed_csi(m, n) = csi(3, j); %3 + sqrt(-1) * j;
            n = n + 1;
        end
        m = m + 1;
    end
end

%% Time of Flight (ToF) Sanitization Algorithm, find a linear fit for the unwrapped CSI phase
% csi_matrix -- the CSI matrix whose phase is to be adjusted
% delta_f    -- the difference in frequency between subcarriers
% Return:
% csi_matrix -- the same CSI matrix with modified phase
function [csi_matrix, phase_matrix] = spotfi_algorithm_1(csi_matrix, delta_f, packet_one_phase_matrix)
    %% Time of Flight (ToF) Sanitization Algorithm
    %  Obtain a linear fit for the phase
    %  Using the expression:
    %      argmin{\rho} \sum_{m,n = 1}^{M, N} (\phi_{i}(m, n) 
    %          + 2 *\pi * f_{\delta} * (n - 1) * \rho + \beta)^2
    %
    %  Arguments:
    %  M is the number of antennas
    %  N is the number of subcarriers
    %  \phi_{i}(m, n) is the phase for the nth subcarrier, 
    %      on the mth antenna, for the ith packet
    %  f_{\delta} is the frequency difference between the adjacent
    %      subcarriers
    %  \rho and \beta are the linear fit variables
    %
    % Unwrap phase from CSI matrix
    R = abs(csi_matrix);
    phase_matrix = unwrap(angle(csi_matrix), pi, 2);

    % Parse input args
    if nargin < 3
        packet_one_phase_matrix = phase_matrix;
    end

    % STO is the same across subcarriers....
    % Data points are:
    % subcarrier_index -> unwrapped phase on antenna_1
    % subcarrier_index -> unwrapped phase on antenna_2
    % subcarrier_index -> unwrapped phase on antenna_3
    fit_X(1:30, 1) = 1:1:30;
    fit_X(31:60, 1) = 1:1:30;
    fit_X(61:90, 1) = 1:1:30;
    fit_Y = zeros(90, 1);
    for i = 1:size(phase_matrix, 1)
        for j = 1:size(phase_matrix, 2)
            fit_Y((i - 1) * 30 + j) = packet_one_phase_matrix(i, j);
        end
    end

    % Linear fit is common across all antennas
    result = polyfit(fit_X, fit_Y, 1);
    tau = result(1);

    for m = 1:size(phase_matrix, 1)
        for n = 1:size(phase_matrix, 2)
            % Subtract the phase added from sampling time offset (STO)
            %phase_matrix(m, n) = packet_one_phase_matrix(m, n) + (2 * pi * delta_f * (n - 1) * tau);
            phase_matrix(m, n) = packet_one_phase_matrix(m, n) - (n - 1) * tau;
        end
    end
    
    % Reconstruct the CSI matrix with the adjusted phase
    csi_matrix = R .* exp(1i * phase_matrix);
end