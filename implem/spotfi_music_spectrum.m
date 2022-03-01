%% Run MUSIC algorithm with SpotFi method including ToF and AoA
% x                -- the signal matrix
% antenna_distance -- the distance between the antennas in the linear array
% frequency        -- the frequency of the signal being localized
% sub_freq_delta   -- the difference between subcarrier frequencies
% data_name        -- the name of the data file being operated on, used for labeling figures
% Return:
% spectrum         -- the computed music spectrum
% aoas             -- the angle-of-arrivals at which the spectrum was evaluated
% tofs             -- the time-of-flights at which the spectrum was evaluated
function [spectrum, aoas, tofs, num_computed_paths] = spotfi_music_spectrum(csi, ...
        antenna_distance, frequency, sub_freq_delta, data_name)
    
    if nargin == 4
        data_name = '-';
    end

    % Useful variables
    n_antennas = size(csi, 1);
    n_subcarriers = size(csi, 2);

    % Sanitize ToFs with Algorithm 1
    % sanitized_csi = spotfi_algorithm_1(csi, sub_freq_delta);
    % csi = sanitized_csi;

    % Acquire smoothed CSI matrix
    smoothed_sanitized_csi = smooth_csi(csi);
    x = smoothed_sanitized_csi;
    
    % Data covarivance matrix
    R = x * x';
    R = R / (size(x, 1) * size(x, 2));
    % Find the eigenvalues and eigenvectors of the covariance matrix
    [eigenvectors, eigenvalue_matrix] = eig(R);
    % Find max eigenvalue for normalization
    max_eigenvalue = -1111;
    for ii = 1:size(eigenvalue_matrix, 1)
        if eigenvalue_matrix(ii, ii) > max_eigenvalue
            max_eigenvalue = eigenvalue_matrix(ii, ii);
        end
    end

    for ii = 1:size(eigenvalue_matrix, 1)
        eigenvalue_matrix(ii, ii) = eigenvalue_matrix(ii, ii) / max_eigenvalue;
    end
    
    % Sorting the eigen values and eigen vectors
    [d,ind] = sort(diag(eigenvalue_matrix));
    eigenvalue_matrix = eigenvalue_matrix(ind, ind);
    eigenvectors = eigenvectors(:,ind);

    % Find the largest decrease ratio that occurs between the last 10 elements (largest 10 elements)
    % and is not the first decrease (from the largest eigenvalue to the next largest)
    % Compute the decrease factors between each adjacent pair of elements, except the first decrease
    % start_index = size(eigenvalue_matrix, 1) - 2;
    start_index = size(eigenvalue_matrix, 1) - 1;
    end_index = start_index - 10;
    decrease_ratios = zeros(start_index - end_index + 1, 1);
    k = 1;
    for ii = start_index:-1:end_index
        temp_decrease_ratio = eigenvalue_matrix(ii + 1, ii + 1) / eigenvalue_matrix(ii, ii);
        decrease_ratios(k, 1) = temp_decrease_ratio;
        k = k + 1;
    end
    
    [max_decrease_ratio, max_decrease_ratio_index] = max(decrease_ratios);    
    index_in_eigenvalues = size(eigenvalue_matrix, 1) - max_decrease_ratio_index;
    num_computed_paths = size(eigenvalue_matrix, 1) - index_in_eigenvalues + 1;
    num_computed_paths = 5;
    
    % Estimate noise subspace
    column_indices = 1:(size(eigenvalue_matrix, 2) - num_computed_paths);
    eigenvectors = eigenvectors(:, column_indices);
    En_squared = eigenvectors * eigenvectors';
    
    % Peak search
    theta = -90:1:90;
    tau = 0:(1 * 10^-9):(334 * 10^-9);
    Pmusic = zeros(length(theta), length(tau));
    for ii = 1:length(theta)
        for jj = 1:length(tau)
            steering_vector = compute_steering_vector(theta(ii), tau(jj), ...
                    frequency, sub_freq_delta, antenna_distance, n_antennas, n_subcarriers);
            PP = steering_vector' * En_squared * steering_vector;
            Pmusic(ii, jj) = abs(1 /  PP);
        end
    end

    % Convert to decibels
    Pmusic = to_db(Pmusic);

    if true
        % Theta (AoA) & Tau (ToF) 3D Plot
        % figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off')
        mesh(tau, theta, Pmusic);
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Power (dB)')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end

    if false
        % Theta (AoA)
        figure_name_string = sprintf('%s: Number of Paths: %d', data_name, num_computed_paths);
        figure('Name', figure_name_string, 'NumberTitle', 'off')
        plot(theta, Pmusic(:, 1), '-k');
        xlabel('Angle, \theta')
        ylabel('Spectrum function P(\theta, \tau)  / dB')
        title('AoA Estimation as a function of theta')
        grid on
    end

    spectrum = Pmusic;
    aoas = theta;
    tofs = tau;
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
function steering_vector = compute_steering_vector(theta, tau, freq, sub_freq_delta, ant_dist, n_antennas, n_subcarriers)
    n_sensors_per_ant = n_subcarriers / 2;

    steering_vector = zeros((n_antennas-1) * n_sensors_per_ant, 1);
    k = 1;
    base_element = 1;
    for ii = 1:n_antennas-1
        for jj = 1:n_subcarriers/2
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
%                   Each column in the matrix includes data from 2 antennas and (n_subcarriers/2) subcarriers each.
function smoothed_csi = smooth_csi(csi)
    n_antennas = size(csi, 1);
    n_subcarriers = size(csi, 2);
    n_sensors_per_ant = n_subcarriers / 2;
    smoothed_csi = zeros((n_antennas-1) * n_sensors_per_ant, n_subcarriers);

    begin_ant_index = 1;
    begin_subcarrier_index = 1;
    for ii=1:size(smoothed_csi,1)
        smoothed_csi(ii, 1:n_sensors_per_ant) = csi(begin_ant_index, ...
            begin_subcarrier_index:begin_subcarrier_index+n_sensors_per_ant-1);
        smoothed_csi(ii, n_sensors_per_ant+1:n_subcarriers) = csi(begin_ant_index+1, ...
            begin_subcarrier_index:begin_subcarrier_index+n_sensors_per_ant-1);

        begin_subcarrier_index = begin_subcarrier_index + 1;
        if begin_subcarrier_index > n_sensors_per_ant
            begin_ant_index = begin_ant_index + 1;
            begin_subcarrier_index = 1;
        end
    end
end

function power_db = to_db(power)
    power_db = 10*log10(power);
end