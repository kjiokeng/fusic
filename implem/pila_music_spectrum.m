function [spectrum, aoas, tofs, n_paths] = pila_music_spectrum(csi, ...
        antenna_distance, frequency, sub_freq_delta, data_name)
    
    if nargin == 4
        data_name = '-';
    end

    % Useful variables
    n_antennas = size(csi, 1);
    n_subcarriers = size(csi, 2);
    n_paths = n_antennas;
    n_paths = 4;

    % Sanitize ToFs with Algorithm 1
    % sanitized_csi = spotfi_algorithm_1(csi, sub_freq_delta);
    % csi = sanitized_csi;


    % Construct the matrix H
    H = zeros(n_antennas * n_subcarriers, 1);
    for antenna_ind=1:n_antennas
        begin_index = (antenna_ind-1)*n_subcarriers+1;
        H(begin_index : begin_index + n_subcarriers - 1) = csi(antenna_ind,:);
    end

    % Compute the covariance matrix R
    R = H * ctranspose(H);

    % Spatial smoothing on R
    % N_sub1 = 2;
    % N_sub2 = 30;
    % L1 = 3 - N_sub1 + 1;
    % L2 = 30 - N_sub2 + 1;
    
    % R_smoothed = zeros(30, 30);
    % for m=1:L1
    %     for n=1:L2
    %         % row_start = m+n-1
    %         row_start = 30*(m-1)+1
    %         col_start = 30*(n-1)+1
    %         sub_R = R(row_start:row_start+29, col_start:col_start+29);
    %         R_smoothed = R_smoothed + sub_R;
    %     end
    % end
    % R_smoothed = R_smoothed / (L1 * L2);
    % R = R_smoothed

    [eigenvectors, eigenvalue_matrix] = eig(R);
    % Sorting the eigen values and eigen vectors
    [d,ind] = sort(diag(eigenvalue_matrix));
    eigenvalue_matrix = eigenvalue_matrix(ind, ind);
    eigenvectors = eigenvectors(:,ind);

    % Signal subspace
    Es = eigenvectors(:,size(eigenvectors,2)-n_paths+1:size(eigenvectors,2));
    % Noise subspace
    En = eigenvectors(:,1:size(eigenvectors,2)-n_paths);
    En_squared = En * ctranspose(En);

    % Music spectrum computation
    theta = -90:1:90; 
    % tau = 0:(2.0 * 10^-9):(334 * 10^-9);
    tau = 0:(2.0 * 10^-9):(334 * 10^-9);
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
    Pmusic = 10 * log10(Pmusic);

    if false
        % Theta (AoA) & Tau (ToF) 3D Plot
        figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off')
        mesh(tau, theta, Pmusic)
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Power (dB)')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end

    % if true
    %     % Theta (AoA)
    %     figure_name_string = sprintf('%s: Number of Paths: %d', data_name, n_paths);
    %     figure('Name', figure_name_string, 'NumberTitle', 'off')
    %     plot(theta, to_db(Pmusic(:, 1)), '-k')
    %     xlabel('Angle, \theta')
    %     ylabel('Spectrum function P(\theta, \tau)  / dB')
    %     title('AoA Estimation as a function of theta')
    %     grid on
    % end

    spectrum = Pmusic;
    aoas = theta;
    tofs = tau;
end

function steering_vector = compute_steering_vector(theta, tau, freq, sub_freq_delta, ant_dist, n_antennas, n_subcarriers)
    speed_of_light = 3e8;

    theta = theta * pi / 180;
    base_freq = freq - ((n_subcarriers+1)/2) * sub_freq_delta;

    steering_vector = zeros(n_antennas * n_subcarriers, 1);
    for antenna_ind=1:n_antennas
        begin_index = (antenna_ind-1)*n_subcarriers;
        for subcarrier_index = 1:n_subcarriers
            freq_i = base_freq + (subcarrier_index-1) * sub_freq_delta;
            lambda_i = speed_of_light / freq_i;
            delta_psi = 2*pi*(subcarrier_index-1)*sub_freq_delta*tau + 2*pi*ant_dist*(antenna_ind-1)*sin(theta)/lambda_i;

            steering_vector(begin_index+subcarrier_index) = exp(-1i * delta_psi);
        end
    end
end

function power_db = to_db(power)
    power_db = 10*log10(power);
end
