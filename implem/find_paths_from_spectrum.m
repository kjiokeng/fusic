%% Find paths from music spectrum
% Pmusic           -- the music spectrum
% theta            -- the values of theta (aoas) at which the spectrum was evaluated
% tau              -- the values of tau (tofs) at which the spectrum was evaluated
% n_paths          -- the number of paths to find

% Return:
% paths            -- the identified most dominant paths
%                     It is a matrix of size n_paths*3
%                     where each row represents a path in the form [aoa, tof, power]
function [paths] = find_paths_from_spectrum(Pmusic, theta, tau, n_paths)
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        pkg load image;
    end

    % Minimum value of the spectrum
    min_pmusic_value = min(min(Pmusic));

    % Save only the peaks
    binary_peaks_pmusic = imregionalmax(Pmusic);
    Pmusic = Pmusic .* binary_peaks_pmusic;
    Pmusic(Pmusic==0) = -Inf;

    % Find n_paths greatest values and their indices
    [sorted_peaks, ind] = sort(Pmusic(:), 'descend');  % Sort the values in descending order
    dominant_peaks_values = sorted_peaks(1:n_paths);
    dominant_peaks_ind = ind(1:n_paths);
    [dominant_peaks_ind_x dominant_peaks_ind_y] = ind2sub(size(Pmusic), dominant_peaks_ind);

    % Get the dominant peaks only
    dominant_peaks = min_pmusic_value * ones(size(Pmusic));
    dominant_peaks(dominant_peaks_ind_x, dominant_peaks_ind_y) = Pmusic(dominant_peaks_ind_x, dominant_peaks_ind_y);

    if false
        figure('Name', 'BINARY Peaks over AoA & ToF MUSIC Spectrum', 'NumberTitle', 'off')
        mesh(tau, theta, double(dominant_peaks))
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Spectrum Peaks')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end


    paths_aoas = theta(dominant_peaks_ind_x);
    paths_tofs = tau(dominant_peaks_ind_y);
    
    paths = zeros(n_paths, 3);
    paths(:,1) = paths_aoas;
    paths(:,2) = paths_tofs;
    paths(:,3) = dominant_peaks_values';
end