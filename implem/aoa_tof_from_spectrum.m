%% Estimate aoas and tofs from precomputed music spectrum
% Pmusic           -- the music spectrum
% theta            -- the values of theta (aoas) at which the spectrum was evaluated
% tau              -- the values of tau (tofs) at which the spectrum was evaluated
% Return:
% estimated_aoas   -- the angle of arrivals that gave peaks from running MUSIC, as a vector
% estimated_tofs   -- the time of flights that gave peaks on the estimated_aoas from running music.
%                         This is a matrix with dimensions [length(estimated_aoas, ), length(tau)].
%                         The columns are zero padded at the ends to handle different peak counts 
%                           across different AoAs.
%                         I.E. if there are three AoAs then there will be three rows in 
%                           estimated_tofs
function [estimated_aoas, estimated_tofs] = aoa_tof_from_spectrum(Pmusic, theta, tau)
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        pkg load image;
    end
    binary_peaks_pmusic = imregionalmax(Pmusic);
    
    % Get AoAs that have peaks
    % fprintf('Future estimated aoas\n')
    aoa_indices = linspace(1, size(binary_peaks_pmusic, 1), size(binary_peaks_pmusic, 1));
    aoa_peaks_binary_vector = any(binary_peaks_pmusic, 2);
    estimated_aoas = theta(aoa_peaks_binary_vector);

    aoa_peak_indices = aoa_indices(aoa_peaks_binary_vector);
    
    % Get ToFs that have peaks
    time_peak_indices = zeros(length(aoa_peak_indices), length(tau));
    % AoA loop (only looping over peaks in AoA found above)
    for ii = 1:length(aoa_peak_indices)
        aoa_index = aoa_peak_indices(ii);
        binary_tof_peaks_vector = binary_peaks_pmusic(aoa_index, :);
        matching_tofs = tau(binary_tof_peaks_vector);
        
        % Pad ToF rows with -1s to have non-jagged matrix
        negative_ones_for_padding = -1 * ones(1, length(tau) - length(matching_tofs));
        time_peak_indices(ii, :) = horzcat(matching_tofs, negative_ones_for_padding);
    end

    
    if true
        figure('Name', 'BINARY Peaks over AoA & ToF MUSIC Spectrum', 'NumberTitle', 'off')
        mesh(tau, theta, double(binary_peaks_pmusic))
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Spectrum Peaks')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end

    if true
        % Theta (AoA) & Tau (ToF) 3D Plot
        figure('Name', 'Selective AoA & ToF MUSIC Peaks, with only peaked AoAs', 'NumberTitle', 'off')
        mesh(tau, estimated_aoas, Pmusic(aoa_peak_indices, :))
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Spectrum Peaks')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end
    
    if true
        % Tau (ToF)
        for ii = 1:length(estimated_aoas)
            figure_name_string = sprintf('ToF Estimation as a Function of Tau w/ AoA: %f', ...
                    estimated_aoas(ii));
            figure('Name', figure_name_string, 'NumberTitle', 'off')
            plot(tau, Pmusic(ii, :), '-k')
            xlabel('Time of Flight \tau / seconds')
            ylabel('Spectrum function P(\theta, \tau)  / dB')
            title(figure_name_string)
            grid on
        end
    end
    
    % Set return values
    % AoA is now a column vector
    estimated_aoas = transpose(estimated_aoas);
    % ToF is now a length(estimated_aoas) x length(tau) matrix, with -1 padding for unused cells
    estimated_tofs = time_peak_indices;
end