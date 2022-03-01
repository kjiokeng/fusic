%% Identifies and correct the most dominant paths followed by the signal
% This algorithm first of all performs 2D MUSIC algoritm to find the most dominant paths that the signal followed,
% then correct the computed values with an accurate distance measurement

function [corrected_paths] = find_and_correct_paths(csi, accurate_tof, ...
	antenna_distance, frequency, sub_freq_delta)
	
	% Useful variables
	n_antennas = size(csi, 1);
	n_subcarriers = size(csi, 2);
	num_computed_paths = n_antennas;

	% Compute the music spectrum with spotfi method
	[spectrum, aoas, tofs, num_computed_paths] = spotfi_music_spectrum(csi, ...
	        antenna_distance, frequency, sub_freq_delta);

	% Compute the music spectrum with pila method
	% [spectrum, aoas, tofs, num_computed_paths] = pila_music_spectrum(csi, ...
 %        antenna_distance, frequency, sub_freq_delta);

	% Find most dominants paths from music spectrum
	n_paths = num_computed_paths;
	paths = find_paths_from_spectrum(spectrum, aoas, tofs, n_paths);

	% Find path that has been used : the path that has effectively been followed by the signal
	used_path = find_used_path(paths);

	% The result
	corrected_paths = paths;

	% Do not perform the correction if tof is negative (erronous)
	if accurate_tof > 0
		% Compute the tof error of music (due to packet detection delay)
		used_path_tof = used_path(2);
		music_tof_error = used_path_tof - accurate_tof;
		
		% Correct all the distance measurements
		corrected_paths(:,2) = corrected_paths(:,2) - music_tof_error;
	end
	
end