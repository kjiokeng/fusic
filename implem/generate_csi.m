function csi = generate_csi(aoas, tofs, n_antennas, antenna_distance, n_subcarriers, freq, sub_freq_delta)
	% csi matrix initialized as empty
	% csi = zeros(n_antennas, n_subcarriers);
	csi = ones(n_antennas, n_subcarriers);

	% Speed of light
	% c = 3e8;
	c = 299792458;

	% Number of signals to be generated
	n_sig = length(aoas);

	for i=1:n_sig
		aoa = aoas(i);
		tof = tofs(i);

		% Generate csi matrix for this aoa and tof
		csi_of_path = ones(n_antennas, n_subcarriers);

		% Phase shift due to aoa
		aoa_phase_shift = exp(-2i * pi * antenna_distance * sind(aoa) * freq / c);
		for j=2:n_antennas
			csi_of_path(j,1) = csi_of_path(j-1,1) * aoa_phase_shift;
		end

		% Phase shift due to tof
		tof_phase_shift = exp(-2i * pi * sub_freq_delta * tof);
		for j=2:n_subcarriers
			csi_of_path(:,j) = csi_of_path(:,j-1) * tof_phase_shift;
		end

		csi = csi + csi_of_path;
		% csi = csi .* csi_of_path;
	end
end