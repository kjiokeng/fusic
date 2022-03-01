function valid = is_valid_csi(csi_entry)
	rssi_threshold = 5;
	valid = 1;

	if (csi_entry.rssi_a < rssi_threshold) || (csi_entry.rssi_b < rssi_threshold) || (csi_entry.rssi_c < rssi_threshold)
		fprintf("Too weak signal on one or more antennas\n")
		valid = 0;
		return
	end

	csi = get_scaled_csi(csi_entry);
	csi = squeeze(csi);
	csi_entry;
	csi = csi(csi_entry.perm,:);

	angles = angle(csi)';
	% angles = unwrap(angles);
	angles = rad2deg(angles');
	diffs(1,:) = angles(2,:) - angles(1,:);
	diffs(2,:) = angles(3,:) - angles(1,:);

	ind2 = abs(diffs(1,:))<90;
	ind3 = abs(diffs(2,:))<90;

	vals2 = diffs(1,:);
	vals2 = vals2(ind2);

	vals3 = diffs(2,:);
	vals3 = vals3(ind3);
	% diffs = diffs(diffs>-100);

	delta2 = mean(vals2)
	delta3 = mean(vals3)
	
	subplot(2,1,1)
	plot(angles')
	title("Unwrapped phase of antennas")
	legend("Ant 1", "Ant 2", "Ant 3")


	subplot(2,1,2)
	plot(diffs', '-x')
	hold on
	plot(vals2, '-s')
	plot(vals3, '-s')
	hold off
	title("Relative phase between antennas")
	legend("Ant 2-1", "Ant 3-1")
	xlabel("Subcarrier index")
	pause

	valid = 0;

end