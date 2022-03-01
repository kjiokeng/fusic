% Close all opened windows
close all;

% Set physical layer parameters (frequency, subfrequency spacing, and antenna spacing
% antenna_distance = 0.025;
% % frequency = 5 * 10^9;
% frequency = 5.785 * 10^9;
% % frequency = 5.32 * 10^9;
% sub_freq_delta = (40 * 10^6) / 30;

% antenna_distance = 0.03;
% frequency = 5230 * 10^6;
% sub_freq_delta = (40 * 10^6) / 32;
% sub_freq_delta = 2 * 312.5 * 10^3; % For 20 MHz
% sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz
% sub_freq_delta = 40*4*10^6/114; % For 40 MHz


% Firmin
% frequency = 5230 * 10^6;
% antenna_distance = 3e8/frequency/2;
% antenna_distance = 0.025;
% sub_freq_delta = (40 * 10^6) / 32;
% sub_freq_delta = 2 * 312.5 * 10^3; % For 20 MHz
% sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz
% sub_freq_delta = 40*4*10^6/114; % For 40 MHz

% Anechoic chamber
antenna_distance = 0.03;
frequency = 5190 * 10^6;
sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz
% sub_freq_delta = 1.403508 * 10^6;

% Indoor office
antenna_distance = 0.029;
frequency = 5190 * 10^6;
sub_freq_delta = 4 * 312.5 * 10^3; % For 40 MHz


% Create a random csi matrix
%csi = (rand(3, 30) + i*rand(3,30)) / sqrt(2);

% Read csi from csv file
% csi=csvread("csi.txt");
% csi=csi(1:3,:);

% Read csi from linux-80211n-csitool output file
addpath('../linux-80211n-csitool-supplementary/matlab/');
% csi_file = 'test-data/csi-files/5ghz/27-02-19/csi-_2m-2m.dat';
% csi_file = 'test-data/csi-files/5ghz/27-02-19/csi-0m-8m.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/_2m-2m.2.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/1m-0d.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/4m-los_45d-nlos.dat';
% csi_file = 'test-data/csi-files/anechoic/5130MHz-40MHz/22-05-19/a0m45.3.dat';
% csi_file = 'test-data/csi-files/29-05-19/csi-4m-0d.dat';
% csi_file = 'test-data/csi-files/29-05-19/csi-_2m-4m.dat';
% csi_file = 'test-data/csi-files/29-05-19/csi-7m-0d-ref-45d.dat';
% csi_file = 'test-data/csi-files/03-06-19/office/4m-0d-ref-4m.dat';
% csi_file = 'test-data/csi-files/04-06-19/5m-0d-ref-+_3.5m.dat';
% csi_file = 'test-data/csi-files/04-06-19/3m-0d-ref-3m.dat';
% csi_file = 'test-data/csi-files/04-06-19/7.50m-0d-ref-+_6.75m.dat';
% csi_file = 'test-data/csi-files/calibration/csi-calibration.3.dat';
% csi_file = 'test-data/csi-files/office/test.2.dat';
% csi_file = 'test-data/csi-files/office/65v.dat';

% Expected results for dir_name='04-06-19/'
dir_name = '04-06-19/';
expected_results = {
					[0 -60],
					[0 -34 55],
					[0 63],
					[0 -54 54],
					[0 63 -38.5],
					[0 -57],
					[0 45],
					[0 -38.5],
					[0 -61 61],
					[0 45],
					[0 -56]
				};

% Expected results for dir_name='06-06-19/'
% dir_name = '06-06-19/';
% expected_results = {
% 					[25 54],
% 					[0 38],
% 					[0 62],
% 					[0 43 -43],
% 					[0 51 -32],
% 					[9 -47 38],
% 					[-25 -38 47],
% 					[0]
% 				};

% Expected results for dir_name='29-05-19'
% expected_results = {
% 					[0],
% 					[0],
% 					[0 -45],
% 					[0 -27],
% 					[0 -39],
% 					[0 -51]
% 				};

% Expected results for dir_name='14-06-19'
% expected_results = {
% 					[0],
% 					[0 -70],
% 					[0 -70],
% 					[0 -70]
% 				};


% Expected results for dir_name='24-06-19/'
% dir_name = '24-06-19/';
% situations = [
% 	0 0.5 -2 0;
% 	0 1.5 -2 0;
% 	0 1 -2 0;
% 	0 2.5 -2 0;
% 	0 2 -2 0;
% 	0 3.5 -1.5 0;
% 	0 3.5 -1 0;
% 	0 3.5 -2.5 0;
% 	0 3.5 -2 0;
% 	0 3.5 -2 0;
% 	0 3 -2 0;
% 	-0.5 0.5 -2 0;
% 	-0.5 1.5 -2 0;
% 	-0.5 1 -2 0;
% 	-0.5 2.5 -2 0;
% 	-0.5 2 -2 0;
% 	-0.5 3 -2 0;
% 	-1.5 0.5 -2 0;
% 	-1.5 1.5 -2 0;
% 	-1.5 1 -2 0;
% 	-1.5 2.5 -2 0;
% 	-1.5 2 -2 0;
% 	-1.5 3 -2 0;
% 	-1.5 6 -2 0;
% 	-1.5 6 -2 0;
% 	-1 0.5 -2 0;
% 	-1 1.5 -2 0;
% 	-1 1 -2 0;
% 	-1 2.5 -2 0;
% 	-1 2 -2 0;
% 	-1 3 -2 0;
% ];
% expected_results = {};
% for k=1:length(situations)
% 	posC = situations(k, 1:2);
% 	posMur = situations(k, 3:4);
% 	[~, theta0,theta1] = getPaths(posC, posMur);
% 	expected_results(end+1) = {[theta0, theta1]};
% end


% Expected results for dir_name='25-06-19/'
dir_name = '25-06-19/';
situations = [
	0 10 -6.5 0 7.3 0;
	0 11 -6.5 0 7.3 0;
	0 12 -6.5 0 7.3 0;
	0 13 -6.5 0 7.3 0;
	0 14 -6.5 0 7.3 0;
	0 15 -6.5 0 7.3 0;
	0 1 -6.5 0 7.3 0;
	0 2 -6.5 0 7.3 0;
	2 3 -6.5 0 7.3 0;
	0 3 -6.5 0 7.3 0;
	0 4 -6.5 0 7.3 0;
	4 3 -6.5 0 7.3 0;
	0 5 -6.5 0 7.3 0;
	0 6 -6.5 0 7.3 0;
	0 7 -6.5 0 7.3 0;
	0 8 -6.5 0 7.3 0;
	0 9 -6.5 0 7.3 0;
	-2 3 -6.5 0 7.3 0;
	-3 3 -6.5 0 7.3 0;
	-4 3 -6.5 0 7.3 0;
	-5 3 -6.5 0 7.3 0;
	0 3 -6.5 0 7.3 0;
];


expected_results = {};
for k=1:length(situations)
	posC = situations(k, 1:2);
	posMur = situations(k, 3:4);
	[~, theta0,theta1] = getPaths(posC, posMur);

	posMur = situations(k, 5:6);
	[~, ~,theta11] = getPaths(posC, posMur);

	expected_results(end+1) = {[theta0, theta1, theta11]};
end


% directory = strcat('../ftm/2nd/ftm-csi/', dir_name);
directory = strcat('test-data/csi-files/', dir_name);
csi_files = dir(directory);
csi_files = csi_files(3:end); % We ignore '.' and '..' files

num_files = length(csi_files)
mkdir(strcat('results/', dir_name));
mkdir(strcat('results/with_calib', dir_name));

tic
for calib=1:2
	for k=1:num_files
		figure;
		hold on;

		csi_file = csi_files(k).name;
		% fprintf("%s\n", csi_file);
		% continue
		csi_trace = read_bf_file(strcat(directory, csi_file));

		% For parallel computing
		D = parallel.pool.DataQueue;
		D.afterEach(@(val) scatter(val(:,1), val(:,2), 'x'));

		start_ind=1;
		num_packets=200;
		num_packets=min(num_packets, length(csi_trace)-start_ind+1);
		parfor ii=start_ind:start_ind+num_packets-1
			packet_nb = ii - start_ind + 1;
			fprintf("File nb: %d/%d, %s, packet nb: %d/%d\n", k, num_files, csi_file, packet_nb, num_packets);
			% val = val + 1;
			csi_entry = csi_trace{ii};
			% valid_csi = is_valid_csi(csi_entry);
			% % valid_csi = 1;
			% if ~valid_csi
			% 	continue;
			% end

			csi = get_scaled_csi(csi_entry);
			csi = squeeze(csi(1,:,:));

			% csi(csi_entry.perm) = csi(1:3,:);

			% Calibration for card 1 (long card)
			% csi(2,:) = csi(2,:) * exp(j * deg2rad(-80)); 
			% csi(3,:) = csi(3,:) * exp(j * deg2rad(-20)); 

			% Calibration for card 2 (short card)
			% csi(2,:) = csi(2,:) * exp(j * deg2rad(-93)); 
			% csi(3,:) = csi(3,:) * exp(j * deg2rad(-197));

			% Calibration for card 2 (short card) (5th floor)
			% csi(2,:) = csi(2,:) * exp(j * deg2rad(-49)); 
			% csi(3,:) = csi(3,:) * exp(j * deg2rad(109)); 

			% Calibration for card 2 (short card) (wired channel)
			if calib==2
				% csi(2,:) = csi(2,:) * exp(j * deg2rad(-12)); 
				csi(2,:) = csi(2,:) * exp(j * deg2rad(-24)); 
				csi(3,:) = csi(3,:) * exp(j * deg2rad(-67)); 
			end

			% SpotFi's phase sanitization
			% csi = spotfi_algorithm_1(csi, sub_freq_delta);

			% covmat = csi * csi';
			% % mat = mat + covmat;
			% [doas,spec,specang] = musicdoa(covmat,1);
			% doas
			% plot(specang,spec)
			% hold on
			% xlabel('Arrival Angle (deg)')
			% ylabel('Magnitude (dB)')
			% title('MUSIC Spectrum')
			% pause
			% continue
			

			% Generate csi matrix
			% aoas = [45, -20, -80];
			% aoas = [45, -20, -80];
			% tofs = [4e-8, 8e-8, 12e-8, 16e-8, 24e-8, 32e-8, 36e-8, 40e-8];
			% n_antennas = 3;
			% n_subcarriers = 30;
			% csi = generate_csi(aoas, tofs, n_antennas, antenna_distance, n_subcarriers, frequency, sub_freq_delta);

			
			% csi = spotfi_algorithm_1(csi, sub_freq_delta);

			% angles = angle(csi);
			% angles = rad2deg(angle(csi));
			% figure('Name', "Phase of CSI", 'NumberTitle', 'off')
			% plot(angles');
			% xlabel('Subcarrier index')
			% ylabel('Phase / deg')
			% title('Phase of CSI')
			% legend("Antenna 1", "Antenna 2", "Antenna 3")
			% pause
			% error('###########CSI');


			accurate_tof = -3e-8;
			corrected_paths = find_and_correct_paths(csi, accurate_tof, ...
				antenna_distance, frequency, sub_freq_delta);

			% disp('############# Estimated distances : (in meters) #############');
			val = [corrected_paths(:,1), 3e8 * corrected_paths(:,2)];			% res(packet_nb) = 1;
			% res{packet_nb} = val;
			% scatter(val(:,1), val(:,2), 'x');
			send(D, val);

			% [music_spectrum, aoas, power_db] = music_aoas(csi, frequency, antenna_distance);
			% [aoas', power_db]
			% [music_spectrum, aoas, power_db] = music_aoas(csi_trace, frequency, antenna_distance);

			% [music_spectrum, tofs, power_db] = music_tofs(csi, frequency, sub_freq_delta);

			% fprintf('############ ii =%d/%d\n', ii, length(csi_trace))
			% pause
			% error('Description');
			% close all;
		end

		% Plot the results as a scatter plot		
		% for packet_nb=1:num_packets
		% 	r = res(packet_nb);
		% 	scatter(r(:,1), r(:,2), 'x');
		% end

		% Add expected results as vertical lines
		if k <= length(expected_results)
			expected_angles = expected_results{k};
			for ang_ind=1:length(expected_angles)
				exp_value = expected_angles(ang_ind);
				plot([exp_value exp_value], ylim, '--r');
			end
		end

		% Add title and labels, and save to file
		xlabel('AoA (in °)');
		ylabel('Distance (in m)')
		xlim([-90 90])
		set(gca,'xtick',-90:10:90);

		image_title = csi_file;
		image_name = strcat('results/', dir_name, csi_file, '.png');
		if calib==2
			image_title = strcat(csi_file, ' with calib');
			image_name = strcat('results/', dir_name, 'with_calib/', csi_file, '-calib', '.png');
		end

		title(image_title);
		saveas(gcf, image_name);
		close all;
		fprintf('##############################################################################\n');
		fprintf('####### Saved superresolution results to %s ########\n', image_name);
		fprintf('##############################################################################\n\n\n');
	end
end
toc

% mat = mat / val;
% [doas,spec,specang] = musicdoa(mat,2);
% plot(specang,spec)
% xlabel('Arrival Angle (deg)')
% ylabel('Magnitude (dB)')
% title('MUSIC Spectrum')
% doas