#!/bin/bash
if [[ $# -ne 1 ]]; then
	echo "Usage : $0 <ground_truth_range (in cm)>"
	exit 1
fi

# The ground truth range
ground_truth_range=$1;

# Set the log file of the entire script
log_file="ftm-range-results-${ground_truth_range}cm"
exec 1>$log_file
# exec 2>&1 # Redirect stderr to stdout

# Config file name
config_file="conf"

# AP MAC addresses
mac_address_2g=04:92:26:89:2b:e0
mac_address_5g=04:92:26:89:2b:e4

# AP emitting center frequencies
cf_2g=2462
cf_5g=5500

# AP second center frequencies (for 40MHz and 80MHz bw)
cf1_2g_40m=2452
cf1_5g_40m=5510
cf1_5g_80m=5530

# Bandwidths
bws=(20 40 80)

# Number of samples per burst
spbs=(2 5 10)

# Number of time each measurement should be repeated
repeats=(1 5)

# The wlan interface name
wlan_inteface="wlp1s0"



### Execution
echo "#### Ground truth : $ground_truth_range cm"

# The command to launch the measurement
command="sudo iw $wlan_inteface measurement ftm_request $config_file"

# Loop over the number of repetions
for repeat_index in ${!repeats[@]}; do
	repeat=${repeats[repeat_index]};

	# Loop over all the bandwidths
	for bw_index in ${!bws[@]}; do
		bw=${bws[bw_index]};
		bw_param_string="bw=$bw";

		# Set cf1 parameter according to the bw
		cf1_param_string_2g="";
		cf1_param_string_5g="";
		if [[ $bw -eq 40 ]]; then
			cf1_param_string_2g="cf1=$cf1_2g_40m";
			cf1_param_string_5g="cf1=$cf1_5g_40m";
		elif [[ $bw -eq 80 ]]; then
			cf1_param_string_5g="cf1=$cf1_5g_80m";
		fi

		# Loop over spb values
		for spb_index in ${!spbs[@]}; do
			spb=${spbs[spb_index]};

			if [[ -f $config_file ]]; then
				rm $config_file;
			fi

			# Write the content of config file
			if [[ $bw -ne 80 ]]; then
				echo "$mac_address_2g bw=$bw cf=$cf_2g $cf1_param_string_2g spb=$spb" >> $config_file;
			fi
			echo "$mac_address_5g bw=$bw cf=$cf_5g $cf1_param_string_5g spb=$spb" >> $config_file;

			# Print the measurement parameters
			echo "#### repeat=$repeat bw=$bw spb=$spb"

			# Repeat the measurement
			for i in `seq $repeat`; do
		    	$command;
		    	sleep 5;
		    done			
		done
	done
done

echo "#### End of experiment ###"
