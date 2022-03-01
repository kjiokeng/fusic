import re
import datetime
import operator

# Files to be treated
ftm_measures_filename = "ftm-report-exp-add-obs5.txt"
signal_power_filename = "signal-power-exp-add-obs5.txt"
obstacle_positions_filename = "obstacle-positions.txt"


def format_results(measures_filename, power_filename, positions_filename):
	# Result dictionary
	result = {}

	# Parse the ftm measures
	# Regex to look for time and associated computed distance
	regex_pattern = re.compile("\[(.*)\] ([a-zA-Z: 0-9,]*) distance: ([-0-9]+) cm\n");
	base_time = None # Time to count as t=0:00
	with open(measures_filename, 'r') as file:
		for line in file:
			if not line.endswith("cm\n"):
				continue

			regex_search_result = regex_pattern.search(line)
			if regex_search_result:
				time = regex_search_result.group(1)
				time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
				# Set the "base time" if not yet set
				if not base_time:
					base_time = time
				distance = regex_search_result.group(3)
				distance = float(distance) / 100 # Conversion from string to float, then from centimeters to meters
				if distance != 0: # Remove unsuccessful measurements
					result[time] = {"distance": distance}


	# Parse signal power measurements
	regex_pattern = re.compile("\[(.*)\] \tsignal: ([-0-9]+) dBm\n");
	with open(power_filename, 'r') as file:
		for line in file:
			regex_search_result = regex_pattern.search(line)
			if regex_search_result:
				time = regex_search_result.group(1)
				time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
				power = regex_search_result.group(2)
				power = float(power) # Conversion from string to float
				result_item = result.get(time)
				if result_item:
					result_item["power"] = power

	# Parse obstacle positions
	regex_pattern = re.compile("([-0-9]+):([-0-9]+)\t([-0-9]+)\n");
	with open(positions_filename, 'r') as file:
		for line in file:
			regex_search_result = regex_pattern.search(line)
			if regex_search_result:
				minutes = regex_search_result.group(1)
				seconds = regex_search_result.group(2)
				time_delta = datetime.timedelta(minutes=float(minutes), seconds=float(seconds))
				position = regex_search_result.group(3)
				position = float(position) # Conversion from string to float

				time = base_time + time_delta
				result_item = result.get(time)
				if result_item:
					result_item["position"] = position

	# Sort the result by time
	result = sorted(result.items(), key=operator.itemgetter(0))
	epoch_datetime=datetime.datetime(1970,1,1);
	sep = ";"
	for key, value in result:
		distance = value.get("distance")
		power = value.get("power")
		position = value.get("position")
		print key.strftime("%H:%M:%S") + sep + (epoch_datetime+(key-base_time)).strftime("%M:%S") \
		  + sep + str(power) + sep + str(position) + sep + str(distance)


if __name__ == '__main__':
	format_results(measures_filename=ftm_measures_filename, power_filename=signal_power_filename, positions_filename=obstacle_positions_filename)
