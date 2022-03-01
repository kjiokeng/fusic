import re
import datetime
import operator
import sys
from os import listdir
from os.path import isfile, join, isdir
import argparse
import numpy as np
import matplotlib.pyplot as plt

# File to be treated
ftm_measures_filename = "ftm-report-exp-add-obs5.txt"

# Column separator (CSV-like format)
sep = ";"

def format_results(measures_filename, dir_path="", print_details=True,
	plot=False, remove_outliers=False, begin=10, end=10):
	# Result dictionary
	result = {}
	dist_mean = 0
	dist_std = 0
	rssi_mean = 0
	rssi_std = 0

	# Parse the ftm measures
	# Regex to look for time and associated computed distance
	regex_pattern = re.compile("\[(.*)\] ([a-zA-Z: 0-9,]*) distance: ([-0-9]+) cm, rssi: ([-0-9]+) dBm\n");
	with open(join(dir_path, measures_filename), 'r') as file:
		for line in file:
			if not line.endswith("dBm\n"):
				continue

			regex_search_result = regex_pattern.search(line)
			if regex_search_result:
				time = regex_search_result.group(1)
				time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
				distance = regex_search_result.group(3)
				distance = float(distance) / 100 # Conversion from string to float, then from centimeters to meters
				rssi = regex_search_result.group(4)
				rssi = float(rssi)
				if distance != 0: # Remove unsuccessful measurements
					result[time] = {"distance": distance, "rssi": rssi}


	

	# Sort the result by time
	base_time = None # Time to count as t=0:00
	result = sorted(result.items(), key=operator.itemgetter(0))
	epoch_datetime=datetime.datetime(1970,1,1)
	distance_mean = 0
	rssi_mean = 0;
	i = 0
	n = 0
	total = len(result)

	timestamps = []
	distances = []
	rssis = []

	for key, value in result:
		i += 1

		if i<=begin or i>=total-end:
			continue

		# Set the "base time" if not yet set
		if not base_time:
			base_time = key

		distance = value.get("distance")
		rssi = value.get("rssi")

		# Outlier detection : ignore values which are more than 2 times higher than the mean
		if (distance * n > distance_mean * 2) and remove_outliers:
			continue

		relative_time = (epoch_datetime+(key-base_time)).strftime("%M:%S")
		# print str(i) + " " + str(distance)
		if print_details:
			print key.strftime("%H:%M:%S") + sep + relative_time \
			  + sep + str(rssi) + sep + str(distance)
			
		timestamps.append(relative_time)
		distances.append(distance)
		rssis.append(rssi)

		distance_mean += distance
		rssi_mean += rssi
		n += 1

	if n!=0:
		dist_array = np.array(distances)
		dist_mean = np.mean(dist_array)
		dist_std = np.std(dist_array)

		rssi_array = np.array(rssis)
		rssi_mean = np.mean(rssi_array)
		rssi_std = np.std(rssi_array)

		print measures_filename + "\n\t" \
			 + "Dist. Mean: " + str(dist_mean) + "\tStd: " + str(dist_std) \
			 + "\n\trssi. Mean: " + str(rssi_mean) + "\tStd: " + str(rssi_std) + "\n"


		if plot:
			plt.plot(distances)
			# plt.xticks(range(len(timestamps)), timestamps, size='small')
			plt.xlabel('Time (s)')
			plt.ylabel('Computed distance')
			bottom, top = plt.ylim()
			plt.ylim(0, top*1.25)
			plt.show()

	return {'dist_mean': dist_mean, 'dist_std': dist_std,
			'rssi_mean': rssi_mean, 'rssi_std': rssi_std}




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Parses and analyses ftm measurement output.')
	parser.add_argument('-p', '--plot', default=False, action="store_true",
		help='Plot the graph of computed distance vs time')
	parser.add_argument('-d', '--details', default=False, action="store_true",
		help='Print execution details')
	parser.add_argument('-o', '--outliers', default=False, action="store_true",
		help='Remove outlier values')
	parser.add_argument('-b', '--begin', default=5, type=int,
		help='Will ignore the BEGIN first measurements.')
	parser.add_argument('-e', '--end', default=5, type=int,
		help='Will ignore the END last measurements.')
	parser.add_argument('filename',
		help='The file to process. If filename is a directory, then all the files contained in that directory are going to be processed')
	args = parser.parse_args()

	if isdir(args.filename):
		dir_path = args.filename
		files = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]
		result = []
		for f in sorted(files):
			res = format_results(measures_filename=f, dir_path=dir_path, begin=args.begin, end=args.end,
				print_details=args.details, plot=args.plot, remove_outliers=args.outliers)
			sys.stdout.flush()

			if not f.startswith("summary"):
				res['file'] = f;
				# res['number'] = int(f.replace("pos", "").replace(".txt", ""))
				result.append(res)

		print "\n\n########## Global summary (csv) ##########"
		# result = sorted(result, key=lambda k: k['number'])
		print 'file,dist_mean,dist_std,rssi_mean,rssi_std'
		for res in result:
			print res['file'] + "," + str(res['dist_mean']) + "," + str(res['dist_std']) \
			+ "," + str(res['rssi_mean']) + "," + str(res['rssi_std'])

	else:
		format_results(measures_filename=args.filename, begin=args.begin, end=args.end,
		print_details=args.details, plot=args.plot, remove_outliers=args.outliers)

