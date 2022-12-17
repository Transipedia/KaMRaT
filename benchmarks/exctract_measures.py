from sys import argv


def parse_file(file):
	values = {}

	with open(file) as fp:
		for line in fp:
			line = line.strip()

			# Wall time in s
			if "Elapsed (wall clock) time (h:mm:ss or m:ss)" in line:
				time = line[line.rfind(' ')+1:]
				time = [int(x) for x in time.split(':')]
				
				values['walltime'] = time[-1] + 60 * time[-2]
				if len(time) == 3:
					values['walltime'] += 3600 * time[0]

			# User time
			if "User time (seconds)" in line:
				values['usertime'] = float(line[line.rfind(' ')+1:])

			# Sys time
			if "System time (seconds)" in line:
				values['systime'] = float(line[line.rfind(' ')+1:])

			# Memory peak
			if "Maximum resident set size (kbytes)" in line:
				values['mempeak'] = int(line[line.rfind(' ')+1:]) // 1024

	return values


if __name__ == "__main__":
	values = parse_file(argv[1])
	values["waittime"] = max(0, values["walltime"] - values["usertime"] - values["systime"])
	values["waittime"] = f'{values["waittime"]:.2f}'

	filename = argv[1]
	filename = filename[filename.rfind('/')+1 : filename.rfind('.')]
	to_print = ["walltime", "usertime", "systime", "waittime", "mempeak"]
	to_print = "\t".join([str(values[x]) for x in to_print])
	print(f'{filename}\t{to_print}')
