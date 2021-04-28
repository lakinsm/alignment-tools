#!/usr/bin/env python3

import sys
import os
import glob


def load_summary_file(infile_path):
	seen_timestamps = set()
	keep_files = set()
	output_data = tuple()
	with open(infile_path, 'r') as f:
		data = f.read().split('\n')
		for line in data[1:]:
			if not line:
				continue
			timestamp, status, samplename = line.split()
			seen_timestamps.add(timestamp)
			output_data += ((timestamp, status, samplename),)
			if status == 'UPDATED':
				keep_files.add(timestamp)
			elif status != 'NO_CHANGE':
				sys.stderr.write('Error: invalid status encountered: {}'.format(status))
				raise ValueError
	return seen_timestamps, keep_files, output_data


def check_delete(infile_path):
	this_delete_status = True
	with open(infile_path, 'r') as f:
		data = f.read().split('\n')
		for line in data[1:]:
			if not line:
				continue
			_, _, _, _, status, _, _ = line.split()
			if status == 'UPDATED':
				this_delete_status = False
	if this_delete_status:
		print('DELETE\t\t{}'.format(infile_path))


def clean_backup_logs(log_dir, summary_file_path):
	summary_file_present = False
	if os.path.exists(summary_file_path):
		summary_file_present = True
		seen, keep, summary_data = load_summary_file(summary_file_path)
	with open(summary_file_path, 'w') as sout:
		sout.write('Timestamp\tStatus\tSamplename\n')
		if summary_file_present:
			for timestamp, status, samplename in summary_data:
				sout.write('{}\t{}\t{}\n'.format(timestamp, status, samplename))
		for log_file in glob.glob(log_dir + '/*.log'):
			if os.stat(log_file).st_size == 0:
				continue
			timestamp = log_file.split('/')[-1].split('_')[0]
			if summary_file_present:
				if timestamp in keep:
					continue
				if timestamp in seen:
					check_delete(log_file)
			this_delete_status = True
			sample_status = {}
			with open(log_file, 'r') as f:
				data = f.read().split('\n')
				for line in data[1:]:
					if not line:
						continue
					_, _, samplename, _, status, _, _ = line.split()
					if samplename not in sample_status:
						sample_status[samplename] = False
					if status == 'UPDATED':
						this_delete_status = False
						sample_status[samplename] = True
			if this_delete_status:
				print('DELETE\t{}'.format(log_file))
			for samplename, status in sample_status.items():
				output_code = 'NO_CHANGE'
				if status:
					output_code = 'UPDATED'
				sout.write('{}\t{}\t{}\n'.format(
					timestamp,
					output_code,
					samplename
				))


if __name__ == '__main__':
	clean_backup_logs(sys.argv[1], sys.argv[2])
