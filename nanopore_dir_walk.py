#!/usr/bin/env python3

import sys
import os
import fnmatch
import argparse
import datetime
import re

if sys.version_info[0] < 3:
	raise Exception('Python 3 or a more recent version is required.')


MINKNOW_REGEX = re.compile(r'[0-9]{8}_[0-9]{4}_[0-9A-Z\-]+_([A-Z]{3}[0-9]+)_')


def find_fastq_pass(root_dir):
	"""
	Take as input an absolute path root directory for Nanopore instrument data and output all paths to fastq_pass
	folders within that root directory.
	:param root_dir: STR, absolute directory path to instrument-level root directory
	:return: tuple of absolute paths to all fastq_pass directories in root_dir
	"""
	fq_pass = tuple()
	for root, dirnames, _ in os.walk(root_dir):
		for d in fnmatch.filter(dirnames, 'fastq_pass'):
			fq_pass += (os.path.join(root, d),)
	return fq_pass


def find_ont_sample_flowcell(fpath):
	"""
	Return the samplename and flowcell name as specified by the input directory path.  This uses regular expressions
	to identify the flowcell and finds the samplename from the parent directories to the flowcell ID.
	:param fpath: STR, directory path for a fastq_pass folder
	:return: TUPLE of STR, (samplename, flowcell_id)
	"""
	flowcell_id = MINKNOW_REGEX.search(fpath).group(1)
	sample1, sample2 = fpath.split(flowcell_id)[0].split('/')[-3:-1]
	if sample2 == 'no_sample':
		samplename = sample1
	elif sample1 == sample2:
		samplename = sample2
	else:
		samplename = sample1
	return samplename, flowcell_id


def check_file_match(root_source, root_dest, fq_pass, write_text_log=None):
	"""
	Compare FASTQ and Nanopore flowcell metadata files from source to destination, and update destination files if
	not matching from source.  This is a destructive/overwrite operation that uses a temporary intermediate file.
	This operation will never delete existing data (if it is absent from source but present in destination).  A logfile
	is optionally written by default with date/time of file updates or status checks.
	:param root_source: STR, absolute directory path to instrument-level root directory for source files
	:param root_dest: STR, absolute directory path to instrument-level root directory for destination files
	:param fq_pass: TUPLE, an array of absolute paths to all fastq_pass directories in the source root directory
	:param write_text_log: NONE/STR, write a tab-delimited log of date/time, source filepath, dest filepath, and status
	:return: None
	"""
	log_handle = None
	if write_text_log:
		log_handle = open(write_text_log, 'a')
		log_handle.write('\n\n{}\tBegin file checking and copy, source: {}, destination: {}\n'.format(
			datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
			root_source,
			root_dest
		))

	n_existing = 0
	n_updated = 0
	for fq_path in fq_pass:
		samplename, flowcell_id = find_ont_sample_flowcell(fq_path)
		stdout_prefix = '{} Checking fastq_pass'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
		sys.stdout.write('{} (updated={}, existing={}), sample: {}, flowcell: {}\r'.format(
			stdout_prefix,
			n_updated,
			n_existing,
			samplename,
			flowcell_id
		))
		for this_root, _, filenames in os.walk(fq_path):
			for f in filenames:
				source_path = os.path.join(this_root, f)
				dest_path = '{}/{}'.format(root_dest.rstrip('/'), source_path.split(root_source)[-1])
				# Check if exists
				dest_isfile = os.path.isfile(dest_path)
				if not dest_isfile:
					n_updated += 1
				if (n_existing + n_updated) % 5 == 0:
					sys.stdout.write('{} (updated={}, existing={}), sample: {}, flowcell: {}\r'.format(
						stdout_prefix,
						n_updated,
						n_existing,
						samplename,
						flowcell_id
					))
		sys.stdout.write('{} (updated={}, existing={}), sample: {}, flowcell: {}\n'.format(
			stdout_prefix,
			n_updated,
			n_existing,
			samplename,
			flowcell_id
		))
	if log_handle:
		log_handle.close()


parser = argparse.ArgumentParser('nanopore_dir_walk.py')
parser.add_argument('root_dir', default=None, type=str,
                    help='Root instrument directory to check for updated files')
parser.add_argument('dest_dir', default=None, type=str,
                    help='Destination directory to check for existing files to update')


if __name__ == '__main__':
	args = parser.parse_args()
	fq_pass_paths = find_fastq_pass(os.path.realpath(args.root_dir))
	check_file_match(os.path.realpath(args.root_dir), os.path.realpath(args.dest_dir), fq_pass_paths)
