#!/usr/bin/env python3

import sys
import os
import fnmatch
import argparse
import datetime

if sys.version_info[0] < 3:
	raise Exception('Python 3 or a more recent version is required.')


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

	for fq_path in fq_pass:
		for _, _, filenames in os.walk(fq_path):
			for f in filenames:
				print(f)

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
	check_file_match(args.root_dir, args.dest_dir, fq_pass_paths)
