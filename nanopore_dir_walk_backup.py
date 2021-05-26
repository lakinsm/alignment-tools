#!/usr/bin/env python3

import sys
import os
import fnmatch
import argparse
import datetime
import re
import hashlib
import shutil
import glob
import time

if sys.version_info[0] < 3:
	raise Exception('Python 3 or a more recent version is required.')


MINKNOW_REGEX = re.compile(r'[0-9]{8}_[0-9]{4}_[0-9A-Z\-]+_([A-Z]{3}[0-9]+)_')
OS_WALK_EXCLUDES = {'fast5', 'fast5_pass', 'fast5_fail', 'fastq_fail'}


def os_walk_condition(r, s):
	"""
	os.walk() is a breadth-first search of directories and files.  This function returns a boolean of whether we want
	to explore subdirectories or not, since some of these subdirectories (OS_WALK_EXCLUDES) contain many files that
	take a lot of time to walk.  All branches in the file tree of length greater than 2 whose basedir is in
	OS_WALK_EXCLUDES are truncated and excluded from further search during the BFS walk.
	:param r: STR, the root directory path for input/source files
	:param s: STR, path to the directory currently under consideration for next walk (in this breadth level)
	:return: BOOL, True if include in walk, False otherwise
	"""
	status = False
	root_base = r.split('/')[-1]
	branch = root_base + '/' + s.split(r)[-1].lstrip('/')
	branch_split = branch.split('/')
	if len(branch_split) > 2:
		if branch_split[-1] not in OS_WALK_EXCLUDES:
			status = True
	else:
		status = True
	return status


def find_fastq_pass(root_dir):
	"""
	Take as input an absolute path root directory for Nanopore instrument data and output all paths to fastq_pass
	folders within that root directory.
	:param root_dir: STR, absolute directory path to instrument-level root directory
	:return: tuple of absolute paths to all fastq_pass directories in root_dir
	"""
	fq_pass = tuple()
	for root, dirs, _ in os.walk(root_dir + '/', topdown=True):
		dirs[:] = [d for d in dirs if os_walk_condition(root_dir, os.path.join(root, d))]
		for d in fnmatch.filter(dirs, 'fastq_pass'):
			fq_pass += (os.path.join(root, d),)
	return fq_pass


def find_ont_metadata(flowcell_dir):
	"""
	Take as input an absolute path flowcell directory for Nanopore instrument data and output all paths to ONT metadata
	files within that root directory.
	:param flowcell_dir: STR, absolute directory path to flowcell-level root directory
	:return: tuple of absolute paths to all fastq_pass directories in root_dir
	"""
	ont_files = tuple()
	for root, dirs, files in os.walk(flowcell_dir + '/', topdown=True):
		dirs[:] = [d for d in dirs if os_walk_condition(flowcell_dir, os.path.join(root, d))]
		for f in files:
			if fnmatch.fnmatch(f, '*.csv') or fnmatch.fnmatch(f, '*.tsv') or fnmatch.fnmatch(f, '*.txt') or \
				fnmatch.fnmatch(f, '*.md') or fnmatch.fnmatch(f, '*.pdf'):
				ont_files += (os.path.join(root, f),)
	return ont_files


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


def print_file_status(prefix, updated, existing, sname, fname, terminal=False):
	"""
	Print file parsing status information to stdout to update the user.
	:param prefix: STR, prefix containing the start date/time
	:param updated: INT, number of updated files
	:param existing: INT, number of existing files (not updated)
	:param sname: STR, samplename of current fastq_pass directory
	:param fname: STR, flowcell ID of current fastq_pass directory
	:param terminal: BOOL, add line ending
	:return: None
	"""
	sys.stdout.write('{}\tupdated: {}\texisting: {}\tflowcell: {}\tsample: {}'.format(
		prefix,
		updated,
		existing,
		fname,
		sname
	))
	if terminal:
		sys.stdout.write('\n')
	else:
		sys.stdout.write('\r')


def file_as_bytes(fhandle):
	"""
	Read file and return as bytes.
	:param fhandle: FILE OBJECT, open file handle from python's open() in read and binary mode
	:return: BYTES OBJECT, binary file as output from python's read()
	"""
	with fhandle:
		return fhandle.read()


def update_dest_file_robust(spath, dpath, source_sha256, n_tries=3, sec_wait=30):
	"""
	Copy source to destination temporary file, then move into real destination file once copy is complete and sha256
	sum has been verified.
	:param spath: STR, source file path
	:param dpath: STR, destination file path
	:param source_sha256: STR, source SHA256 file hash sum
	:param n_tries: INT, number of times to try copy to temp file with subsequent hash
	:param sec_wait: INT, number of seconds to wait between tries
	:return: None
	"""
	parent_dir = '/'.join(dpath.split('/')[:-1])
	parent_dir_status = os.path.isdir(parent_dir)
	if not parent_dir_status:
		os.makedirs(parent_dir, exist_ok=True)

	this_tries = 0
	while this_tries < n_tries:
		shutil.copy2(spath, dpath + '_temp')
		temp_dest_hash = hashlib.sha256(file_as_bytes(open(dpath + '_temp', 'rb'))).hexdigest()
		if temp_dest_hash == source_sha256:
			shutil.move(dpath + '_temp', dpath)
			break
		else:
			sys.stderr.write('\nCopy hashsum failed for source file: {}, destination file: {}\n'.format(spath, dpath))
			dest_temp_isfile = os.path.isfile(dpath + '_temp')
			if dest_temp_isfile:
				os.remove(dpath + '_temp')
			this_tries += 1
			time.sleep(sec_wait)
	if this_tries == n_tries:
		sys.stderr.write('\nRETRIES COUNT EXCEEDED, Copy hashsum failed for source file: {}, destination file: {}\n'.format(spath, dpath))
		raise ValueError


def check_file_match(root_source, root_dest, fq_pass, write_text_log=None, non_interactive=False):
	"""
	Compare FASTQ and Nanopore flowcell metadata files from source to destination, and update destination files if
	not matching from source.  This is a destructive/overwrite operation that uses a temporary intermediate file.
	This operation will never delete existing data (if it is absent from source but present in destination).  A logfile
	is optionally written by default with date/time of file updates or status checks.
	:param root_source: STR, absolute directory path to instrument-level root directory for source files
	:param root_dest: STR, absolute directory path to instrument-level root directory for destination files
	:param fq_pass: TUPLE, an array of absolute paths to all fastq_pass directories in the source root directory
	:param write_text_log: NONE/STR, write a tab-delimited log of date/time, source filepath, dest filepath, and status
	:param slurm: BOOL, don't use carriage returns so that SLURM logs are properly formatted.  Use with SLURM scheduler.
	:return: None
	"""
	log_handle = None
	if write_text_log:
		log_handle = open(write_text_log, 'w')
		log_handle.write('{}\tBegin file checking and copy, source: {}, destination: {}\n'.format(
			datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
			root_source,
			root_dest
		))

	n_existing = 0
	n_updated = 0
	for fq_path in fq_pass:
		samplename, flowcell_id = find_ont_sample_flowcell(fq_path)
		stdout_prefix = '{} Checking fastq_pass'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
		print_file_status(stdout_prefix, n_updated, n_existing, samplename, flowcell_id, non_interactive)

		# Metadata files
		flowcell_dir = '/'.join(fq_path.rstrip('/').split('/')[:-1])
		ont_metadata_files = find_ont_metadata(flowcell_dir)
		print(ont_metadata_files)
		sys.exit()
		for f in ont_metadata_files:
			fstatus = 'UPDATED'
			dest_path = '{}/{}'.format(root_dest.rstrip('/'), f.split(root_source)[-1].lstrip('/'))

			# Check if exists
			dest_isfile = os.path.isfile(dest_path)
			dest_temp_isfile = os.path.isfile(dest_path + '_temp')
			if dest_temp_isfile:
				os.remove(dest_path + '_temp')
			source_hash = hashlib.sha256(file_as_bytes(open(f, 'rb'))).hexdigest()
			if not dest_isfile:
				update_dest_file_robust(f, dest_path, source_hash)
				n_updated += 1
			else:
				# Check SHA256 hash sum
				dest_hash = hashlib.sha256(file_as_bytes(open(dest_path, 'rb'))).hexdigest()
				if source_hash == dest_hash:
					fstatus = 'EXISTING'
					n_existing += 1
				else:
					update_dest_file_robust(f, dest_path, source_hash)
					n_updated += 1

			# Logging
			if (n_existing + n_updated) % 5 == 0:
				print_file_status(stdout_prefix, n_updated, n_existing, samplename, flowcell_id, non_interactive)
			if log_handle:
				# datetime, sample, flowcell, status, source, dest
				log_handle.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
					datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
					samplename,
					flowcell_id,
					fstatus,
					f,
					dest_path
				))

		# FASTQ files
		for this_root, _, filenames in os.walk(fq_path):
			for f in filenames:
				fstatus = 'UPDATED'
				source_path = os.path.join(this_root, f)
				dest_path = '{}/{}'.format(root_dest.rstrip('/'), source_path.split(root_source)[-1].lstrip('/'))

				# Check if exists
				dest_isfile = os.path.isfile(dest_path)
				dest_temp_isfile = os.path.isfile(dest_path + '_temp')
				if dest_temp_isfile:
					os.remove(dest_path + '_temp')
				source_hash = hashlib.sha256(file_as_bytes(open(source_path, 'rb'))).hexdigest()
				if not dest_isfile:
					update_dest_file_robust(source_path, dest_path, source_hash)
					n_updated += 1
				else:
					# Check SHA256 hash sum
					dest_hash = hashlib.sha256(file_as_bytes(open(dest_path, 'rb'))).hexdigest()
					if source_hash == dest_hash:
						fstatus = 'EXISTING'
						n_existing += 1
					else:
						update_dest_file_robust(source_path, dest_path, source_hash)
						n_updated += 1

				# Logging
				if (n_existing + n_updated) % 5 == 0:
					print_file_status(stdout_prefix, n_updated, n_existing, samplename, flowcell_id, non_interactive)
				if log_handle:
					# datetime, sample, flowcell, status, source, dest
					log_handle.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
						datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
						samplename,
						flowcell_id,
						fstatus,
						source_path,
						dest_path
					))
		print_file_status(stdout_prefix, n_updated, n_existing, samplename, flowcell_id, terminal=True)
	if log_handle:
		log_handle.close()


parser = argparse.ArgumentParser('nanopore_dir_walk_backup.py')
parser.add_argument('source_dir', default=None, type=str,
                    help='Source instrument root directory to check for updated files')
parser.add_argument('dest_dir', default=None, type=str,
                    help='Destination backup directory to check for existing files to update')
parser.add_argument('-l', '--logfile', default=None, type=str,
                    help='Path to output log file')
parser.add_argument('-n', '--non_interactive', default=False, action='store_true',
                    help='Flag indicated no carriage returns so that logs are formatted properly')


if __name__ == '__main__':
	args = parser.parse_args()
	fq_pass_paths = find_fastq_pass(os.path.realpath(args.source_dir))
	check_file_match(os.path.realpath(args.source_dir),
	                 os.path.realpath(args.dest_dir),
	                 fq_pass_paths,
	                 write_text_log=args.logfile,
	                 non_interactive=args.non_interactive)
