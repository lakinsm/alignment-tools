#!/usr/bin/env python3

import sys
import os
import fnmatch
import argparse

if sys.version_info[0] < 3:
	raise Exception('Python 3.5 or a more recent version is required.')


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


parser = argparse.ArgumentParser('nanopore_dir_walk.py')
parser.add_argument('root_dir', default=None, required=True, type=str,
                    help='Root instrument directory to check for updated files')
parser.add_argument('dest_dir', default=None, required=True, type=str,
                    help='Destination directory to check for existing files to update')


if __name__ == '__main__':
	args = parser.parse_args()
	fq_pass_paths = find_fastq_pass(os.path.realpath(args.root_dir))
	print(fq_pass_paths)
