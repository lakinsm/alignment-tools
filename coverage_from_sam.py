#!/usr/bin/env python3

import sys
import numpy as np
import queue
from queue import PriorityQueue
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class SamParser(object):
	"""
	Line-by-line parsing of Sequence Alignment Map (SAM) UTF-8 encoded files.  Stores alignment information for
	genome lengths from the provided reference database using the SAM headers, if specified.  Outputs
	aligned reads if specified and only if the aligned reads do not match a provided filter file (headers of reference
	to exclude).  Aligning reads, if output, are output as they are parsed.
	"""

	def __init__(self, sam_path, screen_headers=None, output_reads=None, capture_ref_len=True):
		"""
		Initialize object and data structures for storing coverage and coverage over time, if specified.
		:param sam_path: STR, path to input SAM file
		:param screen_headers: LIST, reference headers to ignore if encountered in SAM file alignments
		:param output_reads: STR, output aligned reads passing filter/screen to this filepath in FASTA format
		:param capture_ref_len: BOOL, store reference genome lengths from SAM headers
		"""
		self.sam_path = sam_path
		self.ref_len = {}
		if screen_headers:
			self.filter = set(screen_headers)
		self.output_reads = output_reads
		self.output_handle = None
		if output_reads:
			self.seen_reads = set()
		self.handle = None

		# Read SAM header information to retrieve database headers and sequence lengths
		self._open()
		self.line = self.handle.readline()
		if self.line:
			while self.line[0] == '@':
				if capture_ref_len and self.line.startswith('@SQ'):
					entries = self.line.split('\t')
					db_header = entries[1].split(':')[-1]
					db_seq_len = int(entries[2].split(':')[-1])
					self.ref_len[db_header] = db_seq_len
				self.line = self.handle.readline()

	def __iter__(self):
		return self

	def __next__(self):
		"""
		Iterator for SAM file handle, yields aligned read entries. The start_idx yielded from the SAM file is
		one-indexed and must have 1 subtracted from it to translate it into zero-indexing.
		:yield: (query, query_reverse_bool, target, target_start_idx, CIGAR)
		"""
		if not self.line:
			self._close()
			raise StopIteration
		entries = self.line.split('\t')
		sam_flag = int(entries[1])
		while sam_flag & 4 != 0:
			self.line = self.handle.readline()
			if not self.line:
				self._close()
				raise StopIteration
			entries = self.line.split('\t')
			sam_flag = int(entries[1])
		query_header = entries[0]
		query_seq = entries[9]
		target_header = entries[2]
		target_start = int(entries[3])
		cigar = entries[5]
		if sam_flag & 16 != 0:  # if reverse
			reverse = True
		else:
			reverse = False
		if self.output_reads:
			if query_header not in self.seen_reads:
				self.seen_reads.add(query_header)
				self.output_handle.write('>{}\n{}\n'.format(
					query_header,
					query_seq
				))
		self.line = self.handle.readline()
		return query_header, reverse, target_header, target_start, cigar

	def _open(self):
		self.handle = open(self.sam_path, 'r')
		if self.output_reads:
			self.output_handle = open(self.output_reads, 'w')

	def _close(self):
		self.handle.close()
		if self.output_reads:
			self.output_handle.close()


def parse_cigar(s, t_idx, rev):
	"""
	Parse SAM CIGAR alignment string and return indices to which the read aligned.
	:param s: STR, CIGAR string
	:param t_idx: INT, zero-based index for target start position
	:param rev: BOOL, read aligned in reverse orientation if true
	:return: tuple of integers, zero-indexed indices to which the read aligned
	"""
	ret = ()
	num = ''
	c_idx = 0
	while c_idx < len(s):
		if s[c_idx].isdigit():
			num += s[c_idx]
		else:
			op = s[c_idx]
			if rev:
				if op == 'M' or op == '=':
					ret += tuple(range(t_idx - int(num) + 1, t_idx + 1))
					t_idx -= int(num)
				elif op == 'D':
					t_idx -= int(num)
				elif op == 'N':
					t_idx -= int(num)
				elif op == 'X':
					ret += tuple(range(t_idx - int(num) + 1, t_idx + 1))
					t_idx -= int(num)
			else:
				if op == 'M' or op == '=':
					ret += tuple(range(t_idx, t_idx + int(num)))
					t_idx += int(num)
				elif op == 'D':
					t_idx += int(num)
				elif op == 'N':
					t_idx += int(num)
				elif op == 'X':
					ret += tuple(range(t_idx, t_idx + int(num)))
					t_idx += int(num)
			num = ''
		c_idx += 1
	return ret


def find_top_targets(cov_dict, n=5):
	"""
	Using a min priority queue, store top 5 coverage results and return a list of genome names to target.
	:param cov_dict: DICT, dictionary of coverage over indices from SamParser {target: {idx: cov}}
	:param n: INT, number of top coverage genomes to return
	:return: tuple of genome names (target keys) in order of coverage, descending
	"""
	ret = ()
	tops = PriorityQueue()
	for k, v in cov_dict.items():
		tops.put((-sum(v.values()), k))
	for _ in range(n):
		try:
			ret += (tops.get(block=False)[1],)
		except queue.Empty:
			break
	return ret


def rolling_average(data, window_size, interval, max_len):
	"""
	Generate a vector of averages for each window of size window_size at fixed intervals interval.
	:param data: TUPLE, (x, y) pairs of integers, where x is the genome idx and y is the depth over that idx
	:param window_size: INT, size of each window over which to compute the averages
	:param interval: INT, fixed interval to shift the window
	:param max_len: INT, maximum index for the end window
	:return: tuple of rolling averages
	"""
	ret = ()
	local_data = np.array(data)
	local_data.sort(axis=0)
	if window_size > max_len:
		window_size = np.ceil(max_len / 50)
		interval = np.ceil(window_size / 2)
	window_start = 0
	window_end = window_size
	while window_end <= max_len:
		idxs = np.where((local_data[:, 0] < window_end) * (local_data[:, 0] >= window_start))
		ret += (np.sum(local_data[idxs, 1]) / float(window_size),)
		window_start += interval
		window_end += interval
	return ret


def plot_cov(cov_dict, target, target_len, window_size, pdf_handle):
	"""
	Using matplotlib.pyplot, produce coverage plots for the top n targets and output to PDF, one per page.
	:param cov_dict: DICT, dictionary of coverage over indices from SamParser {target: {idx: cov}}
	:param target: STR, reference genome name to plot, must be present in cov_dict
	:param target_len: INT, reference genome length of target
	:param window_size: INT, window size over which to produce rolling averages
	:param pdf_handle: FILE HANDLE, pdf file handle object for saving figures
	:return: None
	"""
	y = rolling_average([(k, v) for k, v in cov_dict[target].items()], window_size, np.ceil(window_size / 2),
	                    target_len)
	x = np.array(range(int(window_size) * len(y)))
	plt.figure(figsize=(15, 10))
	plt.fill_between(x, y, 0,
	                 facecolor='blue',
	                 color='black',
	                 alpha=0.3)
	plt.xlabel('Genomic Location')
	plt.ylabel('Coverage (Rolling Average)')
	plt.title('Rolling Average Coverage, Reference Sequence {}'.format(target))
	plt.ticklabel_format(style='sci', axis='x')
	pdf_handle.savefig()
	plt.close()


def write_coverage(cov_dict, ref_len_dict, output_csv_path):
	"""
	Write a long-format .csv file containing genomic coverage percentages at each timepoint for all targets.
	:param cov_dict: DICT, dictionary of coverage over indices from worker {barcode: {target: set(idxs)}
	:param ref_len_dict: DICT, dictionary of reference genome lengths {target: length}
	:param output_csv_path: STR, filepath of the output .csv file
	:return: None
	"""
	target_info = tuple()
	for target, idx_dict in cov_dict.items():
		genome_len = ref_len_dict[target]
		bases_aligned = sum(idx_dict.values())
		bases_covered = len(set(range(genome_len)).intersection(idx_dict.keys()))
		target_info += ((bases_aligned, bases_covered, genome_len, target),)
	sorted_targets = sorted(target_info, key=lambda x: x[0], reverse=True)
	for bases_aligned, bases_covered, genome_len, target in sorted_targets:
		sys.stdout.write('\nTarget: {}\nGenome Length: {}\nAverage Coverage: {}'
		                 '\nPercent Bases Covered: {}\n\n'.format(
			target,
			genome_len,
			float(bases_aligned) / float(genome_len),
			100 * float(bases_covered) / float(genome_len)
		))
	with open(output_csv_path, 'w') as out:
		for target, idx_dict in cov_dict.items():
			out.write('{}'.format(target))
			for idx in range(ref_len_dict[target]):
				try:
					out.write(',{}'.format(idx_dict[idx]))
				except KeyError:
					out.write(',0')
			out.write('\n')


def worker(infile):
	ref_cov = {}
	sam_parser = SamParser(infile)
	for query, q_reverse, target, t_start, cigar in sam_parser:
		idxs = parse_cigar(cigar, t_start - 1, q_reverse)
		if target not in ref_cov:
			ref_cov[target] = {}
		for i in idxs:
			if i in ref_cov[target]:
				ref_cov[target][i] += 1
			else:
				ref_cov[target][i] = 1
	return ref_cov


if __name__ == '__main__':
	# 1. SAM input file, 2. output CSV file, 3. output PDF file
	this_sam_parser = SamParser(sys.argv[1])
	res = worker(sys.argv[1])
	write_coverage(res, this_sam_parser.ref_len, sys.argv[2])

	top_targets = find_top_targets(res)
	n_windows = 500  # Number of points to graph rolling avg for coverage plots
	with PdfPages(sys.argv[3]) as pdf:
		for t in top_targets:
			# k = 2n / (m + 1)
			window_size = np.ceil(float(2 * this_sam_parser.ref_len[t]) / float(n_windows + 1))
			window_size = np.max((window_size, 1))
			plot_cov(res, t, this_sam_parser.ref_len[t], window_size, pdf)
