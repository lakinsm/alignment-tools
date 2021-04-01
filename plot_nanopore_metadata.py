#!/usr/bin/env python3

import sys
import numpy as np
import queue
from queue import PriorityQueue
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_cov(cov_dict, target, target_len, window_size, pdf_handle):
	# y = rolling_average([(k, v) for k, v in cov_dict[target].items()], window_size, np.ceil(window_size / 2),
	#                     target_len)
	x = np.ceil(window_size / 2) * np.array(range(len(y)))
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


def parse_seqfile(infile):
	ret = {}
	with open(infile, 'r') as f:
		f.readline()
		line = f.readline()
		while line:
			entries = line.split()
			fq_name = entries[0]
			if '_fail_' not in fq_name:
				sectime = float(entries[6])
				if fq_name not in ret:
					ret[fq_name] = sectime
				elif ret[fq_name] < sectime:
					ret[fq_name] = sectime
			line = f.readline()
	return ret


def parse_throughfile(infile):
	ret = tuple()
	with open(infile, 'r') as f:
		f.readline()
		line = f.readline()
		while line:
			entries = line.split(',')
			mintime = int(entries[0])
			bcr_pass = int(entries[2])
			bcr_fail = int(entries[3])
			bcbp = int(entries[8])
			ret += ((mintime, bcr_pass, bcr_fail, bcbp),)
			line = f.readline()
	return ret


if __name__ == '__main__':
	seqfile_data = parse_seqfile(sys.argv[1])
	through_data = parse_throughfile(sys.argv[2])
	print(seqfile_data)
	print('\n')
	print(through_data)
