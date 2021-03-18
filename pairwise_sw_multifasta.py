#!/usr/bin/env python3

from skbio.core.ssw import align_striped_smith_waterman
import sys


def fasta_parse(infile):
	with open(infile, 'r') as fasta_file:
		# Skip whitespace
		while True:
			line = fasta_file.readline()
			if line is "":
				return  # Empty file or premature end of file?
			if line[0] is ">":
				break
		while True:
			if line[0] is not ">":
				raise ValueError("Records in FASTA should begin with '>'")
			header = line[1:].rstrip()
			all_lines = []
			line = fasta_file.readline()
			while True:
				if not line:
					break
				if line[0] is ">":
					break
				all_lines.append(line.rstrip())
				line = fasta_file.readline()
			yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
			if not line:
				return  # Stop Iteration


def pairwise_ssw(fasta, outfile):
	headers = tuple(sorted(fasta.keys()))
	with open(outfile, 'w') as out:
		out.write('query_header,target_header,aln_score_opt,aln_score_subopt,qstart,qend,tstart,tend,tend_subopt,cigar\n')
		for i, query_header in enumerate(headers):
			for j in range(i, len(headers)):
				target_header = headers[j]
				res = align_striped_smith_waterman(fasta[query_header], fasta[target_header])
				out.write('{},{},{},{},{},{},{},{},{},{}\n'.format(
					query_header,
					target_header,
					res['optimal_alignment_score'],
					res['suboptimal_alignment_score'],
					res['query_begin'],
					res['query_end'],
					res['target_begin'],
					res['target_end_optimal'],
					res['target_end_suboptimal'],
					res['cigar']
				))


if __name__ == '__main__':
	D = {k: v for k, v in fasta_parse(sys.argv[1])}
	pairwise_ssw(D, sys.argv[2])
