#!/usr/bin/env python3

import sys


def load_vcf_tab(infile, target):
	with open(infile, 'r') as f:
		data = f.read().split('\n')
		samplenames = []
		headers = data[0].split()[4:]
		for i in range(0, len(headers), 3):
			samplenames.append(headers[i].split('.')[0])
		sys.stdout.write('VariantID,{}\n'.format(
			','.join(samplenames)
		))
		target_idx = samplenames.index(target)
		for line in data[1:]:
			entries = line.split()
			var_id = '|'.join(entries[0:4]).replace(',', '_')
			obs_var = entries[3].split(',')[0]
			compare_var_present = False
			target_var_present = False
			obs_gts = []
			for i in range((len(entries) - 4) // 3):
				this_idx = 4 + (i * 3)
				if i == target_idx and entries[this_idx] != '.':
					target_var_present = True
				elif i != target_idx and entries[this_idx] != '.':
					compare_var_present = True
				obs_gts.append(entries[this_idx].count(obs_var))
			if not (compare_var_present and target_var_present):
				continue
			sys.stdout.write('{},{}\n'.format(
				var_id,
				','.join([str(x+1) for x in obs_gts])
			))


if __name__ == '__main__':
	load_vcf_tab(sys.argv[1], sys.argv[2])
