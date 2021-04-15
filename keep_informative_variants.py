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
			var_id = '|'.join(entries[0:4])
			obs_vars = []
			uniq_vars = {}
			compare_var_present = False
			target_var_present = False
			for i in range((len(entries) - 4) // 3):
				this_idx = 4 + (i * 3)
				obs_vars.append(entries[this_idx])
				if i == target_idx and entries[this_idx] != '.':
					target_var_present = True
				elif i != target_idx and entries[this_idx] != '.':
					compare_var_present = True
			if not (compare_var_present and target_var_present):
				continue
			for var in obs_vars:
				uniq_vars.setdefault(var, len(uniq_vars) + 1)
			sys.stdout.write('{},{}\n'.format(
				var_id,
				','.join([str(uniq_vars[i]) for i in obs_vars])
			))


if __name__ == '__main__':
	load_vcf_tab(sys.argv[1], sys.argv[2])
