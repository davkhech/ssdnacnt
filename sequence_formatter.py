"""
Utility to make dna structure readable for amber ff.

davkhech
"""
import argparse
import re


def parse_args(*argument_array):
	parser = argparse.ArgumentParser()
	parser.add_argument('input')
	parser.add_argument('--output', default='output.pdb')
	return parser.parse_args(*argument_array)


mappings = {
	'O1P': 'OP1',
	'O2P': 'OP2',
	'C1*': 'C1\'',
	'C2*': 'C2\'',
	'C3*': 'C3\'',
	'C4*': 'C4\'',
	'C5*': 'C5\'',
	'O3*': 'O3\'',
	'O4*': 'O4\'',
	'O5*': 'O5\'',
}


def main(args):
	input_file_name = args.input
	output_file_name = args.output
	acid_name = re.compile('  [ATGC] A ')
	atoms = mappings.keys()
	with open(input_file_name) as input_file, open(output_file_name, 'w') as output_file:
		for line in input_file:
			out_line = line
			match = acid_name.search(out_line)
			if match:
				coord = match.start() + 1
				out_line = out_line[:coord] + 'D' + out_line[coord + 1:]
			for atom in atoms:
				if atom in out_line:
					indx = out_line.index(atom)
					out_line = out_line[:indx] + mappings[atom] + out_line[indx + 3:]
			output_file.write(out_line)


if __name__ == '__main__':
	main(parse_args())

