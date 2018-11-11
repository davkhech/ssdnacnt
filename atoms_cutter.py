"""
Utility to remove atoms with given indices from gromacs topology file.

davkhech
"""

import argparse

NUMBER_OF_INDICES = {
        'atoms': 1,
        'bonds': 2,
        'pairs': 2,
        'angles': 3,
        'dihedrals': 4,
}


def parse_args(*argument_array):
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('range', nargs='+')
    parser.add_argument('--output', default='output.itp')
    return parser.parse_args(*argument_array)


def remove_from_list(remove_atoms):
    def _remover(check_list):
        for atom in remove_atoms:
            if atom in check_list:
                return True
        return False
    return _remover


def renumber_one(atoms):
    def _renumber(atom):
        x = 0
        for a in atoms:
            if atom > a:
                x += 1
            else:
                return atom - x
        return atom - x
    return _renumber


def renumber(line, indeces, atoms):
    if indeces == 0:
        return line
    renumbering = line[:6 * indeces]
    line_split = list(filter(lambda x: x != '', renumbering.split(' ')))
    line_atoms = list(map(int, line_split))
    renumberer = renumber_one(atoms)
    renumbering = list(map(renumberer, line_atoms))
    new_nums = ' '.join(['%5d' % num for num in renumbering])
    return new_nums + line[6 * indeces + 1:]


def main(args):
    current_num_of_indices = 0
    input_file_name = args.input
    output_file_name = args.output
    remove_atoms = []
    for r in args.range:
        r = list(map(int, r.split('-')))
        remove_atoms += r if len(r) == 1 else list(range(r[0], r[1] + 1))
    atom_checker = remove_from_list(remove_atoms)

    with open(input_file_name) as input_file, open(output_file_name, 'w') as output_file:
        for line in input_file:
            # ignore comments or empty lines
            if line[0] == ';' or line[0] == '\n':
                output_file.write(line)
                continue
            # keep system's definitions 
            if '[' in line:
                current_num_of_indices = NUMBER_OF_INDICES.get(line.split(' ')[1], 0)
                output_file.write(line)
                continue

            check = list(filter(lambda x: x != '', line.split(' ')))
            check = list(map(int, check[:current_num_of_indices]))
            if atom_checker(check):
                continue

            output_file.write(renumber(line, current_num_of_indices, remove_atoms))


if __name__ == '__main__':
    main(parse_args())

