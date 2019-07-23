import argparse
import math
import matplotlib.pyplot as plt
import numpy as np


accept_residues = ('CNT', 'DA', 'DT', 'DG', 'DC')
default_cutoff = 0.45


def parse_args(*argument_array):
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('definition', default=1, type=int)
    parser.add_argument('--ignore-h', action='store_true')
    parser.add_argument('--cutoff', type=float, default=default_cutoff)
    return parser.parse_args(*argument_array)


def calculate_distance(c1, c2):
    return math.sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)


def calculate_nearby_atom_counts(cnt_bucket, dna_partition, cutoff):
    return sum([
        any([calculate_distance(cnt_atom, dna_atom) < cutoff
             for cnt_atom in cnt_bucket])
        for dna_atom in dna_partition
    ])


def calculate_nearby_atom_counts_by_partitions(cnt_bucket, dna_partitions, cutoff):
    return [calculate_nearby_atom_counts(cnt_bucket, dna_partition, cutoff) for dna_partition in dna_partitions]


def repartition_based_on_definition(bucket, definition):
    max_index = -1
    residue_index_map = {}
    partitions = []
    if definition == 1:
        for elem in bucket:
            if elem[0] not in residue_index_map:
                partitions.append([])
                max_index += 1
                residue_index_map[elem[0]] = max_index
            partitions[max_index].append(elem[1])
    else:
        partitions.append([])
        for elem in bucket:
            partitions[0].append(elem[1])
    return partitions


def calculate_q(cnt_bucket, dna_bucket, cutoff, definition):
    dna_partitions = repartition_based_on_definition(dna_bucket, definition)
    counts = calculate_nearby_atom_counts_by_partitions(cnt_bucket, dna_partitions, cutoff)
    if definition == 1:
        return sum([count > 0 for count in counts]) / len(counts)
    else:
        return counts[0] / len(dna_partitions[0])


def process_file(input_file_name, ignore_h):
    cnt_bucket = []
    dna_bucket = []
    with open(input_file_name) as input_file:
        input_file.readline()
        input_file.readline()
        for line in input_file:
            ignore = True
            components = list(filter(lambda x: x != '', line.split(' ')))
            for res in accept_residues:
                if components[0].find(res) != -1:
                    ignore = False
            if ignore or (ignore_h and components[1].find('H') != -1):
                continue
            coordinates = (float(components[3]), float(components[4]), float(components[5]))
            if components[0].find('CNT') != -1:
                cnt_bucket.append(coordinates)
            else:
                dna_bucket.append((components[0], coordinates))
    return cnt_bucket, dna_bucket


def main(args):
    input_file_name = args.input
    definition = args.definition
    cutoff = args.cutoff
    cutoff_array = list(np.arange(0.1, 2, 0.01))
    qs = []
    cnt_bucket, dna_bucket = process_file(input_file_name, args.ignore_h)
    for cutoff in cutoff_array:
        q = calculate_q(cnt_bucket, dna_bucket, cutoff, definition)
        qs.append(q)
    derivative_qs = []
    for ind in range(2, len(qs)):
        derivative_qs.append((qs[ind] - qs[ind - 2]) / 0.0001)
    # axes = plt.gca()
    # axes.set_ylim([0, 1])
    plt.plot(cutoff_array[2:], derivative_qs, marker='o')
    plt.show()


if __name__ == '__main__':
    main(parse_args())
