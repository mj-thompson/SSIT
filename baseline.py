from os.path import join
import numpy as np


class PReads:
    def __init__(self, rp):
        self.reads = rp
        self.locations = []
        self.read0r = False
        self.read1r = False


def calc_differences(read1, read2):
    count = 0
    if len(read1) != len(read2):
        return np.inf
    for i in range(len(read1)):
        if read1[i] != read2[i]:
            count += 1
    return count


def read_reference(ref_fn):
    f = open(ref_fn, 'r')
    first_line = True
    output_reference = ''
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        output_reference += line  # We append each line to the output reference string.
    return output_reference


if __name__ == "__main__":
    folder = 'hw2undergrad_E_2'
    f_base = '{}_chr_1'.format(folder)
    reference_fn = join(folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    reads = join(folder, "saved.txt")
    read_cont = []
    both_rev = []
    with open(reads, 'r') as opfile:
        for line in opfile:
            if line[0] == '>':
                continue
            cont = line.strip()
            cont = cont.split(',')
            x = PReads([cont[0], cont[1]])
            x.locations = [int(cont[2]), int(cont[3])]
            if cont[4] == "True":
                x.read0r = True
            else:
                x.read0r = False
            if cont[5] == "False":
                x.read1r = False
            else:
                x.read1r = True
            read_cont.append(x)
    '''
    read_cont_10x = read_cont[:2000]
    read_cont_5x = read_cont[:1000]
    for read in read_cont:
        if -1 in read.locations:
            continue
        for i in range(len(read.reads)):
            loc = read.locations[i]
            diff = calc_differences(read.reads[i], reference[read.locations[i]:read.locations[i] + 50])
            diffr = calc_differences(read.reads[i][::-1], reference[read.locations[i]:read.locations[i] + 50])
            if diffr < diff:
                if i == 0:
                    read.read0r = True
                else:
                    read.read1r = True
                if read.read1r and read.read0r or (read.locations[0] != -1 and read.locations[1] != -1 and not read.read0r and not read.read1r):
                    both_rev.append(read)
    '''
    for read in read_cont:
        if read.read0r and read.read1r or (not read.read0r and not read.read1r and read.locations[0] != -1 and read.locations[1] != -1):
            both_rev.append(read)
    op = join(folder, "output_baseline.txt")
    with open(op, 'w') as opfile:
        for read in both_rev:
            if read.locations[0] > read.locations[1]:
                if read.locations[0] - read.locations[1] in range(0, 200):
                    opfile.write(read.reads[0] + ',' + str(read.locations[0]) + '\n' + read.reads[1] + "," + str(
                        read.locations[1]) + '\n')
            else:
                if read.locations[1] - read.locations[0] in range(0, 200):
                    opfile.write(read.reads[0] + ',' + str(read.locations[0]) + '\n' + read.reads[1] + "," + str(
                        read.locations[1]) + '\n')
        if len(both_rev) == 0:
            for read in read_cont_10x:
                if read.locations[0] == -1 and read.locations[1] != -1:
                    opfile.write(read.reads[0] + ',' + str(read.locations[0]) + '\n' + read.reads[1] + "," + str(
                        read.locations[1]) + '\n')
                if read.locations[0] != -1 and read.locations[1] == -1:
                    opfile.write(read.reads[0] + ',' + str(read.locations[0]) + '\n' + read.reads[1] + "," + str(
                        read.locations[1]) + '\n')
