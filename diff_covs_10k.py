from os.path import join
import numpy as np


class PReads:
    def __init__(self, rp):
        self.reads = rp
        self.locations = []
        self.read0r = False
        self.read1r = False

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

def findinvs_10(ref, reads, folder):
    list_unmapped = []
    list_locations = []
    for i in range(len(ref)):
        list_locations.append(int(i))
    for read in reads:
        for i in range(len(read.locations)):
            if read.locations[i] != -1:
                for j in range(read.locations[i], read.locations[i] + 50):
                    if j in list_locations:
                        list_locations.remove(j)
            elif read.locations[i] == -1:
                list_unmapped.append(read.reads[i])
    current_run = 0
    poss_invs = []
    for i in range(len(list_locations)-1):
        if list_locations[i] + 1 == list_locations[i+1]:
            current_run += 1
        else:
            if current_run > 15:    #min inv size
                loc_to_append = list_locations[i - current_run +1]
                seq_to_append = reference[loc_to_append:loc_to_append+current_run]
                poss_invs.append((int(list_locations[i - current_run + 1]), reference[list_locations[i - current_run + 1]:list_locations[i - current_run + 1]+current_run]))
                current_run = 0
            current_run = 0
    less50 = [x for x in poss_invs if len(x[1]) < 51]
    call_invs = []
    for x in less50:
        for item in list_unmapped:
            if item.find(x[1]) != -1:
                call_invs.append(x) #if it is in our reads
            if item[::-1].find(x[1]) != -1:
                call_invs.append(x)
    output_invs = join(folder, "output_invs10x.txt")
    with open(output_invs, 'w') as opi:
        opi.write(">INV" + '\n')
        for inv in poss_invs:
            print inv
            opi.write(str(inv[1]) + ',' + str(inv[0]) + '\n')

def findinvs_5(ref, reads, folder):
    list_unmapped = []
    list_locations = []
    for i in range(len(ref)):
        list_locations.append(i)
    for read in reads:
        for i in range(len(read.locations)):
            if read.locations[i] != -1:
                for j in range(read.locations[i], read.locations[i] + 50):
                    if j in list_locations:
                        list_locations.remove(j)
            elif read.locations[i] == -1:
                list_unmapped.append(read.reads[i])
    current_run = 0
    poss_invs = []
    for i in range(len(list_locations)-1):
        if list_locations[i] + 1 == list_locations[i+1]:
            current_run += 1
        else:
            if current_run > 15:    #min inv size
                loc_to_append = list_locations[i - current_run +1]
                seq_to_append = reference[loc_to_append:loc_to_append+current_run]
                poss_invs.append((int(list_locations[i - current_run + 1]), reference[list_locations[i - current_run + 1]:list_locations[i - current_run + 1]+current_run]))
                current_run = 0
            current_run = 0
    print poss_invs
    print type(poss_invs[0])
    print type(poss_invs[0][0])
    print type(poss_invs[0][1])
    less50 = [x for x in poss_invs if len(x[1]) < 51]
    call_invs = []
    for x in less50:
        for item in list_unmapped:
            if item.find(x[1]) != -1:
                call_invs.append(x) #if it is in our reads
            if item[::-1].find(x[1]) != -1:
                call_invs.append(x)
    output_invs = join(folder, "output_invs5x.txt")
    with open(output_invs, 'w') as opi:
        opi.write(">INV" + '\n')
        for inv in poss_invs:
            opi.write(str(inv[1]) + ',' + str(inv[0]) + '\n')


read_cont = []
folder = 'practice_W_3'
f_base = '{}_chr_1'.format(folder)
reference_fn = join(folder, 'ref_{}.txt'.format(f_base))
reference = read_reference(reference_fn)
reads = join(folder, "saved.txt")
with open(reads, 'r') as opfile:
    for line in opfile:
        cur = line.strip()
        cur = cur.split(',')
        x = PReads([cur[0], cur[1]])
        x.locations = [int(cur[2]), int(cur[3])]
        read_cont.append(x)
print "total reads: " + str(len(read_cont))
#for 10k is 3k, for 1mil it's 300,000
#try 5x and 10x coverage
#try the list index thing, make a list of size len(genome)
#remove the indices + 50
#filter
tenxcov = int(10*len(reference)/50)
fivexcov = int(5*len(reference)/50)
findinvs_10(reference, read_cont[:tenxcov], folder)
findinvs_5(reference, read_cont[:fivexcov], folder)

