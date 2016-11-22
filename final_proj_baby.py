from os.path import join
import numpy as np
list_locations = []

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

def make_hash(ref):
    hash_ref = {}
    kmer_len = 12
    for i in range(len(ref)-kmer_len+1):
       kmer = ref[i:i+kmer_len]
       if kmer in hash_ref:
           hash_ref[kmer].append(i)
       else:
           hash_ref[kmer] = [i]
    return hash_ref

def split_read(read):
    splitsequences = []
    kmer_len = 12
    numsegments = int(round(len(read)/kmer_len))
    for i in range(0, len(read)-numsegments, kmer_len):
        splitsequences.append(read[i:i+kmer_len])
    return splitsequences

def read_reads(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    readcont = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
        x = PReads(paired_end_reads)
        readcont.append(x)
    return readcont


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
    for i in range(len(output_reference)):
        list_locations.append(i)
    return output_reference

def trivial_algorithm(pread_container, ref):
    read_alignment_locations = []
    output_read_pairs = []
    perfect_matches = []
    max_mismatches=3
    ref_dict = make_hash(ref)
    for pread in pread_container:
        pairs_to_append = []
        locations_container = []
        best_location = -1
        best_read = -1
        for single_read in pread.reads:
            best_location = -1
            best_read = -1
            split_reads = split_read(single_read)
            for cur_seg in split_reads:
                if cur_seg in ref_dict:
                    pos_in_ref = ref_dict[cur_seg]
                    for pos in pos_in_ref:
                        pos_in_read = single_read.find(cur_seg)
                        true_pos = pos-pos_in_read
                        current_best = calc_differences(single_read,ref[true_pos:true_pos+len(single_read)])
                        if (true_pos < 0) or (true_pos + len(single_read) > len(ref)):
                            continue
                        else:
                            if current_best > max_mismatches:
                                break
                            if best_location == -1 or calc_differences(single_read, ref[true_pos:true_pos+len(single_read)]) < current_best:
                                best_location = true_pos
                                best_read = single_read

            reverse_read = single_read[::-1]
            rev_split_reads = split_read(reverse_read)
            for cur_seg in rev_split_reads:
                if cur_seg in ref_dict:
                    pos_in_ref = ref_dict[cur_seg]
                    for pos in pos_in_ref:
                        pos_in_read = reverse_read.find(cur_seg)
                        true_pos = pos-pos_in_read
                        current_best = calc_differences(reverse_read,ref[true_pos:true_pos+len(single_read)])
                        if (true_pos < 0) or (true_pos + len(single_read) > len(ref)):
                            continue
                        else:
                            if current_best > max_mismatches:
                                break
                            if best_location == -1 or calc_differences(reverse_read, ref[true_pos:true_pos+len(single_read)]) < current_best:
                                best_location = true_pos
                                best_read = reverse_read
            pread.locations.append(best_location)
            for i in range(best_location, best_location+50):
                if i in list_locations:
                    list_locations.remove(i)
            if best_read == reverse_read:
                if single_read == pread.reads[0]:
                    pread.read0r = True
                elif single_read == pread.reads[1]:
                    pread.read1r = True
            pairs_to_append.append(best_read)
        #output_read_pairs.append(pairs_to_append)
        #read_alignment_locations.append(locations_container)
    return read_alignment_locations, output_read_pairs

if __name__ == "__main__":
    folder = 'practice_W_3'
    f_base = '{}_chr_1'.format(folder)
    reads_fn = join(folder, 'reads_{}.txt'.format(f_base))
    input_reads = read_reads(reads_fn)
    input_reads = input_reads[:1000]
    reference_fn = join(folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    alignments, reads = trivial_algorithm(input_reads, reference)
    output_fn = join(folder, "saved.txt")
    '''
    with open(output_fn, 'w') as output_file:
        for Read in input_reads:
            str1 = ',' + str(Read.read0r) + ',' + str(Read.read1r)
            str0 = Read.reads[0] + ',' + Read.reads[1] + ',' + str(Read.locations[0]) + ',' + str(Read.locations[1])
            toapp = str0 + str1
            output_file.write(toapp)
            output_file.write('\n')
    SNPstart = False
    snp_locs = []

    with open("snps_practice_E_1_chr_1.txt", 'r') as f:
        for line in f:
            if line.strip() == '>SNP':
                SNPstart = True
                continue
            if line.strip() == '>STR':
                SNPstart = False
                break
            if SNPstart:
                snp_locs.append(int(line.strip().split(',')[2]))
    #print "printingininginging SNPs"
    #for snp in snp_locs:
    #    print snp
    #for Read in input_reads:
    #    if Read.locations[0] in range(7500, 7600) or Read.locations[1] in range(7500, 7600):
    #        print Read.locations
    #        print Read.reads


        #if Read.read0r == Read.read1r:
        #    print Read.locations
        #    print Read.reads
        #if 1537 in Read.locations:
        #    print Read.read0r, Read.read1r
        #    print Read.locations
        #    print Read.reads
        #if Read.locations[0] == -1 and Read.locations[1] != -1:
        #    print Read.reads
        #    print Read.locations
        #if Read.locations[0] != -1 and Read.locations[1] == -1:
        #    print Read.reads
        #    print Read.locations
    #print "reference:"
    #print reference[1354:1404]
    #print list_locations
    '''
    current_run = 0
    poss_invs = []
    print list_locations

    for i in range(len(list_locations)-1):
        #if list_locations[i] == 7530:
        #    print "hey!"
        if list_locations[i] + 1 == list_locations[i+1]:
            current_run += 1
        else:
            if current_run > 15:    #min inv size
                loc_to_append = list_locations[i - current_run +1]
                seq_to_append = reference[loc_to_append:loc_to_append+current_run]
                poss_invs.append((list_locations[i - current_run + 1], reference[list_locations[i - current_run + 1]:list_locations[i - current_run + 1]+current_run]))
                current_run = 0
            current_run = 0
    toremove = []
    '''
    for poss in poss_invs:
        for snp in snp_locs:
            if snp in range(int(poss[0]),int(poss[0])+len(poss[1])):
                toremove.append(poss)
    #print "possible"
    #print poss_invs
    #print "removing:"
    #print toremove
    '''
    call_invs = [x for x in poss_invs if x not in toremove]
    output_invs = join(folder, "output_invs_5x.txt")
    with open(output_invs, 'w') as opi:
        for inv in call_invs:
            print inv
            opi.write(inv[1] + ',' + str(inv[0]) + '\n')
    print "call inversions:"
    print call_invs
    print "poss invs"
    print poss_invs
    print "len inv" + str(len(poss_invs))
    #for i in range(len(call_invs)):
    #    print poss_invs[i][::-1]
        #print "len of poss: " + str(len(poss_invs[i][1]))
    #print "num of poss_inv = " + str(len(poss_invs))
    #print list_locations
    #for i in range(len(poss_invs)):
    #    for Read in input_reads:
    #        if poss_invs[i][0] in range()
    #print "counter: " + str(counter)
    #print "total reads: " + str(len(input_reads))
    #sortaligns = sorted(alignments, key=lambda x: x[0])

