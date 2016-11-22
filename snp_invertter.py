
snp_cont = []
SNPs = True
with open("SNPpies.txt", 'r') as opfile:
    for line in opfile:
        if line[0] == '>':
            continue
        cur = line.strip()
        snp_loc = cur.split(',')[2]
        snp_cont.append(int(snp_loc))
streaks = []
with open("snp_consecs.txt") as opfile:
    for line in opfile:
        cur = line.strip()
        snp_loc = cur.split(',')[0]
        snp_amnt = cur.split(',')[1]
        streaks.append((int(snp_loc), int(snp_amnt)))
boundaries = []
#print streaks
#print snp_cont
for streak in streaks:
    mark = snp_cont.index(streak[0])
    markend = snp_cont[mark + streak[1]]
    boundaries.append((snp_cont[mark], markend))
print boundaries
ref = ''
with open("ref_hw2undergrad_E_2_chr_1.txt", 'r') as opfile:
    first_line = True
    for line in opfile:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        ref += line
poss_invs = []
for bound in boundaries:
    if bound[0] == 131987:
        print "what"
    if bound[1] - bound[0] < 15 or bound[1] - bound[0] > 40:
        continue
    refseg = ref[bound[0]:bound[1]]
    poss_invs.append((refseg, bound[0]))
with open("poss_sols.txt", 'w') as opfile:
    for poss in poss_invs:
        if poss[1] == 131987:
            print "hah?"
        opfile.write(poss[0] + "," + str(poss[1]) + "\n")



"""
#snp_cont = snp_cont[::-1]
start = -1
current_run = 0
run_start = True
potentials = []
for i in range(len(snp_cont)-1):
    if run_start:
        start = snp_cont[i]
        if snp_cont[i] in range(snp_cont[i+1]-16, snp_cont[i+1]):
            run_start = False
            current_run += 1
    else:
        if snp_cont[i] in range(snp_cont[i+1]-16, snp_cont[i+1]):
            current_run += 1
        else:
            if current_run > 0:
                potentials.append((start, current_run))
            current_run = 0
            run_start = True
with open("snp_consecs.txt", 'w') as opfile:
    for i in range(len(potentials)):
        if potentials[i][1] >= 5:
            opfile.write(str(potentials[i][0]) + ',' + str(potentials[i][1]) + '\n')
#131987 10
#201092 7
#253628 7
#332727 5
"""