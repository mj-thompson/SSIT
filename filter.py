inv_container = []
present = []
hashdict = {}
is_1mil = True

with open("poss_sols.txt", 'r') as file:
    for line in file:
        if line[0] == '>':
            continue
        cur = line.strip()
        p1 = line.split(',')[0]
        p2 = line.split(',')[1]
        inv = (str(p1), int(p2))
        inv_container.append(inv)
totalsize = 0
for inv in inv_container:
    totalsize += len(str(inv[0]))
print "total size of inversions: " + str(totalsize)
#search for all that are less than 50, then look for them in reads
#ACAAGACTATTAGCGTA,850233
#TTGATACGCCGGCAACGACTCCCTTC,769738
#GCTGACTAGACTATGCTCCAAA,428412
#get rid of them if they're all within 10kof each other

less50 = [x for x in inv_container if len(x[0])<51]

with open("saved.txt", 'r') as file:
    for line in file:
        cont = line.split(',')
        if int(cont[2]) == -1:
            hashdict[cont[0]] = -1
        if int(cont[3]) == -1:
            hashdict[cont[1]] = -1
lookup = hashdict.keys()
for x in less50:
    for item in lookup:
        if item.find(x[0]) != -1:
            present.append(x) #if it is in our reads
        if item[::-1].find(x[0]) != -1:
            present.append(x)
#THESE ARE CANDIDATES TO REMOVE!!!
#that will filter out some of the false positives
#print candidates
#printr = set(removes)
save = [x for x in less50 if x not in present]
print "present in reads: "
print set(present)
print "not in reads: "
print set(save)


'''
if there are a ton of false positives, it's likely that there is enough data in reads that the inversion
exists in the donor, it just had too many mismatches and didn't map,
but if there aren't too many false positives, the holes in the sequence more likely weren't well covered
by the reads, and we should filter them out rather than kep them
'''
#if there's greater than 5% discovery rate check the reads as a filter
#if is_1mil:
with open("inv_filtered.txt", 'w') as opi:
    opi.write(">INV" + '\n')
    for inv in set(present):
        opi.write((inv[0]) + ',' + str(inv[1]) + '\n')
#else:
with open("inv_filtered1.txt", 'w') as opi:
    opi.write(">INV" + '\n')
    for inv in set(save):
        opi.write(str(inv[0]) + ',' + str(inv[1]) + '\n')
