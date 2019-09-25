# Usage:
# python val2dot.py <Rmaps> <Valouev alignment results>

rmapfile = "ecoli-2000.valouev"
valouevfile="all-pairs-s21-7.txt"

id2name=dict()
name2id=dict()

fp = open(rmapfile, "r")
i = 0
j = 0
for line in fp:
    if i == 0:
        name = line.rstrip()
        id2name[j] = name
        name2id[name]=j
        j += 1
    i += 1
    if i > 2:
        i = 0

fp.close()

print("graph {")

fp = open(valouevfile, "r")
i = 0
for line in fp:
    if i == 0:
        cols = line.split()
        id1 = name2id[cols[0]]
        id2 = name2id[cols[1]]
        if id1 < id2:
            print(str(id1) + " -- " + str(id2) + ";")
        else:
            print(str(id2) + " -- " + str(id1) + ";")
    i += 1
    if i > 2:
        i = 0

print("}")
