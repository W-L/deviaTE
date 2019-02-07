import sys
# generate pgd files for internal deletion validation
# one argument: headerfile

header = open(sys.argv[1])
te_seq = header.read().split('\n')[1].split('=')[1].split('"')[1]  # extract sequence
te_len = len(te_seq)
header.close()

freqs = range(0, 21, 2)
del_start = 1100
del_stop = [del_start + 100, del_start + 500, del_start + 1000, del_start + 2000]


for f in freqs:

    for stop in del_stop:
        t1 = f
        t2 = 20 - f

        pos = [del_start, stop]
        intdel = '1' + str(pos).replace(', ', '..') + '+'  # build deletion string

        truncline = '50000 ' + '{} '.format(intdel) * t1 + '{} '.format('1') * t2  # use different frequencies

        fn = f / 20.0 * 100  # only for filename - put freq and pos in name
        pgd = open('pgds/' + str(int(fn)) + '_freq_{}_{}.pgd'.format(pos[0], pos[1]), 'w+')

        header = open(sys.argv[1])
        for line in header:
            pgd.write(line)
        pgd.write(truncline + '\n')

        pgd.close()
        header.close()
        #print(truncline)


