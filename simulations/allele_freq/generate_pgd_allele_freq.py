#!/usr/bin/env python3

# generate the pgd files for the allele frequency simulation
# header file and divergence as args

import sys
header_file = sys.argv[1]
div = sys.argv[2]

freqs = range(0, 21)

x = 0

while x < 21:
    ind_1 = 20 - x
    ind_2 = x

    indiv = '{} '.format('1+' + str(div) + '%') * ind_1 + '{} '.format('2+' + str(div) + '%') * ind_2
    freq_line = '50000 {}'.format(indiv)

    fn = x / 20.0 * 100
    pgd = open('pgds/' + str(int(fn)) + '_freq_{}_div.pgd'.format(div), 'w+')

    header = open(header_file, 'r')
    for line in header:
        pgd.write(line)
    header.close()

    pgd.write(freq_line + '\n')
    pgd.close()

    print(freq_line)
    x += 1

