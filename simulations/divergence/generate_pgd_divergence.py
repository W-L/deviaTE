#!/usr/bin/env python3

# this script uses the population genome definition header-file (pgd_header)
# to produce the actual definition files with the different divergence levels

percentages = range(0, 31)

for i in percentages:
    perc_line = '50000 {}+{}%'.format('1', i)

    pgd = open('pgds/' + str(i) + '_div.pgd', 'w+')
    header = open('pgd_header', 'r')
    for line in header:
        pgd.write(line)
    header.close()

    pgd.write(perc_line + '\n')
    pgd.close()
