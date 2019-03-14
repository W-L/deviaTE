from gtfIO import GTFReader, GTFWriter
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--input", dest="input", help="A gtf file containing the RepeatMasked gtf annotation")
parser.add_option("--minlen", dest="minleng", help="minimum length")
parser.add_option("--output", dest="output", help="A gtf output file")
parser.add_option("--maxdiv", dest="maxdiv", help="minimum length", default=99999999999999999)
(options, args) = parser.parse_args()

minleng = int(options.minleng)
w = GTFWriter(options.output)

for e in GTFReader(options.input):
    leng = (e.end - e.start) + 1
    if(leng >= minleng):
        if(e.score < float(options.maxdiv)):
            w.write(e)
w.close()







