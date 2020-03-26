
from VCFParser import VCFParser # import del mio parser
filename = './Test/myExample.vcf'
p = VCFParser(filename)
p.parseFile()
