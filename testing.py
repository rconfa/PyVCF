
from VCFParser import VCFParser # import del parser
filename = './Test/myExample.vcf' # file to be parsed
p = VCFParser(filename) # create the istance for parsing the file
p.parseFile() # parsing
