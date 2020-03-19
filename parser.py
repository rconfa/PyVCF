
from VCFParser import VCFParser # import del mio parser

"""
print("parse file: example-4.2")
p = VCFParser('./Test/example-4.2.vcf')
p.parseFile()

print("\n\nparse file: example-4.1-info-multiple-values")
p = VCFParser('./Test/example-4.1-info-multiple-values.vcf')
p.parseFile()

print("\n\nparse file: 1kg.sites")
p = VCFParser('./Test/1kg.sites.vcf')
p.parseFile()


p = VCFParser('./Test/zipped.vcf.gz')
p.parseFile()
"""

p = VCFParser('./Test/myExample.vcf')
p.parseFile()