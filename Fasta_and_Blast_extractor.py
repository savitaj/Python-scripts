#Import Biopython modules after installation "pip install biopython'
from Bio import SeqIO
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 

#Open a fasta file and parse the headers and sequence 
with open("Example.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print("Id: %s" % record.id) 
        print("Name: %s" % record.name) 
        print("Description: %s" % record.description) 
        print("Annotations: %s" % record.annotations) 
        print("Sequence Data: %s" % record.seq)
        
#blast using the gi accession
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastp", "nr", "1796318600") #SARS-CoV-2 envelope

#result_handle2 = NCBIWWW.qblast("blastp", "p", "1798174255") SARS-CoV-2 nucleocapsid
#result_handle2 = NCBIWWW.qblast("blastp", "p", "1796318601") SARS-CoV-2 membrane glycoprotein

with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

from Bio.Blast import NCBIXML
result=open("my_blast.xml","r")
records= NCBIXML.parse(result)
item=next(records)
for alignment in item.alignments:
          for hsp in alignment.hsps:
                 if hsp.expect <0.01:
                         print('****Alignment****')
                         print('sequence:', alignment.title) 
                         print('length:', alignment.length)
                         print('score:', hsp.score)
                         print('gaps:', hsp.gaps)
                         print('e value:', hsp.expect)
                         print(hsp.query[0:90] + '...')
                         print(hsp.match[0:90] + '...')
                         print(hsp.sbjct[0:90] + '...')
