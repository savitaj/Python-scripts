#Import Biopython modules after installation "pip install biopython'
from Bio import SeqIO
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 

#Open a fasta file and parse the headers and sequence 
with open("SARS_CoV-2_UP000464024_Proteome.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print("Id: %s" % record.id) 
        print("Name: %s" % record.name) 
        print("Description: %s" % record.description) 
        print("Sequence Data: %s" % record.seq)
        
#blast using the gi accession
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastp", "nr", "1796318600") #SARS-CoV-2 envelope

#Can do similarly with other SARS-CoV2 proteins
#result_handle2 = NCBIWWW.qblast("blastp", "nr", "1798174255") SARS-CoV-2 nucleocapsid
#result_handle2 = NCBIWWW.qblast("blastp", "nr", "1796318601") SARS-CoV-2 membrane glycoprotein

# Write the blast results to my_blast.xml
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

# Read the blast output in my_blast.xml file and print the alignments
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
                        
# Set Evalue threshold and print header if expect value is less than threshold
E_VALUE_THRESH = 1e-20 
for record in NCBIXML.parse(open("my_blast.xml")): 
    if record.alignments: 
        print("\n") 
        print("query: %s" % record.query[:100]) 
        for align in record.alignments:
            for hsp in align.hsps: 
                if hsp.expect < E_VALUE_THRESH: 
                    print("match: %s " % align.title[:100])
