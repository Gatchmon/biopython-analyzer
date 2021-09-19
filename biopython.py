# from Bio import Entrez
# Entrez.email = "stefanus.bernard@student.i3l.ac.id"
# handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="fasta", retmode="text")
# print(handle.read())

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
# Entrez.email = "stefanus.bernard@student.i3l.ac.id"
# handle = Entrez.efetch(db="nucleotide", id="EU392870", rettype="fasta", retmode="text")
# record = SeqIO.read(handle, "fasta")
# handle.close()
#
# print(record.name)
# print(record.description)
# print(str(record.seq))
#
# Entrez.email = "stefanus.bernard@student.i3l.ac.id"
# handle = Entrez.efetch(db="nucleotide", id="EU392870", rettype="gb", retmode="text")
# record = SeqIO.read(handle, "genbank")
# handle.close()
#
# print(record.name)
# print(record.description)
# print(str(record.seq))

# import os
# from Bio import SeqIO
# from Bio import Entrez
# Entrez.email = "stefanus.bernard@student.i3l.ac.id"
# filename = "EU490707.gbk"
# if not os.path.isfile(filename):
#     #Downloading...
#     net_handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype='gb', retmode='text')
#     out_handle = open(filename, 'w')
#     out_handle.write(net_handle.read())
#     out_handle.close()
#     net_handle.close()
#     print('saved')
#
# print('parsing...')
# record = SeqIO.read(filename, 'genbank')
# print(record)

# s="yaegfkegrfigfrhg;rtihgrtihg;otjgptj"
# for i in range(0,len(s),3):
#     print(s[i:i+3])

# amino = []
# for i in range(0, len(str(record.seq)),3):
#     print(str(record.seq)[i:i+3])
#     amino.append(str(record.seq)[i:i+3])
#
# print(amino)
# print(amino.count('ATG'))


# for i in str(record.seq):
#         print(i[0:3, 1])


# print(str(record.seq))
#
# print(str(record.seq).count('ATG'))
# my_seq = Seq('ATGCTAGCTAGCTAGCTACGAG')
# print(my_seq.complement())
# print(my_seq.reverse_complement())

# valine = ['GTT', 'GTG', 'GTA', 'GTC']
# a = str(record.seq).count('GTT')
# b = str(record.seq).count('GTG')
# c = str(record.seq).count('GTA')
# d = str(record.seq).count('GTC')
#
# print(
# str(record.seq).count('GTT'),
# str(record.seq).count('GTG'),
# str(record.seq).count('GTA'),
# str(record.seq).count('GTC'))
#
# print('Valine inside the whole sequences: ', a + b + c + d)
#
# import matplotlib.pyplot as plt
#
# # Data to plot
# labels = 'Python', 'C++', 'Ruby', 'Java'
# sizes = [215, 130, 245, 210]
# colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
# explode = (0, 0, 0, 0)  # explode 1st slice
#
# # Plot
# plt.pie(sizes, explode=explode, labels=labels, colors=colors,
#         autopct='%1.1f%%', shadow=True, startangle=140)
#
# plt.axis('equal')
# plt.show()

from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio import Align
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline


# seq1 = 'TAGCTATAGCTAGCA'
# seq2 = 'ATG'
# aligner = Align.PairwiseAligner()
# alignments = aligner.align(seq1, seq2)
# print(len(alignments))
#
# for i in alignments:
#     print(i)


#
#
# align = AlignIO.read('dengue.fasta', 'clustal')
# print(align)
# print(pairwise2.format_alignment(*alignments[2]))

# align = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
# print(len(align))
# print(pairwise2.format_alignment(*align[0]))


# seq1 = 'AAAAAAAPPPRRRSSSS'
# seq2 = 'RPAS'
# seq1.upper()
# seq2.upper()

#
#
# alignments = aligner.align(seq1, seq2)
# for x in alignments:
#     print(x)
#
# print('Total alignment: ' +str(len(alignments)))
# print('Score : '+str(score)+'\n')
#
# print('Adjustment by blosum62 and -10 (gap open penalty) and 0.5 (gap extension penalty)\n')
#
# align = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5)
# print(pairwise2.format_alignment(*align[0]))
# print('Total alignment: ' +str(len(align)))

#
Entrez.email = "stefanus.bernard@student.i3l.ac.id"
handle = Entrez.efetch(db="pubmed", id='30625488', rettype='XML', retmode='txt')
print(handle.read())
handle.close()





