from PySide2 import*
from PySide2.QtWidgets import*
import sys
import finalproject

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio import Align
import math
import matplotlib.pyplot as plt
import numpy

app = QApplication(sys.argv)

class Main:
    ex = finalproject.Ui_MainWindow()
    w = QMainWindow()
    ex.setupUi(w)
    w.show()

    def __init__(self):
        self.ex.entrez_fetch.clicked.connect(self.entrez_record)
        self.ex.entrez_info_2.clicked.connect(self.entrez_info)
        self.ex.entrez_save.clicked.connect(self.entrez_save)
        self.ex.importbutton.clicked.connect(self.import_record)
        self.ex.analyzebutton.clicked.connect(self.analyze)
        self.ex.aminoacidbutton.clicked.connect(self.aminodetails)
        self.ex.savebutton.clicked.connect(self.save)
        self.ex.clearbutton.clicked.connect(self.clear)
        self.ex.import_seq1.clicked.connect(self.import_seq1)
        self.ex.import_seq2.clicked.connect(self.import_seq2)
        self.ex.align.clicked.connect(self.align)
        self.ex.pairwise_save.clicked.connect(self.pairwise_save)
        self.ex.pmid_info.clicked.connect(self.entrez_info)
        self.ex.pmid_fetch.clicked.connect(self.pmid)
        self.ex.pmid_save.clicked.connect(self.pmid_save)

# -----------------------------------------Default Button--------------------------------------#

        self.ex.aminoacidbutton.setEnabled(False)
        self.ex.savebutton.setEnabled(False)
        self.ex.entrez_save.setEnabled(False)
        self.ex.pairwise_save.setEnabled(False)
        self.ex.pmid_save.setEnabled(False)

#-----------------------------------------Entrez Nucleotide--------------------------------------#

    def entrez_record(self):
        id = self.ex.entrez_id.toPlainText()
        email = self.ex.entrez_email.toPlainText()
        spinner = self.ex.spinner.currentText()

        try:
            Entrez.email = str(email)
            handle = Entrez.efetch(db='nucleotide', id=str(id), rettype=str(spinner), retmode='text')
            self.ex.entrez_info.setText(handle.read())
            self.ex.entrez_save.setEnabled(True)

        except:
            n = QWidget()
            n.title = "ID not recognized"
            n.left = 10
            n.top = 10
            n.width = 100
            n.height = 100
            fail = "The ID you inserted is not registered in NCBI nucleotide database"
            reply = QMessageBox.information(n, n.title,
                                            fail)

    def entrez_info(self):
        n = QWidget()
        n.title = "Important Information"
        n.left = 10
        n.top = 10
        n.width = 100
        n.height = 100
        info = "NCBI will use this email address provided to contact you if there is a problem"
        reply = QMessageBox.information(n, n.title,
                                        info)

        n = QWidget()
        n.title = "Important Information"
        n.left = 10
        n.top = 10
        n.width = 100
        n.height = 100
        info = "In case of excessive usage, NCBI will attempt to contact a user at the email address provided prior to blocking access to E-utilities"
        reply = QMessageBox.information(n, n.title,
                                        info)

        n = QWidget()
        n.title = "WARNING"
        n.left = 10
        n.top = 10
        n.width = 100
        n.height = 100
        warning = "Please DO NOT use RANDOM email, it is better not to give an email at all !"
        reply = QMessageBox.information(n, n.title,
                                        warning)


    def entrez_save(self):
        id = self.ex.entrez_id.toPlainText()
        info = self.ex.entrez_info.toPlainText()
        filename = QFileDialog.getSaveFileName(self.w, 'Save As', str(id), "All Files (*);;Text  Files(*.txt);;Genbank Files(*.gbk);;FASTA (*.fasta)")
        if filename:
            dir = filename[0].replace('/', '\\')
            file = open(str(dir), 'w')
            file.write(str(info))
            file.close()
        self.ex.entrez_save.setEnabled(False)

# -----------------------------------------Sequence Parsing--------------------------------------#

    def import_record(self):
        filename = QFileDialog.getOpenFileName(self.w, 'Open', "", "All Files (*);;Text  Files(*.txt);;Genbank Files(*.gbk);;FASTA (*.fasta)")
        if filename:
            dir = filename[0]
            file = open(str(dir), 'r')
            self.ex.insertbox.setText(str(dir))

    def analyze(self):
        record = self.ex.insertbox.text()
        if 'fasta' in record:
            for seq_record in SeqIO.parse(record, "fasta"):
                self.ex.id.setText(seq_record.id)
                self.ex.description.setText(seq_record.description)
                self.ex.seqalphabet.setText(repr(seq_record.seq))
                self.ex.recordsequence.setText(str(seq_record.seq))
                self.ex.sequencelength.setText(str(len(seq_record)))
                my_seq = Seq(str(seq_record.seq))
                complement = my_seq.complement()
                rev_complement = my_seq.reverse_complement()
                self.ex.complementsequence.setText(str(complement))
                self.ex.reversecomplement.setText(str(rev_complement))
            self.ex.aminoacidbutton.setEnabled(True)
            self.ex.savebutton.setEnabled(True)
            seq_length = self.ex.sequencelength.toPlainText()
            total_amino = round(int(seq_length)/3)
            self.ex.totalaminoacid.setText(str(total_amino))


        elif 'gbk' in record:
            for seq_record in SeqIO.parse(record, "genbank"):
                self.ex.id.setText(seq_record.id)
                self.ex.description.setText(seq_record.description)
                self.ex.seqalphabet.setText(repr(seq_record.seq))
                self.ex.recordsequence.setText(str(seq_record.seq))
                self.ex.sequencelength.setText(str(len(seq_record)))
                my_seq = Seq(str(seq_record.seq))
                complement = my_seq.complement()
                rev_complement = my_seq.reverse_complement()
                self.ex.complementsequence.setText(str(complement))
                self.ex.reversecomplement.setText(str(rev_complement))
                self.ex.aminoacidbutton.setEnabled(True)
                self.ex.savebutton.setEnabled(True)
            seq_length = self.ex.sequencelength.toPlainText()
            total_amino = round(int(seq_length) / 3)
            self.ex.totalaminoacid.setText(str(total_amino))

        else:
            n = QWidget()
            n.title = "Data Not Recognized"
            n.left = 10
            n.top = 10
            n.width = 100
            n.height = 100
            fail = "Data Not Recognized!"
            reply = QMessageBox.information(n, n.title,
                                         fail)
            self.ex.aminoacidbutton.setEnabled(False)

    def aminodetails(self):
        main = []
        id = self.ex.id.toPlainText()
        record = self.ex.description.toPlainText()
        amino = self.ex.recordsequence.toPlainText()
        for i in range(0, len(amino),3):
            main.append(amino[i:i+3])

        #---------------------------------------#
        P1 = main.count('TTT')
        P2 = main.count('TTC')
        Phenylalanine = (P1+P2)
        # ---------------------------------------#
        L1 = main.count('TTA')
        L2 = main.count('TTG')
        L3 = main.count('CTT')
        L4 = main.count('CTC')
        L5 = main.count('CTA')
        L6 = main.count('CTG')
        Leucine = (L1+L2+L3+L4+L5+L6)
        # ---------------------------------------#
        I1 = main.count('ATT')
        I2 = main.count('ATC')
        I3 = main.count('ATA')
        Isoleucine = (I1+I2+I3)
        # ---------------------------------------#
        M1 = main.count('ATG')
        Methionine = M1
        # ---------------------------------------#
        V1 = main.count('GTT')
        V2 = main.count('GTC')
        V3 = main.count('GTA')
        V4 = main.count('GTG')
        Valine = (V1+V2+V3+V4)
        # ---------------------------------------#
        S1 = main.count('TCT')
        S2 = main.count('TCC')
        S3 = main.count('TCA')
        S4 = main.count('TCG')
        S5 = main.count('AGT')
        S6 = main.count('AGC')
        Serine = (S1+S2+S3+S4+S5+S6)
        # ---------------------------------------#
        P1 = main.count('CCT')
        P2 = main.count('CCC')
        P3 = main.count('CCA')
        P4 = main.count('CCG')
        Proline = (P1+P2+P3+P4)
        # ---------------------------------------#
        T1 = main.count('ACT')
        T2 = main.count('ACC')
        T3 = main.count('ACA')
        T4 = main.count('ACG')
        Threonine = (T1+T2+T3+T4)
        # ---------------------------------------#
        A1 = main.count('GCT')
        A2 = main.count('GCC')
        A3 = main.count('GCA')
        A4 = main.count('GCG')
        Alanine = (A1+A2+A3+A4)
        # ---------------------------------------#
        Y1 = main.count('TAT')
        Y2 = main.count('TAC')
        Tyrosine = (Y1+Y2)
        # ---------------------------------------#
        Stop1 = main.count('TAA')
        Stop2 = main.count('TAG')
        Stop3 = main.count('TGA')
        Stop = (Stop1+Stop2+Stop3)
        # ---------------------------------------#
        H1 = main.count('CAT')
        H2 = main.count('CAC')
        Histidine = (H1+H2)
        # ---------------------------------------#
        Q1 = main.count('CAA')
        Q2 = main.count('CAG')
        Glutamine = (Q1+Q2)
        # ---------------------------------------#
        N1 = main.count('AAT')
        N2 = main.count('AAC')
        Asparagine = (N1+N2)
        # ---------------------------------------#
        K1 = main.count('AAA')
        K2 = main.count('AAG')
        Lysine = (K1+K2)
        # ---------------------------------------#
        D1 = main.count('GAT')
        D2 = main.count('GAC')
        Glutamic_acid = (D1+D2)
        # ---------------------------------------#
        E1 = main.count('GAA')
        E2 = main.count('GAG')
        Aspartic_acid = (E1+E2)
        # ---------------------------------------#
        C1 = main.count('TGT')
        C2 = main.count('TGC')
        Cysteine = (C1+C2)
        # ---------------------------------------#
        T1 = main.count('TGG')
        Tryptophan = T1
        # ---------------------------------------#
        R1 = main.count('CGT')
        R2 = main.count('CGC')
        R3 = main.count('CGA')
        R4 = main.count('CGG')
        R5 = main.count('AGA')
        R6 = main.count('AGG')
        Arginine = (R1+R2+R3+R4+R5+R6)
        # ---------------------------------------#
        G1 = main.count('GGT')
        G2 = main.count('GGC')
        G3 = main.count('GGA')
        G4 = main.count('GGG')
        Glycine = (G1+G2+G3+G4)

        labels = ['Phe - ' +str(Phenylalanine),'Leu - '+str(Leucine),'Ile - '+str(Isoleucine),
                  'Met - '+str(Methionine),'Val - '+str(Valine),'Ser - '+str(Serine),'Pro - '+str(Proline),'Thr - '+str(Threonine),
                  'Ala - '+str(Alanine),'Tyr - '+str(Tyrosine),'Stop - '+str(Stop),'His - '+str(Histidine),'Glu - '+str(Glutamine),
                  'Asp - '+str(Asparagine),'Lys - '+str(Lysine),'Glu - '+str(Glutamic_acid), 'Asp - '+str(Aspartic_acid),
                  'Cys - '+str(Cysteine),'Try - '+str(Tryptophan),'Arg - '+str(Arginine),'Gly - '+str(Glycine)]

        sizes = [Phenylalanine, Leucine, Isoleucine, Methionine, Valine, Serine, Proline, Threonine, Alanine, Tyrosine, Stop, Histidine, Glutamine, Asparagine, Lysine, Glutamic_acid, Aspartic_acid, Cysteine, Tryptophan, Arginine, Glycine]

        #height = [Phenylalanine, Leucine, Isoleucine, Methionine, Valine, Serine, Proline, Threonine, Alanine, Tyrosine,
                 # Stop, Histidine, Glutamine, Asparagine, Lysine, Glutamic_acid, Aspartic_acid, Cysteine, Tryptophan,
                 # Arginine, Glycine]

        #bottomlabels = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Tyr', 'Stop', 'His', 'Glu', 'Asp', 'Lys', 'Glu', 'Asp', 'Cys', 'Try', 'Arg', 'Gly']



        # id_title = self.ex.id.toPlainText()
        # main = ' Amino Acid Details'

        ####-----------------PIE CHART---------------------####
        colors = ['lightgreen', 'lightblue']
        plt.pie(sizes, labels=labels, colors=colors,
                 autopct='%1.1f%%', shadow=True, startangle=140)

        ####-----------------BAR CHART---------------------####
        # y = numpy.arange(len(labels))
        # plt.bar(y, sizes, alpha = 1.0)
        # plt.xticks(y, labels)
        # plt.ylabel = 'Quantity'

        plt.axis('equal')
        plt.title('Amino Acid Details of '+record+' ('+id+')')
        plt.show()

    def save(self):
        id = self.ex.id.toPlainText()
        desc = self.ex.description.toPlainText()
        seqalphabet = self.ex.seqalphabet.toPlainText()
        seq = self.ex.recordsequence.toPlainText()
        length = self.ex.sequencelength.toPlainText()
        totaminoacid = self.ex.totalaminoacid.toPlainText()
        complement = self.ex.complementsequence.toPlainText()
        revcomplement = self.ex.reversecomplement.toPlainText()

        filename = QFileDialog.getSaveFileName(self.w, 'Save As', "Analysis of "+str(id), "All Files (*);;Text  Files(*.txt);;Genbank Files(*.gbk);;FASTA (*.fasta)")
        if filename:
            dir = filename[0].replace('/', '\\')
            file = open(str(dir), 'w')
            file.write('                                                Analysis of '+id+'\n'+'\n')
            file.write("ID : "+id+'\n'+'\n')
            file.write("Description : "+desc+'\n'+'\n')
            file.write("Sequence & Alphabet : "+seqalphabet+'\n'+'\n')
            file.write("Sequence : "+seq+'\n'+'\n')
            file.write("Length : "+length+'\n'+'\n')
            file.write("Total Amino Acid : "+totaminoacid+'\n'+'\n')
            file.write("Complement of Sequence : "+complement+'\n'+'\n')
            file.write("Reverse Complement of Sequence : "+revcomplement+'\n'+'\n')
            file.close()

    def clear(self):
        self.ex.aminoacidbutton.setEnabled(False)
        self.ex.savebutton.setEnabled(False)

# -----------------------------------------Multiple Sequence Alignment--------------------------------------#

    def import_seq1(self):
        filename = QFileDialog.getOpenFileName(self.w, 'Open', "", "All Files (*);;Text  Files(*.txt);;Genbank Files(*.gbk);;FASTA (*.fasta)")
        dir = filename[0].replace('/', '\\')
        file = str(dir)
        self.ex.seq1.setText(file)

    def import_seq2(self):
        filename = QFileDialog.getOpenFileName(self.w, 'Open', "",
                                               "All Files (*);;Text  Files(*.txt);;Genbank Files(*.gbk);;FASTA (*.fasta)")
        dir = filename[0].replace('/', '\\')
        file = str(dir)
        self.ex.seq2.setText(file)

    def align(self):
        sequence_1 = self.ex.seq1.toPlainText()
        sequence_2 = self.ex.seq2.toPlainText()
        try:
            if 'fasta' in sequence_1 and sequence_2:
                seq1 = SeqIO.read(str(sequence_1), 'fasta')
                seq2 = SeqIO.read(str(sequence_2), 'fasta')
                alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
                result = str(alignments[0])
                var_alignment = str(len(alignments))
                list_alignment = str(alignments)
                self.ex.alignment_result.setText(result)
                self.ex.variable_alignment.setText(var_alignment)
                self.ex.list_alignment.setText(list_alignment)
                self.ex.pairwise_save.setEnabled(True)

            elif 'gbk' in sequence_1 and sequence_2:
                seq1 = SeqIO.read(str(sequence_1), 'gb')
                seq2 = SeqIO.read(str(sequence_2), 'gb')
                alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
                result = str(alignments[0])
                var_alignment = str(len(alignments))
                list_alignment = str(alignments)
                self.ex.alignment_result.setText(result)
                self.ex.variable_alignment.setText(var_alignment)
                self.ex.list_alignment.setText(list_alignment)
                self.ex.pairwise_save.setEnabled(True)

            else:
                seq1 = str(sequence_1).upper()
                seq2 = str(sequence_2).upper()
                aligner = Align.PairwiseAligner()
                alignments = aligner.align(seq1, seq2)
                self.ex.alignment_result.setText(str(alignments[0]))
                self.ex.variable_alignment.setText(str(len(alignments)))
                for i in alignments:
                    self.ex.list_alignment.append(str(i))
                self.ex.pairwise_save.setEnabled(True)

        except:
            n = QWidget()
            n.title = "Invalid Input !"
            n.left = 10
            n.top = 10
            n.width = 100
            n.height = 100
            fail = "Check the format file ! ensure the inputed file are correct !"
            reply = QMessageBox.information(n, n.title,
                                            fail)


    def pairwise_save(self):
        result = self.ex.alignment_result.toPlainText()
        var = self.ex.variable_alignment.toPlainText()
        list = self.ex.list_alignment.toPlainText()
        filename = QFileDialog.getSaveFileName(self.w, 'Save As', "",
                                               "All Files (*);;Text  Files(*.txt);;Genbank Files(*.gbk);;FASTA (*.fasta)")
        if filename:
            dir = filename[0].replace('/', '\\')
            file = open(str(dir), 'w')
            file.write('Alignment result : '+'\n')
            file.write(str(result))
            file.write('\n'+'\n')
            file.write('Variable alignment : '+'\n')
            file.write(str(var))
            file.write('\n'+'\n')
            file.write('List of Alignment : '+'\n')
            file.write(str(list))
        self.ex.pairwise_save.setEnabled(False)

# -----------------------------------------PubMed Lounge--------------------------------------#

    def pmid(self):
        id = self.ex.pmid.toPlainText()
        email = self.ex.pmid_email.toPlainText()
        spinner = self.ex.pmid_spinner.currentText()

        try:
            Entrez.email = str(email)
            handle = Entrez.efetch(db='pubmed', id=str(id), rettype=str(spinner), retmode='txt')
            self.ex.pmid_data.setText(handle.read())
            self.ex.pmid_save.setEnabled(True)

        except:
            n = QWidget()
            n.title = "PMID not recognized"
            n.left = 10
            n.top = 10
            n.width = 100
            n.height = 100
            fail = "The PMID you inserted is not registered in PubMed database"
            reply = QMessageBox.information(n, n.title,
                                            fail)

    def pmid_save(self):
        id = self.ex.pmid.toPlainText()
        info = self.ex.pmid_data.toPlainText()
        filename = QFileDialog.getSaveFileName(self.w, 'Save As', "PMID of "+str(id),
                                               "All Files (*);;Text  Files(*.txt);; XML Files(*.xml)")

        if filename:
            dir = filename[0].replace('/', '\\')
            file = open(str(dir), 'w')
            file.write(str(info))
            file.close()
        self.ex.pmid_save.setEnabled(False)


run = Main()
sys.exit(app.exec_())


