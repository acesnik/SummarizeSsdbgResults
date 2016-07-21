#! /usr/bin/python

# Combines GetAveragePositions.py, SeqVarFdrs.py, and CountSeqVars.py to give
# a comprehensive summary of combined seq variant and GPTM approaches.
#
# SAV peptides: to evaluate whether they are tryptic if the first amino acid
# of the SAV peptide is used, the program looks at the pep.all.fasta file to
# reference the preceding amino acid.

import sys
import utility

__author__ = "anthonycesnik"
__date__ = "$Apr 29, 2015 8:24:16 PM$"


# Evaluate arbitrary ends of peptide for being tryptic
def sap_eval(descrips, line, allPep):
  accession = line[13].split('|')[1]
  decoy = accession[:6] == 'DECOY_'
  if decoy: accession = accession[6:]

  if len(descrips) == 16: sap = descrips[7].split(':') #snpefftopeptides
  elif len(descrips) == 12: sap = descrips[5].split(':') #ProteogenomicDBGenerator
  else:
    print "Error: SAV description not recognized\n" + str(line)
    exit(2)

  sapAAPos = int(sap[1][1:-1])
  aaSeq = ''
  # Evaluate arbitrary ends of peptide for being tryptic
  trypBeg, trypEnd = False, int(line[15]) != 67 or line[12][-1] in ['K','R']
  if int(line[14]) != 1: trypBeg = True
  elif not allPep[accession]: print "Did not find the reference sequence for " + accession + "\nPlease check that the pep.all file is the correct version."
  else: 
    aaSeq = allPep[accession]
    if decoy: aaSeq = aaSeq[::-1] if aaSeq[0] != 'M' else aaSeq[0]+aaSeq[:0:-1] #reverse the seq; leave M at beginning if there, as specified in Morpheus paper
    if sapAAPos <= 34 or aaSeq[sapAAPos-34] in ['K','R']: trypBeg = True #Either starts at beg of protein or beg is tryptic
    else: trypBeg = False
  return trypBeg, trypEnd


def __main__():
  USAGE = "python SummarizeMorpheusResults.py outFile pep.all.fasta listOfFolders"
  if len(sys.argv) < 4: print USAGE
  else:
    folderList = sys.argv[3:]
    allPep = utility.read_fasta()

    try:
      print "opening outFile " + sys.argv[1]
      out_file = open(sys.argv[1], 'w')
      out_file.write("\t\t\tcountseqvar\t\t\t\t\t\t\t\tseqvarfdr\t\t\t\t\n")
      out_file.write("folder\tline\tsapD\tTotalPsmsNoRaggedEnds\tTotalUniqPepsNoRaggedEnds\tUniProtPsms\tUniProtUniqPeps\tSapPsms\tSapUniqPeps\tJuncPsms\tJuncUniqPeps\tTotalFdr\tUniProtFdr\tSapFdr\tJuncFdr\n")
      for folder in folderList:
        allTarget,allDecoy,allFdr,uniTarget,uniDecoy,sapTarget,sapDecoy,juncTarget,juncDecoy = 0,0,0,0,0,0,0,0,0 #seqvarfdr
        totalPsms,totalPeps,uniPsms,uniPeps,sapPsms,sapPeps,juncPsms,juncPeps = 0,0,0,0,0,0,0,0 #countseqvars
        lineCt = 0
        invisSap = ['I','L'] #the leucine - isoleucine transition is invisible to mass spec

        #Print the folder and room for line and cutoff info
        out_file.write(folder + "\t\t\t")

        #Count sequence variants in the unique peptides folder
        print "opening " + folder + "/unique_peptides.tsv"
        uniqPeps = open(folder + "/unique_peptides.tsv", 'r')
        with open(folder + "/unique_peptides.tsv", 'r') as uniqPeps:
          for line in uniqPeps:
            line = line.split('\t')
            if line[0][:8] == 'Filename': continue
            elif float(line[30]) < 1 and line[26] == 'True':
                descrips = line[13].split('|')
                descrips = descrips[2].split(' ')
                if descrips[0] in ['pep:sap','pep:sav']:
                  if len(descrips) == 16: sap = descrips[7].split(':') #snpefftopeptides
                  elif len(descrips) == 12: sap = descrips[5].split(':') #ProteogenomicDBGenerator
                  else:
                    print "Error: SAV description not recognized\n" + str(line)
                    exit(2)
                  refAA, altAA = sap[1][0], sap[1][-1]
                  #Evaluate arbitrary ends of peptide for being tryptic
                  trypBeg, trypEnd = sap_eval(descrips, line, allPep)
                  #If it's tryptic, decide whether to count it as a sap peptide
                  if trypBeg and trypEnd:
                    totalPeps += 1
                    if int(line[14]) <= 34 and int(line[15]) >= 34 and not (refAA in invisSap and altAA in invisSap): sapPeps += 1 #must contain sap position and not be an I<->L transition
                elif line[11].find("UniProt:") != -1: uniPeps += 1
                elif descrips[0] == 'pep:splice': juncPeps += 1
                else: totalPeps += 1

        #Work with PSMs folder for each metric
        #BUG SOLVED: not reading the whole psms file because the original seqvarfdr code had a break in it that was being used before the
        #end of morpheus's 1% FDR.
        print "opening " + folder + "/PSMs.tsv"
        with open(folder + "/PSMs.tsv", "r") as psms:
          for line in psms:
            line = line.split('\t')
            if line[0][:8] == 'Filename': continue

            descrips = line[13].split('|')
            descrips = descrips[2].split(' ')

            if float(line[30]) < 1:
              if line[26] == 'True':
                #countseqvars
                if descrips[0] in ['pep:sap','pep:sav']:
                  if len(descrips) == 16: sap = descrips[7].split(':') #snpefftopeptides
                  elif len(descrips) == 12: sap = descrips[5].split(':') #ProteogenomicDBGenerator
                  else:
                    print "Error: SAV description not recognized\n" + str(line)
                    exit(2)
                  refAA, altAA = sap[1][0], sap[1][-1]
                  #Evaluate arbitrary ends of peptide for being tryptic
                  trypBeg, trypEnd = sap_eval(descrips, line, allPep)
                  #If it's tryptic, decide whether to count it as a sap peptide
                  if trypBeg and trypEnd:
                    totalPsms += 1
                    allTarget += 1
                    if int(line[14]) <= 34 and int(line[15]) >= 34 and not (refAA in invisSap and altAA in invisSap): #must contain sap position, and must not have L<->I transition
                        sapPsms += 1
                        sapTarget += 1
                elif line[11].find("UniProt:") != -1:
                  uniPsms += 1
                  uniTarget += 1
                elif descrips[0] == 'pep:splice':
                  juncPsms += 1
                  juncTarget += 1
                else:
                  totalPsms += 1
                  allTarget += 1
              else:
                if descrips[0] in ['pep:sap','pep:sav']:
                  if len(descrips) == 16: sap = descrips[7].split(':') #snpefftopeptides
                  elif len(descrips) == 12: sap = descrips[5].split(':') #ProteogenomicDBGenerator
                  else:
                    print "Error: SAV description not recognized\n" + str(line)
                    exit(2)
                  refAA, altAA = sap[1][0], sap[1][-1]
                  #Evaluate arbitrary ends of peptide for being tryptic
                  trypBeg, trypEnd = sap_eval(descrips, line, allPep)
                  #If it's tryptic, decide whether to count it as a sap peptide
                  if trypBeg and trypEnd:
                    allDecoy += 1
                    if int(line[14]) <= 34 and int(line[15]) >= 34 and not (refAA in invisSap and altAA in invisSap): sapDecoy += 1 #must contain sap position
                elif line[11].find("UniProt:") != -1: uniDecoy += 1
                elif descrips[0] == 'pep:splice': juncDecoy += 1
                else: totalPeps += 1
            else: break

        #write to output: sumAggregate, countseqvar
        out_file.write(str(totalPsms) + '\t' + str(totalPeps) + '\t' + str(uniPsms) + '\t' + str(uniPeps) + '\t' + str(sapPsms) + '\t' + str(sapPeps) + '\t' + str(juncPsms) + '\t' + str(juncPeps) + '\t')

        #seqvar write
        try: allFdr = float(allDecoy)/float(allTarget) * 100
        except ZeroDivisionError: allFdr = 'N/A'
        try: uniFdr = float(uniDecoy)/float(uniTarget) * 100
        except ZeroDivisionError: uniFdr = 'N/A'
        try: sapFdr = float(sapDecoy)/float(sapTarget) * 100
        except ZeroDivisionError: sapFdr = 'N/A'
        try: juncFdr = float(juncDecoy)/float(juncTarget) * 100
        except ZeroDivisionError: juncFdr = 'N/A'
        out_file.write(str(allFdr) + '\t' + str(uniFdr) + '\t' + str(sapFdr) + '\t' + str(juncFdr) + '\n')

    #except(IndexError): print "error processing a line in PSMs file"
    except(IOError): print "Input/output error"

if __name__ == "__main__": __main__()
  

