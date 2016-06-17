#! /usr/bin/python

# Combines GetAveragePositions.py, SeqVarFdrs.py, and CountSeqVars.py to give
# a comprehensive summary of combined seq variant and GPTM approaches.
#
# SAV peptides: to evaluate whether they are tryptic if the first amino acid
# of the SAV peptide is used, the program looks at the pep.all.fasta file to
# reference the preceding amino acid.

__author__ = "anthonycesnik"
__date__ = "$Apr 29, 2015 8:24:16 PM$"

import sys

#Read in fasta file
def read_fasta():
  try:
    pepAll = open(sys.argv[2], 'r')
    allPep = {}
    header, seq = '',''
    starting = True
    count = 0
    for faLine in pepAll:
      if len(faLine.split(' ')) == 1: seq += faLine.strip()
      elif starting:
        header = faLine.split(' ')[0][1:]
        starting = False
      else:
        allPep[header] = seq
        count += 1
        if count % 10000 == 0: print "Read in " + str(count) + " fasta entries"
        seq = ''
        header = faLine.split(' ')[0][1:]
    pepAll.close()  
    return allPep
  except(IOError): print "IOError"

#Evaluate arbitrary ends of peptide for being tryptic
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
  #Evaluate arbitrary ends of peptide for being tryptic
  trypBeg, trypEnd = False, int(line[15]) != 67 or line[12][-1] in ['K','R']
  if int(line[14]) != 1: trypBeg = True
  elif not allPep[accession]: print "Did not find the reference sequence for " + accession + "\nPlease check that the pep.all file is the correct version."
  else: 
    aaSeq = allPep[accession]
    if decoy: aaSeq = aaSeq[::-1] if aaSeq[0] != 'M' else aaSeq[0]+aaSeq[:0:-1] #reverse the seq; leave M at beginning if there, as specified in Morpheus paper
    if sapAAPos <= 34 or aaSeq[sapAAPos-34] in ['K','R']: trypBeg = True #Either starts at beg of protein or beg is tryptic
    else: trypBeg = False
  return trypBeg,trypEnd

USAGE = "python SummarizeMorpheusResults.py outFile pep.all.fasta listOfFolders"
if len(sys.argv) < 4: print USAGE
else:
  folderList = []
  allPep = read_fasta()
  for i in range (3, len(sys.argv)):
    folderList.append(sys.argv[i])

  try:
    print "opening outFile " + sys.argv[1]
    outF = open(sys.argv[1], 'w')
    outF.write("\t\t\tsummaryaggregate\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcountseqvar\t\t\t\t\t\t\t\tavgpos\t\t\t\t\t\t\t\t\t\tseqvarfdr\t\t\t\t\n")
    outF.write("folder\tline\tsapD\tDataset\tProteins\tMS/MS Spectra\tPSM Morpheus Score Threshold\tTarget PSMs\tDecoy PSMs\tPSM FDR (%)\tUnique Peptide Morpheus Score Threshold\tUnique Target Peptides\tUnique Decoy Peptides\tUnique Peptide FDR (%)\tProtein Group Summed Morpheus Score Threshold\tTarget Protein Groups\tDecoy Protein Groups\tProtein Group FDR (%)\tTotalPsmsNoRaggedEnds\tTotalUniqPepsNoRaggedEnds\tUniProtPsms\tUniProtUniqPeps\tSapPsms\tSapUniqPeps\tJuncPsms\tJuncUniqPeps\taveragePosAll\tmedianPosAll\taveragePosSap\tmedianPosSap\taveragePosJunc\tmedianPosJunc\taverageAll/averageSap\tmedianAll/medianSap\taverageAll/averageJunc\tmedianAll/medianJunc\tTotalFdr\tUniProtFdr\tSapFdr\tJuncFdr\n")
    for folder in folderList:
      allTarget,allDecoy,allFdr,uniTarget,uniDecoy,sapTarget,sapDecoy,juncTarget,juncDecoy = 0,0,0,0,0,0,0,0,0 #seqvarfdr
      totalPsms,totalPeps,uniPsms,uniPeps,sapPsms,sapPeps,juncPsms,juncPeps = 0,0,0,0,0,0,0,0 #countseqvars
      juncPos,sapPos,allPos = [],[],[] #getavgpos
      sumAggregate = ""
      lineCt = 0
      invisSap = ['I','L'] #the leucine - isoleucine transition is invisible to mass spec

      #Print the folder and room for line and cutoff info
      outF.write(folder + "\t\t\t")
      
      #Store the aggregate line from the summary file
      print "opening " + folder + "/summary.tsv"
      with open(folder + "/summary.tsv", 'r') as summary:
        for line in summary:
          if line[:9] == 'AGGREGATE': 
            sumAggregate = line.rstrip() #get rid of tab/newlines at end
            break
          else: continue
        if not sumAggregate: sumAggregate = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
        
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
                  allPos.append(int(line[28]))
                  if int(line[14]) <= 34 and int(line[15]) >= 34 and not (refAA in invisSap and altAA in invisSap): #must contain sap position, and must not have L<->I transition
                      sapPsms += 1
                      sapPos.append(int(line[28]))
                      sapTarget += 1
              elif line[11].find("UniProt:") != -1: 
                uniPsms += 1
                uniTarget += 1
              elif descrips[0] == 'pep:splice': 
                juncPsms += 1
                juncPos.append(int(line[28]))
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
      outF.write(sumAggregate + '\t' + str(totalPsms) + '\t' + str(totalPeps) + '\t' + str(uniPsms)+'\t'+str(uniPeps)+'\t'+str(sapPsms)+'\t'+str(sapPeps)+'\t'+str(juncPsms)+'\t'+str(juncPeps)+'\t')
      #avgpos write
      if len(allPos) == 0: outF.write('N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\t') #avg median all
      else:
        averageSap,medianSap,averageJunc,medianJunc = 0,0,0,0
        averageAll = float(sum(allPos)) / float(len(allPos))
        medianAll = allPos[ len(allPos)/2 ] if len(allPos) % 2 != 0 else (allPos[ len(allPos)/2 ] + allPos[ len(allPos)/2 - 1 ]) / 2
        outF.write(str(averageAll) + '\t' + str(medianAll) + '\t')
        
        if len(sapPos) == 0: outF.write('N/A\tN/A\t') #avg median sap
        else: 
          averageSap = float(sum(sapPos)) / float(len(sapPos)) if sapPos else 'N/A'
          medianSap = sapPos[len(sapPos)/2] if len(sapPos) % 2 != 0 else (sapPos[ len(sapPos)/2 ] + sapPos[ len(sapPos)/2 - 1 ]) / 2
          outF.write(str(averageSap) + '\t' + str(medianSap) + '\t') 
        if len(juncPos) == 0: outF.write('N/A\tN/A\t') #avg median junc
        else: 
          averageJunc = float(sum(juncPos)) / float(len(juncPos)) if juncPos else 'N/A'
          medianJunc = juncPos[len(juncPos)/2] if len(juncPos) % 2 != 0 else (juncPos[ len(juncPos)/2 ] + juncPos[ len(juncPos)/2 - 1 ]) / 2
          outF.write(str(averageJunc) + '\t' + str(medianJunc) + '\t')
        if len(sapPos) == 0: outF.write('N/A\tN/A\t') #avg median all / median sap
        else: outF.write(str(averageAll/averageSap) + '\t' + str(float(medianAll)/float(medianSap)) + '\t')
        if len(juncPos) == 0: outF.write('N/A\tN/A\t') #avg median all / median junc
        else: outF.write(str(averageAll/averageJunc) + '\t' + str(float(medianAll)/float(medianJunc)) + '\t')  
      
      #seqvar write
      try: allFdr = float(allDecoy)/float(allTarget) * 100
      except ZeroDivisionError: allFdr = 'N/A'
      try: uniFdr = float(uniDecoy)/float(uniTarget) * 100
      except ZeroDivisionError: uniFdr = 'N/A'
      try: sapFdr = float(sapDecoy)/float(sapTarget) * 100
      except ZeroDivisionError: sapFdr = 'N/A'
      try: juncFdr = float(juncDecoy)/float(juncTarget) * 100
      except ZeroDivisionError: juncFdr = 'N/A'
      outF.write(str(allFdr)+'\t'+str(uniFdr)+'\t'+str(sapFdr)+'\t'+str(juncFdr)+'\n')

  #except(IndexError): print "error processing a line in PSMs file"
  except(IOError): print "Input/output error"

  
    

  

