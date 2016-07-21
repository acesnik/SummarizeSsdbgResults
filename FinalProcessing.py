#! /usr/bin/python

# Processing multiple PSMs.tsv files to count uniprot mods and psms. Then,
# summarize the unique peptides and protein groups.

__author__ = "anthonycesnik"
__date__ = "$Aug 1, 2015 5:29:41 PM$"

import sys
import utility

#Evaluate arbitrary ends of peptide for being tryptic
#def sav_eval(descrips, line, allPep):
#  accession = line[13].split('|')[1]
#  decoy = accession[:6] == 'DECOY_'
#  if decoy: accession = accession[6:]
#  sav = descrips[7].split(':')
#  savAAPos = int(sav[1][1:-1])
#  aaSeq = ''
#  #Evaluate arbitrary ends of peptide for being tryptic
#  trypBeg, trypEnd = False, int(line[15]) != 67 or line[12][-1] in ['K','R']
#  if int(line[14]) != 1: trypBeg = True
#  elif not allPep[accession]: print "Did not find the reference sequence for " + accession + "\nPlease check that the pep.all file is the correct version."
#  else: 
#    aaSeq = allPep[accession]
#    if decoy: aaSeq = aaSeq[::-1] if aaSeq[0] != 'M' else aaSeq[0]+aaSeq[:0:-1] #reverse the seq; leave M at beginning if there, as specified in Morpheus paper
#    if savAAPos <= 34 or aaSeq[savAAPos-34] in ['K','R']: trypBeg = True #Either starts at beg of protein or beg is tryptic
#    else: trypBeg = False
#  return trypBeg,trypEnd
  
#List the uniprot mods in an PSM
aa_abbrev_dict = { 'Phe':'F', 'Leu':'L', 'Ser':'S', 'Tyr':'Y', 'Cys':'C', 
  'Trp':'W', 'Pro':'P', 'His':'H', 'Gln':'Q', 'Arg':'R', 'Ile':'I', 'Met':'M',
  'Thr':'T', 'Asn':'N', 'Lys':'K', 'Val':'V', 'Ala':'A', 'Asp':'D', 'Glu':'E',
  'Gly':'G' }
def get_mods(pep_seq):
  mods = []
  findresult = pep_seq.find('UniProt: ')
  while findresult > -1:
    modsplit = pep_seq[findresult+9:].split(')')
    modname = modsplit[0]
    while len(modsplit) > 1 and (modsplit[1] == "" or (modsplit[1][0] != "." and modsplit[1][0] not in aa_abbrev_dict.values())):
      modname += ")" + modsplit[1]
      del modsplit[1]
    findresult = pep_seq.find('UniProt: ',findresult+1)
    mods.append(modname)
  return mods

#Add the PSM to the hashes
def add_mod(mod, folder, modhash_f, modhash_t):
  if (mod,folder) in modhash_f:
    modhash_t[mod] += 1
    modhash_f[(mod,folder)] += 1
  elif mod in modhash_t:
    modhash_t[mod] += 1    
    modhash_f[(mod,folder)] = 1
  else:
    modhash_f[(mod,folder)],modhash_t[mod] = 1,1
  
#Evaluate arbitrary ends of peptide for being tryptic
#If it's tryptic, decide whether to count it as a sav peptide
def sav_eval(descrips, line, allPep):
  invisSav = ['L','I']
  if len(descrips) == 16: sav = descrips[7].split(':') #snpefftopeptides
  elif len(descrips) == 12: sav = descrips[5].split(':') #ProteogenomicDBGenerator
  else: 
    print "Error: SAV description not recognized\n" + str(line)
    exit(2)
  refAA, savAAPos, altAA = sav[1][0], int(sav[1][1:-1]), sav[1][-1]
  accession = line[13].split('|')[1]
  decoy = accession[:6] == 'DECOY_'
  if decoy: accession = accession[6:]
  aaSeq = ''
  #Evaluate arbitrary ends of peptide for being tryptic
  trypBeg, trypEnd = False, int(line[15]) != 67 or line[12][-1] in ['K','R']
  if int(line[14]) != 1: trypBeg = True
  elif not allPep[accession]: print "Did not find the reference sequence for " + accession + "\nPlease check that the pep.all file is the correct version."
  else: 
    aaSeq = allPep[accession]
    if decoy: aaSeq = aaSeq[::-1] if aaSeq[0] != 'M' else aaSeq[0]+aaSeq[:0:-1] #reverse the seq; leave M at beginning if there, as specified in Morpheus paper
    if savAAPos <= 34 or aaSeq[savAAPos-34] in ['K','R']: trypBeg = True #Either starts at beg of protein or beg is tryptic
    else: trypBeg = False
  return trypBeg and trypEnd and int(line[14]) <= 34 and int(line[15]) >= 34 and not (refAA in invisSav and altAA in invisSav) #must contain sav position and not be an I<->L transition

# list the sav/nsj info from a sav/nsj peptide description
def info(descrips):
  descrips = descrips.split('|')
  info = [descrips[1]] #accession
  descrips = descrips[2].split(' ')
  if descrips[0] in ['pep:sap','pep:sav']:
    if len(descrips) == 16: #snpefftopeptides
      info.append(descrips[7].split(':')[1]) #sap
      info.append(':'.join(descrips[8].split(':')[-2:])) #snp location
      info.append(descrips[9].split(':')[1]) #codon change
      info.append(':'.join(descrips[13].split(':')[-4:])) #chrom-loci
    elif len(descrips) == 12: #ProteogenomicDBGenerator
      info.append(descrips[5].split(':')[1])
      info.append(':'.join(descrips[3].split(':')[-2:]))
      info.append(descrips[4].split(':')[1])
      info.append(':'.join(descrips[9].split(':')[-4:]))
    else: 
      print "Error: pep:sav description not recognized\n" + str(descrips)
      exit(2)
  elif descrips[0] == 'pep:splice':
    if len(descrips) == 8: #ProteogenomicDBGenerator
      info.append(':'.join(descrips[5].split(':')[2:4])) #chrom-strand
      exonDescrips = descrips[1].split(':')
      exonLoci = ':'.join(descrips[5].split(':')[4:]).split(';')
    elif len(descrips) == 14: #TranslateBedSequences-AC
      info.append(':'.join(descrips[10].split(':')[-2:])) #chrom-strand
      exonDescrips = descrips[6].split(':')
      exonLoci = descrips[11][5:].split(';')
    else:
      print "Error: pep:splice description not recognized\n" + str(descrips)
      exit(2)
    info.append(exonDescrips[1]) #exon1name
    info.append(exonDescrips[2].split(';')[0]) #exon1type
    info.append(exonLoci[0]) #exon1loci
    info.append(exonDescrips[2].split(';')[1]) #exon2name
    info.append(exonDescrips[3]) #exon2type
    info.append(exonLoci[1]) #exon2loci
  elif descrips[0] == 'pep:sep':
    info.append(descrips[2].split(',')[0])
    info.append(descrips[2].split(',')[1])
    info.append(descrips[3])
  return info

USAGE = "python FinalProcessing.py outFolder pep.all.fasta listOfFolders"
if len(sys.argv) < 4: print USAGE
else:
  #Infiles
  folderList = []
  allPep = utility.read_fasta() #fasta file
  for i in range (3, len(sys.argv)):
    folderList.append(sys.argv[i])

  try:
    #Outfiles
    print "opening outFile " + sys.argv[1] + "/uniqpep_summary.txt"
    outF_uniqpeps = open(sys.argv[1] + "/uniqpep_summary.txt", 'w')
    print "opening outFile " + sys.argv[1] + "/uniqpep_sav_summary.txt"
    outF_uniqpeps_sav = open(sys.argv[1] + "/uniqpep_sav_summary.txt", 'w')
    print "opening outFile " + sys.argv[1] + "/uniqpep_nsj_summary.txt"
    outF_uniqpeps_nsj = open(sys.argv[1] + "/uniqpep_nsj_summary.txt", 'w')
    print "opening outFile " + sys.argv[1] + "/uniqpep_sep_summary.txt"
    outF_uniqpeps_sep = open(sys.argv[1] + "/uniqpep_sep_summary.txt", 'w')
    print "opening outFile " + sys.argv[1] + "/uniprot_mod_summary.txt"
    outF_uniprot = open(sys.argv[1] + "/uniprot_mod_summary.txt", 'w')
    print "opening outFile " + sys.argv[1] + "/protein_group_summary.txt"
    outF_protgps = open(sys.argv[1] + "/protein_group_summary.txt", 'w')
    print "opening summaryFile " + sys.argv[1] + "/summary.txt"
    outF_summary = open(sys.argv[1] + "/summary.txt", 'w')
    
    #Process UniProt PSMs
    modhash_folder = {}
    modhash_total = {}
    total_targets = {}
    for folder in folderList:
      print "opening " + folder + "/PSMs.tsv"
      with open(folder + "/PSMs.tsv", "r") as psms: 
        for line in psms:
          line = line.split('\t')
          if line[0][:8] == 'Filename': continue
          isTarget = float(line[30]) < 1 and line[26] == 'True'
          isTargetUniProt = float(line[30]) < 1 and line[26] == 'True' and line[11].find("UniProt:") != -1
          if isTarget:
            if folder in total_targets: total_targets[folder] += 1
            else: total_targets[folder] = 1
          if isTargetUniProt: 
            for mod in get_mods(line[11]): 
              add_mod(mod, folder, modhash_folder, modhash_total)
    #output to file
    outF_uniprot.write("modification\t")
    for folder in folderList: outF_uniprot.write(folder + '\t')
    outF_uniprot.write('Total\n')
    modlist = sorted(modhash_total.keys())
    for mod in modlist:
      outF_uniprot.write(mod + '\t')
      for folder in folderList:
        if (mod, folder) in modhash_folder: outF_uniprot.write(str(float(modhash_folder[(mod, folder)]) / float(total_targets[folder])) + '\t')
        else: outF_uniprot.write('0\t')
      outF_uniprot.write(str(modhash_total[mod]) + '\n')
    outF_uniprot.close()
    
    #Process Unique Peps
    allhash = {}
    savhash = {}
    nsjhash = {}
    sephash = {}
    for folder in folderList:
      print "opening " + folder + "/unique_peptides.tsv"
      with open(folder + "/unique_peptides.tsv", "r") as uniqPeps:
        for line in uniqPeps:
          line = line.split('\t')
          if line[0][:8] == 'Filename': continue
          elif float(line[30]) < 1 and line[26] == 'True':
            allhashkey = (line[13],line[11],line[14],line[15],line[16])
            descrips = line[13].split('|')[2].split(' ')
            if descrips[0] in ['pep:sap','pep:sav'] and sav_eval(descrips, line, allPep): 
              if line[13] in savhash and folder not in savhash[line[13]]: savhash[line[13]].append(folder) # protein description : folder list it was found in
              else: savhash[line[13]] = [folder]
              if allhashkey in allhash and folder not in allhash[allhashkey]: allhash[allhashkey].append(folder)
              else: allhash[allhashkey] = [folder]
            elif descrips[0] == 'pep:splice':
              if line[13] in nsjhash and folder not in nsjhash[line[13]]: nsjhash[line[13]].append(folder)
              else: nsjhash[line[13]] = [folder]
              if allhashkey in allhash and folder not in allhash[allhashkey]: allhash[allhashkey].append(folder)
              else: allhash[allhashkey] = [folder]
            elif descrips[0] == 'pep:sep':
              if line[13] in sephash and folder not in sephash[line[13]]: sephash[line[13]].append(folder)
              else: sephash[line[13]] = [folder]
              if allhashkey in allhash and folder not in allhash[allhashkey]: allhash[allhashkey].append(folder)
              else: allhash[allhashkey] = [folder]
            else:
              if allhashkey in allhash and folder not in allhash[allhashkey]: allhash[allhashkey].append(folder)
              else: allhash[allhashkey] = [folder]
          else: continue
    #output to file
    outF_uniqpeps_sav.write('accession\tsav\tsnv_location\tcodon_change\tchrom-loci\tfound_in_these_searches\n')
    outF_uniqpeps_nsj.write('accession\tchrom-strand\texon1name\texon1type\texon1loci\texon2name\texon2type\texon2loci\tfound_in_these_searches\n')
    outF_uniqpeps_sep.write('accession\tstart_site\tbiotype\tchrom-loci\tfound_in_these_searches\n')
    outF_uniqpeps.write('description\tpeptide_sequence\tstart_residue\tstop_residue\tmissed_cleavages\tfound_in_these_searches\n')
    for sav in savhash: 
      savinfo = info(sav)
      for item in savinfo: outF_uniqpeps_sav.write(item + '\t')
      for folder in savhash[sav]: outF_uniqpeps_sav.write(folder + ',')
      outF_uniqpeps_sav.write('\n')
    for nsj in nsjhash:
      nsjinfo = info(nsj)
      for item in nsjinfo: outF_uniqpeps_nsj.write(item + '\t')
      for folder in nsjhash[nsj]: outF_uniqpeps_nsj.write(folder + ',')
      outF_uniqpeps_nsj.write('\n')
    for sep in sephash:
      sepinfo = info(sep)
      for item in sepinfo: outF_uniqpeps_sep.write(item + '\t')
      for folder in sephash[sep]: outF_uniqpeps_sep.write(folder + ',')
      outF_uniqpeps_sep.write('\n')
    for pep in allhash:
      for info in pep: outF_uniqpeps.write(info+'\t')
      for folder in allhash[pep]: outF_uniqpeps.write(folder + ',')
      outF_uniqpeps.write('\n')
    outF_uniqpeps_sav.close()
    outF_uniqpeps_nsj.close()
    outF_uniqpeps.close()
    
    #Process protein groups
    prothash = {}
    for folder in folderList:
      print "opening " + folder + "/protein_groups.tsv"
      with open(folder + "/protein_groups.tsv", "r") as protGps:
        for line in protGps:
          line = line.split('\t')
          if line[0].strip() == 'Protein Description': continue
          elif float(line[14]) < 1 and line[10] == 'True':
            pepDescrips = line[0].split(';; ')
            pepDescrips[-1] = pepDescrips[-1].strip(';') #remove the hanging semicolon
            for pep in pepDescrips:
              descrips = pep.split('|')[2].split(' ')
              if descrips[0] in ['pep:sap','pep:sav'] and pep in savhash: #the unique peptide must have been found in the uniq pep eval
                if pep in prothash and folder not in prothash[pep]: prothash[pep].append(folder)
                else: prothash[pep] = [folder]
              else: #uniprot or nsj pep
                if pep in prothash and folder not in prothash[pep]: prothash[pep].append(folder)
                else: prothash[pep] = [folder]
          else: continue
    #output to file
    outF_protgps.write('id\tfound_in_these_searches\n')
    for pep in prothash:
      outF_protgps.write(pep + '\t')
      for folder in prothash[pep]: outF_protgps.write(folder + ',')
      outF_protgps.write('\n')
    outF_protgps.close()
    
    #Summary
    folders = ""
    for folder in folderList: folders += folder + ","
    outF_summary.write("field\tsavs\tnsjs\tseps\tptmpeps\tallunmodpeps\tfolders\n")
    savuniqline,nsjuniqline,sepuniqline,moduniqline,unmoduniqline = 0,0,0,0,0
    uniqsavpep,uniqnsjpep,uniqseppep,uniqmodpep,uniqunmodpep = 0,0,0,0,0
    modpeps,unmodpeps,savpeps,nsjpeps,seppeps = [],[],[],[],[]
    for item in allhash:
      #Note: I chose to use the full peptide data for these two,
      #but then only the description for the SAV and NSJ peptides because
      #the SAV and NSJ are captured not by different ends, but then additional
      #information about the PTMs may be found by different missed cleavages
      if item[1].find("UniProt:") != -1: 
        uniqmodpep += 1
        if len(allhash[item]) == 1: moduniqline += 1
        #modpeps.append(item)
      else:
        if (item[0].find("pep:sap") != -1 or item[0].find("pep:sav") != -1) and item[0] not in savpeps:
          uniqsavpep += 1
          if len(allhash[item]) == 1: savuniqline += 1
          savpeps.append(item[0])
        elif item[0].find("pep:splice") != -1 and item[0] not in nsjpeps:
          uniqnsjpep += 1
          if len(allhash[item]) == 1: nsjuniqline += 1
          nsjpeps.append(item[0])
        elif item[0].find("pep:sep") != -1 and item[0] not in seppeps:
          uniqseppep += 1
          if len(allhash[item]) == 1: sepuniqline += 1
          seppeps.append(item[0])
        else:
          uniqunmodpep += 1
          if len(allhash[item]) == 1: unmoduniqline += 1
          #unmodpeps.append(item)
    outF_summary.write("uniq_peps\t" + str(uniqsavpep) + '\t' + str(uniqnsjpep) + '\t' + str(uniqseppep) + '\t' + str(uniqmodpep) + '\t' + str(uniqunmodpep) + '\t' + folders + '\n')
    outF_summary.write("uniq_to_folder\t" + str(savuniqline) + '\t' + str(nsjuniqline) + '\t' + str(sepuniqline) + '\t' + str(moduniqline) + '\t' + str(unmoduniqline) + '\t' + folders + '\n')
    outF_summary.write("percent_uniq\t" + (str(float(savuniqline)/float(uniqsavpep)) if uniqsavpep else 'N/A') + '\t' + (str(float(nsjuniqline)/float(uniqnsjpep)) if uniqnsjpep else 'N/A') + '\t' + (str(float(sepuniqline)/float(uniqseppep)) if uniqseppep else 'N/A') + '\t' + (str(float(moduniqline)/float(uniqmodpep)) if uniqmodpep else 'N/A') + '\t' + (str(float(unmoduniqline)/float(uniqunmodpep)) if uniqunmodpep else 'N/A') + '\t' + folders + '\n')
    outF_summary.close()
    
  except(IOError): print "IOError"