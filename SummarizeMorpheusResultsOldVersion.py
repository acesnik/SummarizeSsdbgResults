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
def sap_eval(descriptions, line, allPep):
    accession = line[13].split('|')[1]
    decoy = accession[:6] == 'DECOY_'
    if decoy: accession = accession[6:]

    if len(descriptions) == 16: sap = descriptions[7].split(':') # snpefftopeptides
    elif len(descriptions) == 12: sap = descriptions[5].split(':') # ProteogenomicDBGenerator
    else:
        print "Error: SAV description not recognized\n" + str(line)
        exit(2)

    # Evaluate arbitrary ends of peptide for being tryptic
    first_aa_index = 1
    end_aa_index = 67
    sav_amino_acid_position = int(sap[1][1:-1])
    end_is_tryptic = int(line[15]) != end_aa_index or line[12][-1] in ['K', 'R']

    aaSeq = allPep[accession]
    if decoy: aaSeq = aaSeq[::-1] if aaSeq[0] != 'M' else aaSeq[0] + aaSeq[:0:-1]  # reverse the seq; leave M at beginning if there, as specified in Morpheus paper
    if not allPep[accession]:
        print "Did not find the reference sequence for " + accession + "\nPlease check that the pep.all file is the correct version."
    starts_at_beginning_aa = int(line[14]) == first_aa_index
    contains_start_aa_of_protein = sav_amino_acid_position <= 34
    preceding_amino_acid_is_lys_or_arg = aaSeq[sav_amino_acid_position - 34] in ['K','R']
    beginning_is_tryptic = not starts_at_beginning_aa or contains_start_aa_of_protein or preceding_amino_acid_is_lys_or_arg

    return beginning_is_tryptic, end_is_tryptic


def __main__():
    USAGE = "python SummarizeMorpheusResults.py outFile pep.all.fasta folder1 folder2 ... folderN"
    if len(sys.argv) < 4: print USAGE
    else:
        folder_list = sys.argv[3:]
        fasta = utility.read_fasta()

        try:
            print "opening outFile " + sys.argv[1]
            out_file = open(sys.argv[1], 'w')
            out_file.write("\t\t\tcountseqvar\t\t\t\t\t\t\t\tseqvarfdr\t\t\t\t\n")
            out_file.write("folder\tline\tsapD\tTotalPsmsNoRaggedEnds\tTotalUniqPepsNoRaggedEnds\tUniProtPsms\tUniProtUniqPeps\tSapPsms\tSapUniqPeps\tJuncPsms\tJuncUniqPeps\tTotalFdr\tUniProtFdr\tSapFdr\tJuncFdr\n")
            for folder in folder_list:
                allTarget, allDecoy, allFdr, uniTarget, uniDecoy, savTarget, sapDecoy, juncTarget, juncDecoy = 0, 0, 0, 0, 0, 0, 0, 0, 0 #seqvarfdr
                totalPsms, totalPeps,uniPsms,uniPeps,savPsms,savPeps,juncPsms,juncPeps = 0,0,0,0,0,0,0,0 #countseqvars
                invisSap = ['I', 'L'] #the leucine - isoleucine transition is invisible to mass spec

                # Print the folder and room for line and cutoff info
                out_file.write(folder + "\t\t\t")

                # Count sequence variants in the unique peptides folder
                print "opening " + folder + "/unique_peptides.tsv"
                uniqPeps = open(folder + "/unique_peptides.tsv", 'r')
                with open(folder + "/unique_peptides.tsv", 'r') as uniqPeps:
                    for line in uniqPeps:
                        line = line.split('\t')
                        if line[0][:8] == 'Filename': continue
                        elif float(line[30]) < 1 and line[26] == 'True':
                            description = line[13].split('|')
                            description = description[2].split(' ')
                            if description[0] in ['pep:sap', 'pep:sav']:
                                if len(description) == 16: sap = description[7].split(':') # snpefftopeptides
                                elif len(description) == 12: sap = description[5].split(':') # ProteogenomicDBGenerator
                                else:
                                    print "Error: SAV description not recognized\n" + str(line)
                                    exit(2)
                                reference_aa, alternate_aa = sap[1][0], sap[1][-1]

                                # Evaluate arbitrary ends of peptide for being tryptic, and if it's tryptic, decide whether to count it as a sap peptide
                                beginning_is_tryptic, end_is_tryptic = sap_eval(description, line, fasta)
                                if beginning_is_tryptic and end_is_tryptic:
                                    totalPeps += 1
                                    contains_sav_position = int(line[14]) <= 34 and int(line[15]) >= 34
                                    is_invisible_transition = reference_aa in invisSap and alternate_aa in invisSap
                                    if contains_sav_position and not is_invisible_transition:
                                        savPeps += 1
                            elif line[11].find("UniProt:") != -1: uniPeps += 1
                            elif description[0] == 'pep:splice': juncPeps += 1
                            else: totalPeps += 1

                # Work with PSMs folder for each metric
                print "opening " + folder + "/PSMs.tsv"
                with open(folder + "/PSMs.tsv", "r") as psms:
                    for line in psms:
                        line = line.split('\t')
                        if line[0][:8] == 'Filename': continue

                        description = line[13].split('|')
                        description = description[2].split(' ')

                        if not float(line[30]) < 1: break
                        if line[26] == 'True':
                            # countseqvars
                            if description[0] in ['pep:sap', 'pep:sav']:
                                if len(description) == 16: sap = description[7].split(':') # snpefftopeptides
                                elif len(description) == 12: sap = description[5].split(':') # ProteogenomicDBGenerator
                                else:
                                    print "Error: SAV description not recognized\n" + str(line)
                                    exit(2)
                                reference_aa, alternate_aa = sap[1][0], sap[1][-1]

                                # Evaluate arbitrary ends of peptide for being tryptic, and if it's tryptic, decide whether to count it as a sap peptide
                                beginning_is_tryptic, end_is_tryptic = sap_eval(description, line, fasta)
                                if beginning_is_tryptic and end_is_tryptic:
                                    totalPsms += 1
                                    allTarget += 1
                                    contains_sav_position = int(line[14]) <= 34 and int(line[15]) >= 34
                                    is_invisible_transition = reference_aa in invisSap and alternate_aa in invisSap
                                    if contains_sav_position and not is_invisible_transition:
                                        savPsms += 1
                                        savTarget += 1
                            elif line[11].find("UniProt:") != -1:
                                uniPsms += 1
                                uniTarget += 1
                            elif description[0] == 'pep:splice':
                                juncPsms += 1
                                juncTarget += 1
                            else:
                                totalPsms += 1
                                allTarget += 1
                        else:
                            if description[0] in ['pep:sap', 'pep:sav']:
                                if len(description) == 16: sap = description[7].split(':') # snpefftopeptides
                                elif len(description) == 12: sap = description[5].split(':') # ProteogenomicDBGenerator
                                else:
                                    print "Error: SAV description not recognized\n" + str(line)
                                    exit(2)
                                reference_aa, alternate_aa = sap[1][0], sap[1][-1]

                                # If it's tryptic, decide whether to count it as a sap peptide
                                beginning_is_tryptic, end_is_tryptic = sap_eval(description, line, fasta)
                                if beginning_is_tryptic and end_is_tryptic:
                                    allDecoy += 1
                                    contains_sav_position = int(line[14]) <= 34 and int(line[15]) >= 34
                                    is_invisible_transition = reference_aa in invisSap and alternate_aa in invisSap
                                    if contains_sav_position and not is_invisible_transition:
                                        sapDecoy += 1 # must contain sap position
                            elif line[11].find("UniProt:") != -1: uniDecoy += 1
                            elif description[0] == 'pep:splice': juncDecoy += 1
                            else: totalPeps += 1

                # write to output: sumAggregate, countseqvar
                out_file.write(str(totalPsms) + '\t' + str(totalPeps) + '\t' + str(uniPsms) + '\t' + str(uniPeps) + '\t' + str(savPsms) + '\t' + str(savPeps) + '\t' + str(juncPsms) + '\t' + str(juncPeps) + '\t')

                # seqvar write
                try: allFdr = float(allDecoy)/float(allTarget) * 100
                except ZeroDivisionError: allFdr = 'N/A'
                try: uniFdr = float(uniDecoy)/float(uniTarget) * 100
                except ZeroDivisionError: uniFdr = 'N/A'
                try: sapFdr = float(sapDecoy)/float(savTarget) * 100
                except ZeroDivisionError: sapFdr = 'N/A'
                try: juncFdr = float(juncDecoy)/float(juncTarget) * 100
                except ZeroDivisionError: juncFdr = 'N/A'
                out_file.write(str(allFdr) + '\t' + str(uniFdr) + '\t' + str(sapFdr) + '\t' + str(juncFdr) + '\n')

        # except(IndexError): print "error processing a line in PSMs file"
        except IOError: print "Input/output error"

if __name__ == "__main__": __main__()
  
