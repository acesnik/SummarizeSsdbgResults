#! /usr/bin/python

# Combines GetAveragePositions.py, SeqVarFdrs.py, and CountSeqVars.py to give
# a comprehensive summary of combined seq variant and GPTM approaches.
#
# SAV peptides: to evaluate whether they are tryptic if the first amino acid
# of the SAV peptide is used, the program looks at the pep.all.fasta file to
# reference the preceding amino acid.

import sys
import utility
import os.path

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
            out_file.write("folder\tline\tsapD\tTotalPsmsNoRaggedEnds\tTotalUniqPepsNoRaggedEnds\tUniProtPsms\t"
                           "UniProtUniqPeps\tSapPsms\tSapUniqPeps\tJuncPsms\tJuncUniqPeps\tsav_and_UniProt_pep\t"
                           "junc_and_UniProt_pep\tsav_and_UniProt_psms\tjunc_and_UniProt_psms\tTotalFdr\tUniProtFdr"
                           "\tSapFdr\tJuncFdr\tSavUniFDR\tJuncUniFDR\n")
            for folder in folder_list:
                allTarget, allDecoy, allFdr, uniTarget,uniDecoy, savTarget, sapDecoy, juncTarget, juncDecoy = 0, 0, 0, 0, 0, 0, 0, 0, 0 #seqvarfdr
                sav_and_UniProt_pep, junc_and_UniProt_pep, sav_and_UniProt_decoy, junc_and_UniProt_decoy, sav_and_UniProt_target, junc_and_UniProt_target = 0, 0, 0, 0, 0, 0 # new counters
                totalPsms, totalPeps,uniPsms,uniPeps,savPsms,savPeps,juncPsms,juncPeps = 0,0,0,0,0,0,0,0 #countseqvars
                sav_and_UniProt_psms, junc_and_UniProt_psms = 0, 0
                invisSap = ['I', 'L'] #the leucine - isoleucine transition is invisible to mass spec

                # Print the folder and room for line and cutoff info
                out_file.write(folder + "\t\t\t")

                unique_peptide_file = folder + "/unique_peptides.tsv"
                if not os.path.isfile(unique_peptide_file): unique_peptide_file = folder + "/aggregate.unique_peptides.tsv"
                if not os.path.isfile(unique_peptide_file): print "unique peptides file in " + folder + " is not valid."; exit(2);

                psms_file = folder + "/PSMs.tsv"
                if not os.path.isfile(psms_file): psms_file = folder + "/aggregate.PSMs.tsv"
                if not os.path.isfile(psms_file): print "psms file in " + folder + " is not valid."; exit(2);

                # Count sequence variants in the unique peptides folder
                print "opening " + unique_peptide_file
                with open(unique_peptide_file, 'r') as uniqPeps:
                    for line in uniqPeps:
                        line = line.split('\t')
                        if line[0][:8] == 'Filename': continue
                        elif float(line[30]) < 1 and line[26] == 'True':
                            totalPeps += 1
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
                                    contains_sav_position = int(line[14]) <= 34 and int(line[15]) >= 34
                                    is_invisible_transition = reference_aa in invisSap and alternate_aa in invisSap
                                    if contains_sav_position and not is_invisible_transition:
                                        savPeps += 1
                                        if line[11].find("UniProt:") != -1: sav_and_UniProt_pep += 1
                                else:
                                    totalPeps -= 1
                                    continue
                            if description[0] == 'pep:splice':
                                juncPeps += 1
                                if line[11].find("UniProt:") != -1: junc_and_UniProt_pep += 1
                            if line[11].find("UniProt:") != -1: uniPeps += 1

                # Work with PSMs folder for each metric
                print "opening " + psms_file
                with open(psms_file, "r") as psms:  # May need to change the name.
                    for line in psms:
                        line = line.split('\t')
                        if line[0][:8] == 'Filename': continue

                        description = line[13].split('|')
                        description = description[2].split(' ')

                        if not float(line[30]) < 1: break
                        if line[26] == 'True':
                            # countseqvars
                            allTarget += 1
                            totalPsms += 1
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

                                    contains_sav_position = int(line[14]) <= 34 and int(line[15]) >= 34
                                    is_invisible_transition = reference_aa in invisSap and alternate_aa in invisSap
                                    if contains_sav_position and not is_invisible_transition:
                                        savPsms += 1
                                        savTarget += 1
                                        if line[11].find("UniProt:") != -1:
                                            sav_and_UniProt_psms += 1
                                            sav_and_UniProt_target += 1
                                else:
                                    totalPsms -= 1
                                    continue
                            if description[0] == 'pep:splice':
                                juncPsms += 1
                                juncTarget += 1
                                if line[11].find("UniProt:") != -1:
                                    junc_and_UniProt_psms += 1
                                    junc_and_UniProt_target += 1
                            if line[11].find("UniProt:") != -1:
                                uniPsms += 1
                                uniTarget += 1

                        else:
                            allDecoy += 1
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
                                    contains_sav_position = int(line[14]) <= 34 and int(line[15]) >= 34
                                    is_invisible_transition = reference_aa in invisSap and alternate_aa in invisSap
                                    if contains_sav_position and not is_invisible_transition:
                                        sapDecoy += 1 # must contain sap position
                                        if line[11].find("UniProt:") != -1: sav_and_UniProt_decoy += 1
                                else:
                                    allDecoy -= 1
                                    continue
                            if description[0] == 'pep:splice':
                                juncDecoy += 1
                                if line[11].find("UniProt:") != -1: junc_and_UniProt_decoy += 1
                            if line[11].find("UniProt:") != -1: uniDecoy += 1

                out_file.write(str(totalPsms) + '\t' + str(totalPeps) + '\t' + str(uniPsms) + '\t' + str(uniPeps) +
                               '\t' + str(savPsms) + '\t' + str(savPeps) + '\t' + str(juncPsms) + '\t' + str(juncPeps) +
                               '\t' + str(sav_and_UniProt_pep) + '\t' + str(junc_and_UniProt_pep) + '\t' +
                               str(sav_and_UniProt_psms) + '\t'+str(junc_and_UniProt_psms) + '\t')

                # seqvar write
                try: allFdr = float(allDecoy)/float(allTarget) * 100
                except ZeroDivisionError: allFdr = 'N/A'
                try: uniFdr = float(uniDecoy)/float(uniTarget) * 100
                except ZeroDivisionError: uniFdr = 'N/A'
                try: sapFdr = float(sapDecoy)/float(savTarget) * 100
                except ZeroDivisionError: sapFdr = 'N/A'
                try: juncFdr = float(juncDecoy)/float(juncTarget) * 100
                except ZeroDivisionError: juncFdr = 'N/A'

                try: savUniFDR = float(sav_and_UniProt_decoy) / float(sav_and_UniProt_target) * 100
                except ZeroDivisionError: savUniFDR = 'N/A'
                try: juncUniFDR = float(junc_and_UniProt_decoy) / float(junc_and_UniProt_target) * 100
                except ZeroDivisionError: juncUniFDR = 'N/A'

                out_file.write(str(allFdr) + '\t' + str(uniFdr) + '\t' + str(sapFdr) + '\t' + str(juncFdr) + '\t' +
                               str(savUniFDR) + '\t' + str(juncUniFDR) + '\n')

        # except(IndexError): print "error processing a line in PSMs file"
        except IOError: print "Input/output error"

if __name__ == "__main__": __main__()
  

