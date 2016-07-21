# utility methods for summarizeSsdbgResults

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