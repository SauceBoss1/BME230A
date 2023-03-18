from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
import argparse as ag
import sys
import os

import geneFinder as gf
import nGenome_Alignment as nga
import SequenceComparison as sc

class commandLine:
  def __init__(self, inOpts = None):
    self.parser = ag.ArgumentParser(description='The main pipeline for comparing the similarities between genes')

    self.parser.add_argument('-i', '--input', nargs='?', type=str, default=sys.stdin, action='store', help='input file')
    self.parser.add_argument('-o', '--output', nargs='?', type=str, default=sys.stdout, action='store', help='alignment output file')
    self.parser.add_argument('-l', '--minlen', nargs='?', type=int, default=100, action='store', help='the minimum length for each ORF')
    self.parser.add_argument('-p', '--plot', action='store_true', help='plot or just align')
    self.parser.add_argument('-a', '--align', action='store_true', help='enables seq alignment')
    self.parser.add_argument('-s', '--stringLim', action='store', type=int, default=500, help='output string limit')

    self.parser.add_argument('-v','--version', action='version', version='%(prog)s 0.4')

    if inOpts is None:
      self.args = self.parser.parse_args()
    else:
      self.args = self.parser.parse_args(inOpts)

def main(inOpts = None):
  cL = commandLine(inOpts=inOpts)

  if cL.args.plot:
    ss = sc.SequenceSimilarity(cL.args.input)
    ss.generate_scatter_plot()
    ss.generate_violin_plots()
    ss.generate_histogram_plot()
    ss.generate_orf_similarity_plot()
    ss.generate_orf_similarity_plot_with_dendogram()
    #ss.generate_orf_similarity_histogram()
  else:
    out = sys.stdout
    if cL.args.output is not sys.stdout:
      out = open(cL.args.output, 'w')

    records = list(SeqIO.parse(cL.args.input, 'fasta'))
    headers = [r.id for r in records]
    sequences = [str(r.seq) for r in records]
    combineList = list(zip(headers, sequences))
    maxlen = max(len(record.seq) for record in records)

    if cL.args.align:

      for record in records:
        if len(record.seq) != maxlen:
          sequence = str(record.seq).ljust(maxlen,'.')
          record.seq = Seq.Seq(sequence)
      assert all(len(record.seq) == maxlen for record in records)


      temp_out = '{}_padded.fa'.format(os.path.splitext(cL.args.input)[0])
      with open(temp_out,'w') as f:
        SeqIO.write(records, f, 'fasta')

      alignment = AlignIO.read(temp_out, "fasta")
      print(alignment,file=sys.stderr)

      for align in alignment:
        print(align.seq, file=out)
        #print(align.get_alignment_length(), file=sys.stderr)
        # for record in align:
        #   print(record.seq, file=out)



    # alignSeqs = nga.alignTest(sequences, headers, pipeLine=False)
    # alignSeqs.returnAlignments()

    findorf = gf.AlignAllOrfs(combineList, minLength=cL.args.minlen)
    x = findorf.printAlignedOrfs(outFile=out, ultraAlign=True, strLim=cL.args.stringLim)

    #print(findorf.getOrfs, file=sys.stderr)

    out.close()

if __name__ == '__main__':
  main()