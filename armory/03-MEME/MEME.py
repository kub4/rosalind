#!/usr/bin/env python3
"""
New Motif Discovery
===================

The novel-motif finding tool MEME can be found here:
http://meme.nbcr.net/meme/cgi-bin/meme.cgi

Given: A set of protein strings in FASTA format that
share some motif with minimum length 20.

Return: Regular expression for the best-scoring motif.

Sample Dataset
--------------

>Rosalind_7142
PFTADSMDTSNMAQCRVEDLWWCWIPVHKNPHSFLKTWSPAAGHRGWQFDHNFFVYMMGQ
FYMTKYNHGYAPARRKRFMCQTFFILTFMHFCFRRAHSMVEWCPLTTVSQFDCTPCAIFE
WGFMMEFPCFRKQMHHQSYPPQNGLMNFNMTISWYQMKRQHICHMWAEVGILPVPMPFNM
SYQIWEKGMSMGCENNQKDNEVMIMCWTSDIKKDGPEIWWMYNLPHYLTATRIGLRLALY
>Rosalind_4494
VPHRVNREGFPVLDNTFHEQEHWWKEMHVYLDALCHCPEYLDGEKVYFNLYKQQISCERY
PIDHPSQEIGFGGKQHFTRTEFHTFKADWTWFWCEPTMQAQEIKIFDEQGTSKLRYWADF
QRMCEVPSGGCVGFEDSQYYENQWQREEYQCGRIKSFNKQYEHDLWWCWIPVHKKPHSFL
KTWSPAAGHRGWQFDHNFFSTKCSCIMSNCCQPPQQCGQYLTSVCWCCPEYEYVTKREEM
>Rosalind_3636
ETCYVSQLAYCRGPLLMNDGGYGPLLMNDGGYTISWYQAEEAFPLRWIFMMFWIDGHSCF
NKESPMLVTQHALRGNFWDMDTCFMPNTLNQLPVRIVEFAKELIKKEFCMNWICAPDPMA
GNSQFIHCKNCFHNCFRQVGMDLWWCWIPVHKNPHSFLKTWSPAAGHRGWQFDHNFFQMM
GHQDWGTQTFSCMHWVGWMGWVDCNYDARAHPEFYTIREYADITWYSDTSSNFRGRIGQN

Sample Output
-------------

DLWWCWIPVHK[NK]PHSFLKTWSPAAGHRGWQFDHNFF

"""
import sys
from subprocess import check_output

"""
Because the online web version of MEME is too slow (failed
two attempts because of the timeout), I have downloaded,
compiled and installed a local copy of the MEME package.
This Python script is only a glue that uses the meme command
and then extracts the wanted regular expression from the
output file.
""" 

# the meme command we want to use followed by the rosalind dataset
# -protein  : sequences use protein alphabet
# -minw 20  : minumum motif width = 20 (as requested by rosalind)
# -nostatus : do not print progress reports to terminal  
# -text     : output in text format (do not create a directory)
command = "meme -protein -minw 20 -nostatus -text " + sys.argv[1]

# run command with arguments and return its output as a list
memeout = check_output(command, shell=True).decode().splitlines()

# parse memeout and print the result
for i, line in enumerate(memeout):
  if "regular expression" in line:
    regexline = i + 2
print(memeout[regexline])
