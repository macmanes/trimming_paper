#!/usr/bin/env python

# see http://seqanswers.com/forums/showthread.php?t=16505

import random
import sys

random.seed()
numSeq = int(sys.argv[1])
fileName1 = sys.argv[2]
fileName2 = sys.argv[3]
increment = 0

#if it's a fasta file
if (fileName1.find(".fasta") != -1):
  increment = 2
#else if it's a fastq file
elif (fileName1.find(".fastq") != -1):
  increment = 4
#quit if neither
else:
  sys.stdout.write("This is not a fastA/fastQ file\n")
  sys.exit()

FILE1 = open(fileName1, 'r')
totalReads1 = list()
index = 0
for line in FILE1:
  if(index % increment == 0):
    totalReads1.append(index/increment)
  index += 1
FILE1.close()
if(len(totalReads1) < numSeq):
  sys.stdout.write("You only have "+str(len(totalReads))+" reads!\n")
  sys.exit()
  
if (fileName2):
  FILE2 = open(fileName2, 'r')
  totalReads2 = list()
  index = 0
  for line in FILE2:
    if (index % increment == 0):
      totalReads2.append(index/increment)
    index += 1
  FILE2.close()
  if (len(totalReads1) != len(totalReads2)):
    sys.stdout.write("read counts do not match\n")
    sys.exit()

ttl = len(totalReads1)
random.shuffle(totalReads1)
totalReads1 = totalReads1[0: numSeq]
totalReads1.sort()

FILE1 = open(fileName1, 'r')
if (fileName2):
  FILE2 = open(fileName2, 'r')
listIndex = 0

if (increment == 4):
  OUTFILE1 = open('subsamp_1.fastq', 'w')
  if (fileName2):
    OUTFILE2 = open('subsamp_2.fastq', 'w')
else:
  OUTFILE1 = open('subsamp_1.fasta', 'w')
  if (fileName2):
    OUTFILE2 = open('subsamp_2.fasta', 'w')

for i in range(0, ttl):
  curRead1 = ""
  curRead2 = ""
  for j in range(0, increment):
    curRead1 += FILE1.readline()
    if (fileName2):
      curRead2 += FILE2.readline()
  if (i == totalReads1[listIndex]):
    OUTFILE1.write(curRead1)
    if (fileName2):
      OUTFILE2.write(curRead2)
    listIndex += 1
    if(listIndex == numSeq):
      break
FILE1.close()
if (fileName2):
  FILE2.close()
OUTFILE1.close()
if (fileName2):
  OUTFILE2.close()
