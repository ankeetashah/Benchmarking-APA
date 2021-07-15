import sys

exon_dict = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		line = line.rstrip()
		exon_dict[line] = 0
		

with open(sys.argv[1], 'r') as f:
	for line in f:
		line = line.rstrip()
		if line not in exon_dict:
			with open(sys.argv[3], 'a') as output:
				output.write(line + "\n")




