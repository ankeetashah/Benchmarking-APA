import sys 
def main(options):

	check_strand(options.prefix)
	score(options.prefix)

def check_strand(prefix):
	with open(prefix + ".noMP.restricted.peaks.sort.refined.bed", 'r') as f:
		for line in f:
			line = line.rstrip()
			chrom, s, e, name, score, strand = line.split("\t")
			if "," in strand:
				former, latter = strand.split(",")
				with open(prefix + ".noMP.restricted3.peaks.sort.refined.strand.bed", 'a') as output:
					output.write(chrom + "\t" + s + "\t" + e + "\t" + name + "\t" + score + "\t" + former + "\n")
					output.write(chrom + "\t" + s + "\t" + e + "\t" + name + "\t" + score + "\t" + latter + "\n")
			else:
				with open(prefix + ".noMP.restricted3.peaks.sort.refined.strand.bed", 'a') as output:
					output.write(line + "\n")


def score(prefix):
	with open(prefix + ".noMP.restricted3.peaks.sort.refined.strand.bed", 'r') as f:
		for line in f:
			orientation = ""
			line = line.rstrip()
			chrom, s, e, SRR, score, strand = line.split("\t")
			score = len(SRR.split(","))
			SRR_first_list = SRR.split(",")[0].split(".")	
			if len(SRR_first_list) == 2:
				orientation = SRR_first_list[1] 
			elif len(SRR_first_list) == 4:
				orientation = SRR_first_list[3]
			else: 
				orientation = SRR_first_list[2]
				
			if orientation == "upstream":
				strand = "-"
			else:
				strand = "+" 

			with open(prefix + "..noMP.restricted.peaks.sort.refined.score.bed", 'a') as output:
				output.write(chrom + "\t" + s + "\t" + e + "\t" + SRR + "\t" + str(score) + "\t" + strand + "\n") 


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-p", "--prefix", dest="prefix", default = 'test',
        help="input and output prefix (default test)")

    (options, args) = parser.parse_args()
    
    main(options)