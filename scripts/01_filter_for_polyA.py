def main(options):

    samfile = pysam.AlignmentFile(options.input, "rb")
    polyA = pysam.AlignmentFile(options.output + ".polyA.bam", "wb", template=samfile)
    nopolyA = pysam.AlignmentFile(options.output + ".nopolyA.bam", "wb", template=samfile)
    MP = pysam.AlignmentFile(options.output + ".MP.bam", "wb", template=samfile)
    noMP = pysam.AlignmentFile(options.output + ".noMP.bam", "wb", template=samfile)
    
    if options.fasta == None:
        fa = Faidx('../ref/hg19.fa')
        fa.build('hg19.fa', 'hg19.fa.fai')
        fa.index
    else:
        fa = Faidx(options.fasta)

    for read in samfile.fetch():
        softclipped = read_samfile(read)
        print softclipped
        composition(read, softclipped, options.output, polyA, nopolyA, MP, noMP, fa)

def read_samfile(read):
    softclipped = []
    for i in range(0, len(read.cigartuples)):
        if i == 0:
            if int(read.cigartuples[i][0]) == 4:
                softclipped.append(int(read.cigartuples[i][1]))
            else:
                softclipped.append('NA')
        if i == len(read.cigartuples)-1:
            if int(read.cigartuples[i][0]) == 4:
                softclipped.append(int(read.cigartuples[i][1]))
            else:
                softclipped.append('NA')
    return(softclipped)

def rev_polyA_noMP(softclipped, read, up_down, filename, length, polyA, nopolyA, MP, noMP, fa):
    up_down = "upstream"
	base = "T"
    sequence = "" #from here we know we are rev strand -> meaning base = T
    for i in range(0, softclipped):
        sequence += read.query_sequence[i]
    if softclipped < (length -softclipped):
        quantify = 0
        composition = 0
        count = 0
        check = 0
        if len(sequence) >= 20:
            for k in range(len(sequence)-20, len(sequence)):
                count += 1
                if sequence[k] == base:
                    composition += 1
            if len(sequence) <= 40:
                for k in range(0, len(sequence)-20):
                    count += 1
                    if sequence[k] == base:
                            composition += 1
            else:
                for k in range(len(sequence)-40, len(sequence)-20):
                    count += 1
                    if sequence[k] == base:
                        composition += 1
            check = 0.8
        elif len(sequence) > 0 and len(sequence) < 20:
            for k in range(0, len(sequence)):
                count += 1
                if sequence[k] == base:
                    composition += 1
            check = 0.90
        else:
            return
        if softclipped !=0:
            quantify = float(composition)/float(count)
        else:
            quantify = 0
        if quantify < check:
            if up_down not in read.query_name:
                read.query_name = read.query_name + "." + up_down
            nopolyA.write(read) #this is a non-polyA read  
        else:
            if up_down not in read.query_name:
                read.query_name = read.query_name + "." + up_down
            polyA.write(read)
            GENOMIC = ""
            if int(softclipped) < 10:
                use = 10
            else:
                use = int(softclipped)
            CHROM = read.reference_name
            if int(read.reference_start) - use > 1:
                START = int(read.reference_start) - use
            else:
                return
            END = int(read.reference_end)-1
            GENOMIC = str(fa.fetch(CHROM, START, END)).upper()
            mispriming_count = 0.0
            mispriming_quant = 0.0
            if END-START > use+10:
            	for m in range(use, use+10):
            		mispriming_count += 1
            		if GENOMIC[m] == base:
            			mispriming_quant += 1
            else:
            	for m in range(use, (END-START)):
            		mispriming_count += 1
            		if GENOMIC[m] == base:
            			mispriming_quant += 1
            if mispriming_count != 0:
            	mispriming_quant = mispriming_quant/mispriming_count
            if mispriming_quant < 0.6:
            	mispriming_count = 0.0
            	mispriming_quant = 0.0
            	if use-10 >= 0:
            		for m in range(use-10, use):
            			mispriming_count += 1
            			if GENOMIC[m] == base:
            				mispriming_quant += 1
            	else:
            		for m in range(0, use):
            			mispriming_count += 1
            			if GENOMIC[m] == base:
            				mispriming_quant += 1
            	if mispriming_count != 0:
            		mispriming_quant = mispriming_quant/mispriming_count
            	if mispriming_quant < 0.6:
            		noMP.write(read)
            	else:
            		MP.write(read)
            else:
            	MP.write(read)
 
                         
def fwd_polyA_noMP(softclipped, read, strand, filename, length, polyA, nopolyA, MP, noMP, fa):
    up_down = "downstream":
    base = "A"
    sequence = ""
    for i in range(length-softclipped, length):
    	sequence += read.query_sequence[i]
    if len(sequence) < (length - len(sequence)):
    	quantify = 0
    	composition = 0
    	count = 0
    	check = 0
    	if len(sequence) >= 20:
                if len(sequence) <= 40:
                    for k in range(20, len(sequence)):
                            count += 1
                            if sequence[k] == base:
                                    composition += 1
                else:
                    for k in range(20, 40):
                        count += 1
                        if sequence[k] == base:
                            composition += 1
                for k in range(0, 20):
                    count += 1
                    if sequence[k] == base:
                        composition += 1
                    check = 0.8
        elif len(sequence) > 0 and len(sequence) < 20:
            for k in range(0, len(sequence)):
                count += 1
                if sequence[k] == base:
                    composition += 1
            check = 0.90
        else:
            return
        if softclipped !=0:
            quantify = float(composition)/float(count)
        else:
            quantify = 0
        if quantify < check: #this read does not have a polyA stretch
            if up_down not in read.query_name:
                read.query_name = read.query_name + "." + up_down
            nopolyA.write(read)
        else:
            if up_down not in read.query_name:
                read.query_name = read.query_name + "." + up_down
            polyA.write(read)
            GENOMIC = ""
            if int(softclipped) < 10:
                use = 10
            else:
                use = int(softclipped)
            CHROM = read.reference_name
            if int(read.reference_start) - use > 1:
                START = int(read.reference_start) 
            else:
                return
            END = int(read.reference_end)-1 + use
            GENOMIC = str(fa.fetch(CHROM, START, END)).upper()
            mispriming_count = 0.0
            mispriming_quant = 0.0
            if len(GENOMIC) >  len(GENOMIC)-use+10:
            	for m in range(len(GENOMIC)-use, len(GENOMIC)-use+10):
                    mispriming_count += 1
                    if GENOMIC[m] == base:
                       mispriming_quant += 1
            else:
                for m in range(len(GENOMIC)-use, len(GENOMIC)):
                    mispriming_count += 1
                    if GENOMIC[m] == base:
                       mispriming_quant += 1
            if mispriming_count != 0:
                mispriming_quant = mispriming_quant/mispriming_count
                if mispriming_quant < 0.6: #means no mispriming downstream 
                    mispriming_count = 0.0
                    mispriming_quant = 0.0
                    if len(GENOMIC)-use-10 >= 0:
                        for m in range(len(GENOMIC)-use-10, len(GENOMIC)-use):
                            mispriming_count += 1
                            if GENOMIC[m] == base:
                                mispriming_quant += 1
                    else:
                        for m in range(0, len(GENOMIC)-use):
                            mispriming_count += 1
                            if GENOMIC[m] == base:
                                mispriming_quant += 1
                    if mispriming_count != 0:
                        mispriming_quant = mispriming_quant/mispriming_count
                    if mispriming_quant < 0.6: #means no mispriming upstream 
                        noMP.write(read)
                    else:
                        MP.write(read)
                else:
                    MP.write(read)


def composition(read, softclipped, filename, polyA, nopolyA, MP, noMP, fa):
        compositions = []
        composition = 0
        length = read.query_length
        if  length > 0:
            for j in range(0, len(softclipped)):
                if j == 0 and softclipped[j] != 'NA': #reverse strand
                    rev_polyA_noMP(softclipped[j], read, "upstream", filename,length, polyA, nopolyA, MP, noMP, fa)
                if j == 1 and softclipped[j] != 'NA':
                    fwd_polyA_noMP(softclipped[j], read, "downstream", filename, length, polyA, nopolyA, MP, noMP, fa)


if __name__ == "__main__":
    import subprocess
    import sys
    import pysam
    from pyfaidx import Faidx
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-i", "--inputbam", dest="input",
        help="aligned Iso-Seq reads (pooled if worked with multiple libraries)")
    
    parser.add_option("-o", "--outprefix", dest="output", default = 'test',
        help="output prefix (default test)")

    parser.add_option("-f", "--fasta", dest="fasta", default = '../ref/hg19.fa',
        help="reference genome fasta")

    (options, args) = parser.parse_args()
    
    main(options)


