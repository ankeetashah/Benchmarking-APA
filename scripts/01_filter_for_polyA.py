def main(options):

    samfile = pysam.AlignmentFile(options.input, "rb")
    polyA = pysam.AlignmentFile(options.ouput + ".polyA.bam", "wb", template=samfile)
    nopolyA = pysam.AlignmentFile(options.output + ".nopolyA.bam", "wb", template=samfile)
    MP = pysam.AlignmentFile(options.output + ".MP.bam", "wb", template=samfile)
    noMP = pysam.AlignmentFile(options.ouput + ".noMP.bam", "wb", template=samfile)
    
    if options.fasta == None:
        fa = Faidx('../ref/hg19.fa')
        fa.build('hg19.fa', 'hg19.fa.fai')
        fa.index
    else:
        fa = Faidx(options.fasta)

    softclipped = read_samfile(samfile)
    composition = composition(samfile, softclipped, options.output)

def read_samfile_softclipped(samfile):
    for read in samfile.fetch():
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

def polyA_noMP(softclipped, read, strand, filename):
    base = "A"
    up_down = "upstream"
    if strand == "reverse":
        base = "T"
        up_down = "downstream"

    sequence = ""
    for i in range(0, softclipped):
        sequence += read.query_sequence[i]
    if softclipped[j] < (length -softclipped):
        quantify = 0
        composition = 0
        count = 0
        check = 0
        if len(sequence) >= 20:
            if base == "T":
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

        else:
                if len(sequence) <= 40:
                    for k in range(20, len(sequence)):
                            count += 1
                            if sequence[k] == 'A':
                                    composition += 1
                else:
                    for k in range(20, 40):
                        count += 1
                        if sequence[k] == 'A':
                            composition += 1
                for k in range(0, 20):
                    count += 1
                    if sequence[k] == 'A':
                        composition += 1
                    check = 0.8
        if len(sequence) > 0 and len(sequence) < 20:
            for k in range(0, len(sequence)):
                count += 1
                if sequence[k] == base:
                    composition += 1
            check = 0.90
        else:
            return
        if softclipped[j] !=0:
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
                    #with open(sys.argv[6], 'a') as soft: #metadata for read  -- all reads that have POLYA
                        #soft.write(read.query_name + "\t" + str(len(sequence)) + "\t" + sequence + "\t" + str(length)  + "\n")
            GENOMIC = ""
            if int(softclipped[j]) < 10:
                se = 10
            else:
                use = int(softclipped)
            CHROM = read.reference_name
            if int(read.reference_start) - use > 1:
                START = int(read.reference_start) - use #because softclipping happens at front end
            else:
                return
            END = int(read.reference_end)-1
            if strand == "forward":
                END = int(read.reference_end)-1 + use
            GENOMIC = str(fa.fetch(CHROM, START, END)).upper()
            mispriming_count = 0.0
            mispriming_quant = 0.0
            if strand == "reverse":
                if END-START > use+10: #looking downstream 
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
                    if mispriming_quant < 0.6: #means no mispriming downstream 
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
                        if mispriming_quant < 0.6: #means no mispriming upstream 
                            noMP.write(read)
                            with open(filename + ".noMP.meta.txt", 'a') as soft:
                                soft.write(read.query_name + "\t" + str(len(sequence)) + "\t" + sequence + "\t" + str(length)  + "\n")
                        else:
                            MP.write(read)
                    else:
                        MP.write(read)
            else:
                    if len(GENOMIC) >  len(GENOMIC)-use+10:
                        for m in range(len(GENOMIC)-use, len(GENOMIC)-use+10):
                            mispriming_count += 1
                            if GENOMIC[m] == 'A':
                                mispriming_quant += 1
                    else:
                        for m in range(len(GENOMIC)-use, len(GENOMIC)):
                            mispriming_count += 1                                      
                            if GENOMIC[m] == 'A':
                                mispriming_quant += 1
            if mispriming_count != 0:
                mispriming_quant = mispriming_quant/mispriming_count
                if mispriming_quant < 0.6: #means no mispriming downstream 
                    mispriming_count = 0.0
                    mispriming_quant = 0.0
                    if len(GENOMIC)-use-10 >= 0:
                        for m in range(len(GENOMIC)-use-10, len(GENOMIC)-use):
                            mispriming_count += 1
                            if GENOMIC[m] == 'A':
                                mispriming_quant += 1
                    else:
                        for m in range(0, len(GENOMIC)-use):
                            mispriming_count += 1
                            if GENOMIC[m] == 'A':
                                mispriming_quant += 1
                    if mispriming_count != 0:
                        mispriming_quant = mispriming_quant/mispriming_count
                    if mispriming_quant < 0.6: #means no mispriming upstream 
                        noMP.write(read)
                        with open(filename + ".bed", 'a') as soft:
                            soft.write(read.query_name + "\t" + str(len(sequence)) + "\t" + sequence + "\t" + str(length)  + "\n")
                    else:
                        MP.write(read)
                else:
                        MP.write(read)


def composition(samfile, softclipped, filename):
    for read in samfile.fetch():
        compositions = []
        composition = 0
        length = read.query_length
        if  length > 0:
            for j in range(0, len(softclipped)):
                if j == 0 and softclipped[j] != 'NA': #reverse strand
                    polyA_noMP(softclipped, read, "forward", filename)
                if j == 1 and softclipped[j] != 'NA':
                    polyA_noMP(softclipped, read, "reverse", filename)


if __name__ == "__main__":
    import subprocess
    import sys
    import pysam
    from pyfaidx import Faidx
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-ib", "--inputbam", dest="input",
        help="aligned Iso-Seq reads (pooled if worked with multiple libraries)")
    
    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'test',
        help="output prefix (default test)")

    parser.add_option("-f", "--fasta", dest="fasta", default = '../ref/hg19.fa',
        help="reference genome fasta")

    (options, args) = parser.parse_args()
    
    main(options)


