def main(options):

    ensure_exon_overlap(options.prefix)
    define_PAS(options.input, options.prefix)

def ensure_exon_overlap(filename):
    UTR_reads = {}
    with open(filename + ".noMP.bed.UTR.EXONS", 'r') as f:
        for line in f:
            line = line.rstrip()
            NAME = 0  
            chrom, read_s, read_e, SRR, score, strand, chrom_2, UTR_s, UTR_e, UTR_ID, zero, strand, chrom_3, exon_s, exon_e, exon_ID, zero, strand = line.split("\t")
            if (chrom, read_s, read_e, SRR, score, strand) not in UTR_reads:
                UTR_reads[(chrom, read_s, read_e, SRR, score, strand)] = 0 #indicates not a pass
                SRR_list = SRR.split(".")
                if "SRR" in SRR: 
                    NAME = 1 
                if len(SRR_list) == 2:
                    if  NAME == 0:
                        orientation = SRR_list[1]
            
                elif len(SRR_list) == 4:
                    if NAME == 1:
                        orientation = SRR_list[3]
            
                else:
                    orientation = SRR_list[2]

                if orientation == "upstream": #meaning it is the negative or reverse
                    if int(read_s) < int(UTR_e):
                        if int(read_e) >= int(exon_s):
                            UTR_reads[(chrom, read_s, read_e, SRR, score, strand)] = 1
                            with open(filename + ".noMP.restricted.bed", 'a') as output:
                                output.write(chrom + "\t" + read_s + "\t" + read_e + "\t" + SRR + "\t" + score + "\t" + strand + "\n") 
                else:
                    if int(read_e) > int(UTR_s):
                        if int(read_s) <= int(exon_e):
                            UTR_reads[(chrom, read_s, read_e, SRR, score, strand)] = 1    
                            with open(filename + ".noMP.restricted.bed", 'a') as output:
                                output.write(chrom+ "\t" + read_s + "\t" + read_e + "\t" + SRR + "\t" + score + "\t" + strand + "\n") 
            else:
                if UTR_reads[(chrom, read_s, read_e, SRR, score, strand)] == 0:
                    SRR_list = SRR.split(".")
                    if len(SRR_list) == 2:
                        orientation = SRR_list[1]
                    else:
                        orientation = SRR_list[2] 
            
                    if orientation == "upstream": #meaning it is the negative or reverse
                        if int(read_s) < int(UTR_e):
                            if int(read_e) >= int(exon_s):
                                UTR_reads[(chrom, read_s, read_e, SRR, score, strand)] = 1
                                with open(filename + ".noMP.restricted.bed", 'a') as output:
                                    utput.write(chrom+ "\t" + read_s+ "\t" + read_e+ "\t" + SRR+ "\t" + score+ "\t" + strand + "\n")
                        else:
                            if int(read_e) > int(UTR_s):
                                if int(read_s) <= int(exon_e):
                                    UTR_reads[(chrom, read_s, read_e, SRR, score, strand)] = 1
                                    with open(filename + ".noMP.restricted.bed", 'a') as output:
                                        output.write(chrom+ "\t" +read_s+ "\t" + read_e+ "\t" + SRR+ "\t" + score+ "\t" + strand + "\n")
                                        
    with open(filename + ".noMP.bed.INTRON", 'r') as f:
        for line in f:
            line = line.rstrip()
            chrom, read_s, read_e, SRR, score, strand = line.split("\t")
            with open(filename + ".noMP.restricted.bed", 'a') as output:
                output.write(chrom+ "\t" +read_s+ "\t" + read_e+ "\t" + SRR+ "\t" + score+ "\t" + strand + "\n")


def define_PAS(input_, filename):
    softclipped = {} 
    with open(input_, 'r') as f:
        for line in f:
            line = line.rstrip()
            name, sc, sc_seq, length = line.split("\t")
            if name not in softclipped:
                softclipped[name] = []
                softclipped[name].append((sc, sc_seq, length))
            else:
                print "repeat"
                softclipped[name].append((sc, sc_seq, length))

    with open(filename + ".noMP.restricted.bed", 'r') as f:
        for line in f:
            line = line.rstrip()
            chrom, s, e, name, score, strand = line.split("\t")
            PAS = 0
            if "SRR" not in name:
                if len(name.split(".")) == 2: 
                    if name.split(".")[1] == "upstream":
                        PAS_S = int(s)-1 
                        PAS_E = int(s)-1 + 100 
                    else:
                        PAS_E = int(e)
                        PAS_S = int(e) - 100
                else: 
                    if name.split(".")[2] == "upstream":
                        PAS_S = int(s)-1
                        PAS_E = int(s)-1 + 100
                    else:
                        PAS_E = int(e)
                        PAS_S = int(e) - 100
            else:
                if len(name.split(".")) == 3:
                    if name.split(".")[2] == "upstream":
                        PAS_S = int(s)-1
                        PAS_E = int(s)-1 + 100
                    else:
                        PAS_E = int(e)
                        PAS_S = int(e) - 100
                else: 
                    if name.split(".")[3] == "upstream":
                        PAS_S = int(s)-1
                        PAS_E = int(s)-1 + 100
                    else:
                        PAS_E = int(e)
                        PAS_S = int(e) - 100            

            if name in softclipped:
                ID = name
                with open(filename + ".noMP.restricted.peaks.bed" ,'a') as output:
                    output.write(chrom + "\t" + str(PAS_S) + "\t" + str(PAS_E) + "\t" + ID + "\t" + score + "\t" + strand + "\n")

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-p", "--prefix", dest="prefix", default = 'test',
        help="input and output prefix (default test)")

    parser.add_option("-i", "--input", dest="input",
        help="input metadata file")


    (options, args) = parser.parse_args()
    
    main(options)
        
