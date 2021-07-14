
def main(options):

    if options.annotate == "annotate":
        annotate(options.reference, options.prefix)
    if options.usage == "usage":
        usage(options.prefix)

def annotate(reference, filename):  
    genes = {}
    with open(reference, 'r') as f:
        for line in f:
            line = line.rstrip()
            ID, name, strand= line.split("\t")  
            if ID not in genes:
                genes[ID] = (name, strand)


    with open(filename + ".homer.annotate.BED", 'r') as f:
        next(f)
        for line in f:
            line = line.rstrip()
            chrom = line.split("\t")[1]
            s = str(int(line.split("\t")[2])-1)
            e = line.split("\t")[3]
            strand = line.split("\t")[4]
            score = line.split("\t")[5]
            name = "NA"
            gene = "NA"
            if  "NA" not in line.split("\t")[7]:
                print line.split("\t")[7]
                if line.split("\t")[7] == "Intergenic":
                    gene = line.split("\t")[10]
                if  "exon" in line.split("\t")[7]  or  "intron" in line.split("\t")[7]:
                    gene = line.split("\t")[7].split("(")[1].split(",")[0] 
                if "TTS" in line.split("\t")[7] or "TSS" in line.split("\t")[7]:
                    gene = line.split("\t")[7].split("(")[1].split(")")[0]  
                if gene in genes:
                    name = genes[gene][0]
            peak_ID = gene + "_" + name + "_" + line.split("\t")[7].split("(")[0] 

            if gene in genes:
                print gene 
                if strand == genes[gene][1]:
                    with open(filename + ".FINAL.bed",'a') as output:   
                        output.write(chrom + "\t" + s + "\t" + e + "\t" + peak_ID + "\t" + score + "\t" + strand + "\n")    
                else:
                    with open(filename + ".REMOVE.bed",'a') as output:                      
                        output.write(chrom + "\t" + s + "\t" + e + "\t" + peak_ID + "\t" + score + "\t" + strand + "\n")
            else:
                with open(filename + ".REMOVE.bed",'a') as output:
                    output.write(chrom + "\t" + s + "\t" + e + "\t" + peak_ID + "\t" + score + "\t" + strand + "\n")

def usage(filename):
    genes = {}

    with open(filename + ".FINAL.UNIQ.bed", 'r') as f:
        for line in f:
            line = line.rstrip()
            chrom, s, e, name, score, strand = line.split("\t")
            gene = name.split("_")[2]
            score = int(score)
            if gene not in genes:
                genes[gene] = score
            else:
                genes[gene] += score

    with open(filename + ".FINAL.UNIQ.bed", 'r') as f:
            for line in f:
                    line = line.rstrip()
                    chrom, s, e, name, score, strand = line.split("\t")
                    coverage = score 
                    gene= name.split("_")[2]
                    score = float(score) / float(genes[gene])

                    with open(filename + ".USAGE.bed", 'a') as output:
                        output.write(chrom + "\t" + s + "\t" + e + "\t" + name + "\t" + str(score) + "\t" + strand + "\t" + coverage + "\n") 


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-p", "--prefix", dest="prefix", default = 'test',
        help="input and output prefix (default test)")

    parser.add_option("-r", "--reference", dest="reference",
        help="reference annotations (e.g. Refseq2Gene.txt)")

    parser.add_option("-a", "--annotate", dest="annotate",
        help="fixing annotations mode")

    parser.add_option("-u", "--usage", dest="usage",
        help="calculating PAUs mode")

    (options, args) = parser.parse_args()
    
    main(options)


