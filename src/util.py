"""
Utilities for splicing analysis.
"""

import gc

def build_cvg_table(f):

    table = {}
    f.seek(0)
    for line in f.readlines():
        splitted = line.split("\t")
        key = splitted[3].strip()

        table[key] = {"reads": splitted[-1].strip(),
                      "length": str(abs(int(splitted[1]) - int(splitted[2]))),
                      "start": splitted[1].strip(),
                      "stop": splitted[2].strip()}

    return table

def build_junc_table(junc_file):
    table = {}

    junc_file.seek(0)
    for line in junc_file.readlines():
        splitted = line.split("\t")
        key = splitted[3].strip()

        if splitted[11] == "-1":
            # No reads for this feature
            table[key] = {"reads": 0}
            continue

        if key in table:
            table[key]["reads"] += int(splitted[14].strip())
        else:
            table[key] = {"reads": int(splitted[14].strip())}

    return table

def get_chr_list(wao_file):
    # covert to list later
    chr_table = {}
    wao_file.seek(0)

    for line in wao_file:
        chr_table[line.split("\t")[0].strip()] = None

    print ("\t Built chromosome list.")
    return chr_table.keys()

def build_chr_table(wao_file, chr_name):
    """
    Builds table for given chromosome from
    a -wao intersectBed file, assuming it is
    sorted by chromosome name.
    """

    wao_file.seek(0)
    table = {}
    found = False
    for line in wao_file:
        splitted = line.split("\t")

        if splitted[0].strip() == chr_name:
            found = True
            key = splitted[3].strip()

            if key in table:
                table[key][splitted[13].strip()] = None
            else:
                table[key] = {splitted[13].strip(): None}

        else:
            if found:
                return table

    return table

def build_exon_intron_table(exon_wao, intron_wao, chr_file):
    exon_wao.seek(0)
    intron_wao.seek(0)
    table = {}

    print("\t Building chromosome list....")
    chr_list = get_chr_list(chr_file)

    done = 1
    total = len(chr_list)
    for active_chr in chr_list:
        print("\t Building %d out of %d chromosomes..." % (done, total))
        exon_table = build_chr_table(exon_wao, active_chr)
        intron_table = build_chr_table(intron_wao, active_chr)

        for intron_key in intron_table:
            exon1_key = intron_key.split(":")
            intron_id = exon1_key[1]
            intron_num = int(exon1_key[2])
            exon1_key[0] = "exon"
            exon2_key = exon1_key[:]
            exon1_key = ":".join(exon1_key)
            exon2_key[2] = str(intron_num + 1)
            exon2_key = ":".join(exon2_key)

            intron_reads = intron_table[intron_key]
            exon1_reads = exon_table[exon1_key]
            exon2_reads = exon_table[exon2_key]

            exon_intron = [read for read in exon1_reads if read in intron_reads and read not in exon2_reads]
            intron_exon = [read for read in exon2_reads if read in intron_reads and read not in exon1_reads]

            table["exon-intron:%s:%d" % (intron_id, intron_num)] = {"reads": len(exon_intron)}
            table["intron-exon:%s:%d" % (intron_id, intron_num)] = {"reads": len(intron_exon)}

        del exon_table
        del intron_table
        gc.collect()
        done += 1

    return table

def build_fpkm_table(fpkm_file):
    """
    Takes fpkm file and builds a preprocessing
    table.

    Note that fpkm values should be the 10th column.

    Args:
         fpkm_file: open handle to gene_fpkm file from
                    cufflinks / cuffdiff
    """
    table = {}
    fpkm_file.seek(0)
    for line in fpkm_file.readlines():
        splitted = line.split("\t")
        gene_name = splitted[3].strip()
        fpkm = splitted[9].strip()

        if gene_name in table:
            print (">> gene.fpkm_tracking has duplicate genes! Exiting...")
            exit(1)

        table[gene_name] = {"fpkm": fpkm}

    return table

def build_gene_length_table(gene_file):
    """
    Builds pre-processing table from file
    containing all gene lengths.

    Args:
         gene_file: open handle to file containing
                    all gene lengths.
    """

    table = {}
    gene_file.seek(0)
    for line in gene_file.readlines():
        splitted = line.split("\t")
        gene_name = splitted[0].strip()
        length = splitted[1].strip()

        if gene_name in table:
            print(">> gene length file has duplicate rows! Exiting...")
            exit(1)

        table[gene_name] = {"length": length}

    return table

def build_mrna_table(mrna_file):
    """
    Builds pre-processing table for mrna
    length.

    Args:
        mrna_file: open handle to mrna file.
                   When using hg38_ref, the file
                   is hg38_select_feature_lengths


    """
    table = {}
    mrna_file.seek(0)

    for line in mrna_file.readlines():
        splitted = line.split("\t")
        gene_ID = splitted[0].strip()
        mrna_length = splitted[2].strip()

        if gene_ID in table:
            print( ">> Duplicate gene_ID in %s. Exiting..." % mrna_file.name)
            exit(1)

        table[gene_ID] = {"length": mrna_length}

    return table

'''
- Reference file builders.
'''
def build_exon_exon_junction_ref(exonref_file, output_file):
    """
    Builds exon-exon junctions from a
    reference exon file (i.e hg38_only_exons_refGene_longest.gff3).
    Called once per genome version.
    Outputs in .gff format.

    Args:
         exonref_file: open handle to reference file GFF
                       containing all exonic information.
         output_file: open handle to output file.

    """

    table = {}
    exonref_file.seek(0)

    for line in exonref_file.readlines():
        exonid = line.split("\t")[8].split(";")[0][8:].strip()
        table[exonid] = line.split("\t")[:8]

    for exon_id1, data1 in table.items():
        split = exon_id1.split(":")
        exon_id2 = split[0] + ":" + str(int(split[1]) + 1)

        if exon_id2 in table:
            data2 = table[exon_id2]
        else:
            continue

        if int(data1[3]) < int(data2[3]):
            start = data1[3]
            stop = data2[4]
            final_id = "ID=exon-exon:%s:%s:%d" % (split[0], split[1], int(split[1]) + 1)
        else:
            start = data2[3]
            stop = data1[4]
            final_id = "ID=exon-exon:%s:%d:%s" % (split[0], int(split[1]) + 1, split[1])


        final_row = data1[:3] + [start, stop] + data1[5:8] + [final_id]
        output_file.write("\t".join(final_row) + "\n")

def build_exon_intron_junction_ref(exonref_file, output_file):
    """
    Builds exon-intron junctions from a
    reference exon file (i.e hg38_only_exons_refGene_longest.gff3).
    Called once per genome version.
    Outputs in .gff format.

    Args:
         exonref_file: open handle to reference GFF file
                       containing all exonic information.
         output_file: open handle to output file.
    """

    table = {}
    exonref_file.seek(0)

    for line in exonref_file.readlines():
        exonid = line.split("\t")[8].split(";")[0][8:].strip()
        table[exonid] = line.split("\t")[:8]

    for exon_id1, data1 in table.items():
        split = exon_id1.split(":")
        exon_id2 = split[0] + ":" + str(int(split[1]) + 1)

        if exon_id2 in table:
            data2 = table[exon_id2]
        else:
            continue

        if int(data1[3]) < int(data2[3]):
            final_id1 = "ID=exon-intron:%s:%s" % (split[0], split[1])
            start1 = int(data1[4]) - 1
            stop1 = int(data1[4]) + 2
            final_id2 = "ID=intron-exon:%s:%s" % (split[0], split[1])
            start2 = int(data2[3]) - 1
            stop2 = int(data2[3]) + 2
        else:
            final_id1 = "ID=exon-intron:%s:%s" % (split[0], split[1])
            start1 = int(data2[4]) - 1
            stop1 = int(data2[4]) + 2
            final_id2 = "ID=intron-exon:%s:%s" % (split[0], split[1])
            start2 = int(data2[3]) - 1
            stop2 = int(data2[3]) + 2

        final_row1 = data1[:3] + [str(start1), str(stop1)] + data1[5:8] + [final_id1]
        final_row2 = data2[:3] + [str(start2), str(stop2)] + data2[5:8] + [final_id2]
        output_file.write("\t".join(final_row1) + "\n")
        output_file.write("\t".join(final_row2) + "\n")

'''
- Miscellaneous
'''
def append_id(map_file, bed_file, col):
    """
    Takes coveragebed ouput and appends gene_id name
    to the last column.

    Args:
         map_file: open handle to file containing gene_id -> name mappings
         bed_file: open handle to .bed file to be converted."
         col: 0-indexed column number containing geneid, for eg,
              the column should contain entries like `ID=exon-exon:NR_024540:3:2`

    Returns:
         Open file handle to new .bed file with IDs.
    """
    id_table  = {}
    missing   = 0

    bed_file.seek(0)
    map_file.seek(0)
    lines = map_file.readlines()

    for line in lines:
        id_table[line.split("\t")[0].strip()] = line.split("\t")[1].strip()

    output = open(bed_file.name + "_ID", "w+")
    lines = bed_file.readlines()
    for line in lines:
        if len(line.split("\t")[col].split(":")) < 2:
            continue

        gene_id = line.split("\t")[col].split(":")[1].strip()
        if gene_id in id_table:
            output.write(line.strip() + "\t" + id_table[gene_id] + "\n")
        else:
            output.write(line.strip() + "\t" + "-" + "\n")
            missing +=1

    print( ">> Total missing for %s : %s" % (bed_file.name, missing))
    return output

def check_intron_integrity(intron_line, intron_key,
                           junc_table_key, junc_table,
                           up_exon_key, down_exon_key, exon_table, sample_name):
    """
    Takes a splitted line from the intron file and verifies
    positional information. Mostly for sanity checking.
    """

    intron_start = int(intron_line[1])
    intron_stop = int(intron_line[2])

    junc_start = int(junc_table[junc_table_key]["start"])
    junc_stop = int(junc_table[junc_table_key]["stop"])

    up_exon_start = int(exon_table[up_exon_key]["start"])
    up_exon_stop = int(exon_table[up_exon_key]["stop"])

    down_exon_start = int(exon_table[down_exon_key]["start"])
    down_exon_stop = int(exon_table[down_exon_key]["stop"])

    if not (junc_start == up_exon_start <= up_exon_stop <= intron_start <= intron_stop <=
            down_exon_start <= down_exon_stop == junc_stop):
        print (">> The following assertion failed:")
        print (">> junc_start == up_exon_start <= up_exon_stop <= intron_start <= intron_stop <= down_exon_start <= down_exon_stop == junc_stop")
        print (">> %s == %s <= %s <= %s <= %s <= %s <= %s == %s" % \
       (junc_start , up_exon_start , up_exon_stop , intron_start , intron_stop , down_exon_start , down_exon_stop , junc_stop))
        exit(1)

def compute_total_cvg(refgene_file, output_file):
    """
    Takes intersectBed refGene vs accepted_hits.bam output
    and computes total coverage into `output_file`
    """
    refgene_file.seek(0)
    table = {}

    for line in refgene_file:
        splitted = line.split("\t")
        feature = splitted[7].strip()
        if feature in table:
            table[feature] += int(splitted[10])
        else:
            table[feature] = int(splitted[10])

    total = 0
    for feature in table:
        output_file.write("%s\t%d\n" % (feature, table[feature]))
        total += table[feature]

    intergenic = total - table["gene"]
    output_file.write("intergenic\t%d\n" % intergenic)
    output_file.write("Total\t%d" % total)

def do_analysis(intron_file, exon_table,
                exon_exon_junc_table, exon_intron_junc_table,
                fpkm_table, gene_length_table, mrna_table, out_file, sample_name):
    """
    Takes cvgbed intron file and pre-processing tables
    and  builds the analysis file.

    Args:
    intron_file: open handle to cvgbed intron file.
    exon_table: pre-processed exon table
    intron_table: pre-processed intron table
    out_file: open handle to output file.

    """

    out_file.write("gene_name\tgene_ID:intron_number\tintron_length\t")
    out_file.write("upstream_exon_length\tdownstream_exon_length\t")
    out_file.write("gene_length\tmrna_length\t")
    out_file.write(sample_name + "_intron_reads\t")
    out_file.write(sample_name + "_upstream_exon_reads\t")
    out_file.write(sample_name + "_downstream_exon_reads\t")
    out_file.write(sample_name + "_flanking_junc_reads\t")
    out_file.write(sample_name + "_exon_intron_reads\t")
    out_file.write(sample_name + "_intron_exon_reads\t")
    out_file.write(sample_name + "_gene_FPKM\t")
    out_file.write(sample_name + "_psi_value\n")

    # out_file.write("gene_name\tgene_ID\tintron_number\tintron_reads\tintron_length\t")
    # out_file.write("flanking_junc_reads\texon_intron_reads\tintron_exon_reads\t")
    # out_file.write("downstream_exon_reads\tupstream_exon_reads\t")
    # out_file.write("downstream_exon_length\tupstream_exon_length\tgene_FPKM\tgene_length\tmrna_length\n")

    intron_file.seek(0)
    for line in intron_file.readlines():
        splitted = line.split("\t")

        key = splitted[3].strip().split(":")[1:]
        key = ":".join(key).strip()
        gene_ID = key.split(":")[0].strip()
        gene_name = splitted[-1].strip()

        intron_number = key.split(":")[1].strip()

        # exon-exon numbering is dependent on ordering of
        # coordinates. For eg: intron:1 can belong to
        # exon-exon:1:2 or exon-exon:2:1
        up_exon_key =  "exon:" + gene_ID + ":" + intron_number
        down_exon_key = "exon:" + gene_ID + ":" + str(int(intron_number)+1)
        exon_exon_junc_table_key = "exon-exon:%s:%s:%d" % (gene_ID, intron_number, int(intron_number) + 1)
        exon_intron_key = "exon-intron:%s:%s" % (gene_ID, intron_number)
        intron_exon_key = "intron-exon:%s:%s" % (gene_ID, intron_number)
        if exon_exon_junc_table_key not in exon_exon_junc_table:
            up_exon_key = "exon:" + gene_ID + ":" + str(int(intron_number)+1)
            down_exon_key = "exon:" + gene_ID + ":" + intron_number
            exon_exon_junc_table_key = "exon-exon:%s:%d:%s" % (gene_ID, int(intron_number) + 1, intron_number)
            if exon_exon_junc_table_key not in exon_exon_junc_table:
                print (">> %s does not have flanking junctions! Exiting..." % exon_exon_junc_table_key)
                exit(1)

        intron_reads = splitted[10].strip()
        intron_length = str(abs(int(splitted[1]) - int(splitted[2])))
        flank_junc_reads = str(exon_exon_junc_table[exon_exon_junc_table_key]["reads"])
        exon_intron_reads = str(exon_intron_junc_table[exon_intron_key]["reads"])
        intron_exon_reads = str(exon_intron_junc_table[intron_exon_key]["reads"])
        up_exon_reads = exon_table[up_exon_key]["reads"].strip()
        down_exon_reads = exon_table[down_exon_key]["reads"].strip()
        up_exon_len = exon_table[up_exon_key]["length"].strip()
        down_exon_len = exon_table[down_exon_key]["length"].strip()

        if gene_name not in fpkm_table:
            gene_fpkm = "0"
        else:
            gene_fpkm = fpkm_table[gene_name]["fpkm"].strip()

        gene_length = gene_length_table[gene_ID]["length"].strip()
        mrna_length = mrna_table[gene_ID]["length"].strip()

        if (float(intron_reads) + float(exon_intron_reads) + float(intron_exon_reads)) == 0:
            psi_value = "Inf"
        else:
            psi_value = (float(up_exon_reads) + float(down_exon_reads) + float(flank_junc_reads)) \
                    / (float(intron_reads) + float(exon_intron_reads) + float(intron_exon_reads))

#        check_intron_integrity(splitted, key,
#                               exon_exon_junc_table_key, exon_exon_junc_table,
#                               up_exon_key, down_exon_key, exon_table)

        out = gene_name + "\t" + gene_ID  + ":" + intron_number + "\t" + intron_length + "\t"  \
        + up_exon_len + "\t" + down_exon_len + "\t" \
        + gene_length + "\t" + mrna_length + "\t" \
        + intron_reads + "\t" + up_exon_reads + "\t" + down_exon_reads + "\t" + flank_junc_reads + "\t" \
        + exon_intron_reads + "\t" + intron_exon_reads + "\t" + gene_fpkm + "\t" + str(psi_value) + "\n"

        out_file.write(out)

if __name__ == "__main__":
    build_exon_intron_junction_ref(open("/home/data/hg38_ref/hg38_only_exons_refGene_longest.gff3"), open("./hg38_exon_intron_junctions.bed", "w+"))

