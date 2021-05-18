#!/usr/bin/env python3
"""
Structure variation annotation using annovar
Suitable for sniffles or nanosv output vcf files
"""

import sys
import argparse
import os
import subprocess


import sv_vcf


header_bed_lastcols_list = ["SVTYPE", "SVID", "SVLEN", "RE/SUPP", "INFO"]

FLANK_SIZE = 10 # 断点上下游

def bed_to_ainput(bed, prefix):
    """
    convert bed into annovar input
    """
    prefix = prefix+".ainput"
    out = open(prefix, "w")
    with open(bed, "r") as bedio:
        for line in bedio:
            fields = line.strip().split("\t")[:5]
            fields[1] = str(int(fields[1])+1)
            fields[2] = str(int(fields[2])+1)
            fields.insert(3,"0")
            fields.insert(3,"0")
            new_line = "\t".join(fields)
            print(new_line, file=out)
    out.close()
    return prefix # return ainput path


def bed_to_ainput_point(bed, prefix, side):
    """
    convert bed into annovar input
    """
    prefix = prefix+".%s.ainput" % side
    out = open(prefix, "w")
    if side == 'left':
        idx = 1
    elif side == 'right':
        idx = 2
    else:
        raise RuntimeError("side must be left or right.")

    with open(bed, "r") as bedio:
        for line in bedio:
            fields = line.strip().split("\t")[:5]
            if 'M' in fields[0]:
                fields[1] = str(int(fields[idx]))
                fields[2] = str(int(fields[idx]))
            else:
                fields[1] = str(int(fields[idx]) - FLANK_SIZE)
                fields[2] = str(int(fields[idx]) + FLANK_SIZE)
            fields.insert(3,"0")
            fields.insert(3,"0")
            new_line = "\t".join(fields)
            print(new_line, file=out)
    out.close()
    return prefix # return ainput path


def run_annovar_point(ainput, table_annovar, database, prefix, hg_version, threads):
    protocols = ",".join(['refGene',
    "genomicSuperDups",
    "rmsk",
    "cpgIslandExt",])

    argument = "--neargene 2000,,,"
    operations = "g,r,r,r" #hg19
    run_list = [
        "perl",
        table_annovar,
        ainput,
        database,
        "-buildver", hg_version,
        "-thread", threads,
        "-out", prefix,
        "-otherinfo", "-nastring", ".", "-remove",
        "-protocol", protocols,
        "--operation", operations,
        "--argument", argument
    ]
    print(' '.join(run_list))
    subprocess.run(run_list)
    os.remove(ainput)



def run_annovar(ainput, table_annovar, database, prefix, hg_version, threads):
    print("Processing: %s annotation" % hg_version)
    protocols = ",".join(["refGene",
    "cytoBand",
     "genomicSuperDups",
     "rmsk",
#     "cpgIslandExt",
    "wgRna"])

    operations = "g,r,r,r,r" #hg19
    argument = "--neargene 2000,,,,"
    subprocess.run([
        "perl",
        table_annovar,
        ainput,
        database,
        "-buildver", hg_version,
        "-thread", threads,
        "-out", prefix,
        "-otherinfo", "-nastring", ".", "-remove",
        "-protocol", protocols,
        "--operation", operations,
        "--argument", argument
    ])
    os.remove(ainput)




def table_header(table,prefix,ref_version):
    """
    replace Otherinfo with SVTYPE ID -SVLEN Support_reads_num comments
    """
    out = open(prefix+".{}_multianno.xls".format(ref_version), "w")
    with open(table,"r") as table_io:
        header = table_io.readline()
        fields = header.strip().split("\t") #[,,'Otherinfo1', 'Otherinfo2', 'Otherinfo3', 'Otherinfo4', 'Otherinfo5'], 20200322 changed. old annovar version:just Ohterinfo1
        for i in range(len(header_bed_lastcols_list)):
            fields.pop(-1)
        fields.extend(header_bed_lastcols_list)
        print("\t".join(fields), file=out)
        for line in table_io:
            print(line, file=out, end="")
    out.close()



def get_args():
    parser = argparse.ArgumentParser(description="Structure variation annotation using annovar",usage="usage: %(prog)s [options]")
    parser.add_argument("--bed",
        help="bed file [default %(default)s]", metavar="FILE")
    parser.add_argument("--db",default="/sfs-grand-med-research/database/annovar/humandb",
        help="Database directory [default %(default)s]",metavar="DIR")
    parser.add_argument("--annovar",default="/sfs-grand-med-research/database/annovar/table_annovar.pl",
        help="Table Annovar [default %(default)s]",metavar="DIR")

    parser.add_argument("--prefix",
        help="output prefix [default %(default)s]", metavar="STR")
    parser.add_argument("--ref_version",
        help="ref version/hg19 or hg38 [default %(default)s]", metavar="STR")
    parser.add_argument("--threads",default=1, type=int,help="number of threads [default %(default)s]")

    parser.add_argument("--outdir", default=".", # default current work dir
        help="output directory [default %(default)s]", metavar="DIR")

    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    return parser.parse_args()



def main():
    args = get_args()

    # make dir
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    prefix_with_dir = os.path.join(args.outdir, args.prefix)

    bed = args.bed

    ref_version = args.ref_version

    ainput = bed_to_ainput(bed,  "{}.{}".format(prefix_with_dir, ref_version))
    run_annovar(ainput, args.annovar, args.db, prefix_with_dir, ref_version,str(args.threads))
    # ZH5_ZH6.hg38_multianno.txt
    # output '{}.{}_multianno.txt'.format(prefix_with_dir, ref_version)

    ainput_point_left = bed_to_ainput_point(bed, '{}.{}'.format(prefix_with_dir, ref_version), 'left')
    ainput_point_right = bed_to_ainput_point(bed, '{}.{}'.format(prefix_with_dir, ref_version), 'right')

    run_annovar_point(ainput_point_left, args.annovar, args.db, prefix_with_dir + '.left', ref_version, str(args.threads))
    run_annovar_point(ainput_point_right, args.annovar, args.db, prefix_with_dir + '.right', ref_version, str(args.threads))
     # output '{}.{}_multianno.txt'.format(prefix_with_dir, ref_version)


    #table_header(prefix_with_dir+".{}_multianno.txt".format(args.ref_version), prefix_with_dir,args.ref_version)

   # os.remove(prefix_with_dir+".{}_multianno.txt".format(args.ref_version))


if __name__ == '__main__':
    main()
