import argparse
import os

script_path = os.path.split(os.path.realpath(__file__))[0]

def hicat_human_reads(args):
    reads_dir = args.input_reads_dir
    sample_list = args.reads_sample_file
    workdir = args.reads_output_dir
    thread = args.thread

    hicat_reads_classification = script_path + '/HiCATReads_classification.py'
    hicat_reads_annotation = script_path + '/HiCATReads_annotation.py'
    ref_dir = script_path + '/HiCAT_ref'

    cmd = 'mkdir -p ' + workdir
    os.system(cmd)

    with open(workdir + '/read_count', 'w') as o:
        with open(sample_list, 'r') as f:
            for line in f:
                items = line.strip().split('\t')
                ID = items[0]
                fai_file = reads_dir + '/' + ID + '.fasta.fai'
                read_number = 0
                total_base_number = 0
                with open(fai_file, 'r') as fai:
                    for _ in fai:
                        read_number += 1
                        length = int(_.split('\t')[1])
                        total_base_number += length
                o.writelines(ID + '\t' + str(read_number) + '\t' + str(total_base_number) + '\n')

    read_number_file = workdir + '/read_count'
    with open(read_number_file, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            ID = items[0]
            count = items[1]

            cmd = 'mkdir -p ' + workdir + '/' + ID
            print(cmd)
            os.system(cmd)
            outdir = workdir + '/' + ID
            input_file = reads_dir + '/' + ID + '.fasta'
            cmd = 'python ' + hicat_reads_classification + ' -i ' + input_file + ' -rn ' + count + ' -o ' + outdir + ' -th ' + str(thread)
            print(cmd)
            os.system(cmd)
            cmd = 'python ' + hicat_reads_annotation + ' -i ' + outdir + '/reads_out' + ' -r ' + ref_dir + ' -th ' + str(thread)
            print(cmd)
            os.system(cmd)

            cmd = 'rm -rf ' + outdir + '/split_fa'
            os.system(cmd)
            cmd = 'rm -rf ' + outdir + '/reads_out/split_fasta_*'
            os.system(cmd)

def hicat_human_assembly(args):
    assembly_dir = args.input_assembly_dir
    sample_list = args.assembly_sample_file
    thread = args.thread
    min_similarity = args.min_similarity
    max_hor_len = args.max_hor_len
    show_hor_number = args.show_hor_number
    show_hor_min_repeat_number = args.show_hor_min_repeat_number

    hicat_mini = script_path + '/HiCAT_mini.py'

    alpha_satellite_template = script_path + '/AlphaSat.fa'
    lastz_script = script_path + '/processLastz.py'
    region_script = script_path + '/getCENRegion.py'
    cen_script = script_path + '/prepareAlphaSat.py'

    with open(sample_list, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            ID = items[0]
            assembly_file = assembly_dir + '/' + ID + '.fasta'
            fai_file = assembly_dir + '/' + ID + '.fasta.fai'

            cmd = 'mkdir -p ' + assembly_dir + '/' + ID
            os.system(cmd)
            workdir = assembly_dir + '/' + ID

            cmd = 'python ' + lastz_script + ' -r ' + assembly_file + ' -fi ' + fai_file + ' -w ' + workdir + ' -a ' + alpha_satellite_template
            os.system(cmd)

            cmd = 'python ' + region_script + ' -fi ' + fai_file + ' -w ' + workdir
            os.system(cmd)

            cmd = 'python ' + cen_script + ' -r ' + assembly_file + ' -f ' + workdir + '/lastz_regions.bed' + ' -w ' + workdir
            os.system(cmd)

            censeq_dir = workdir + '/censeq'
            chr_file_list = os.listdir(censeq_dir)
            for i in chr_file_list:
                if i.startswith('chr'):
                    cen_seq = censeq_dir + '/' + i + '/' + i + '.cen.fa'
                    template_file = script_path + '/HiCAT_ref/' + i[3:] + '/out/out_monomer_represent.fa'
                    outdir = censeq_dir + '/' + i
                    cmd = 'python ' + hicat_mini + ' -i ' + cen_seq + ' -t ' + alpha_satellite_template + ' -tm ' + template_file + ' -o ' + outdir + ' -th ' + str(thread) + \
                          ' -ms ' + str(min_similarity) + ' -mh ' + str(max_hor_len) + ' -sp ' + str(show_hor_number) + ' -sn ' + str(show_hor_min_repeat_number)
                    os.system(cmd)

def hicat_human_reads_aggregate(args):
    reads_dir = args.all_sample_dir
    sample_list = args.sample_file
    ref_genome_size = args.ref_genome_size
    rare_ratio = args.rare_ratio
    plot_upper = args.plot_upper

    hicat_reads_aggregate = script_path + '/HiCATReads_sampleAggregate.py'
    pattern_summary = script_path + '/getPatternSummary.py'
    pattern_table = script_path + '/getPatternTable.py'
    normal_HOR = script_path + '/normalHORnumber.py'

    cmd = 'python ' + hicat_reads_aggregate + ' -i ' + reads_dir + ' -s ' + sample_list
    os.system(cmd)

    #pattern summary and pattern table
    cmd = 'python ' + pattern_summary + ' -i ' + reads_dir + ' -s ' + sample_list
    os.system(cmd)
    cmd = 'python ' + pattern_table + ' -i ' + reads_dir
    os.system(cmd)

    #normal HOR number and plot
    cmd = 'python ' + normal_HOR + ' -i ' + reads_dir + ' -n ' + reads_dir + '/read_count' + ' -s ' + sample_list + ' -m ' + reads_dir + '/all.sample.HOR.repeatnumber.matrix.xls' \
          + ' -o ' + reads_dir + ' -ref ' + str(ref_genome_size) + ' -r ' + str(rare_ratio) + ' -u ' + str(plot_upper)
    os.system(cmd)


def hicat_human_assembly_match_to_reads(args):
    pattern_file = args.pattern_file
    sample_list = args.sample_file
    assembly_dir = args.hicat_human_assembly_dir
    show_hor_number = args.show_hor_number
    show_hor_min_repeat_number = args.show_hor_min_repeat_number

    hicat_mini_aggregate = script_path + '/HiCAT_aggregate.py'

    chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                '19', '20', '21', '22', 'X', 'Y']
    with open(sample_list, 'r') as f:
        for line in f:
            ID = line.strip()
            for target_chr in chr_list:
                hicat_assembly_out_dir = assembly_dir + '/' + ID + '/censeq/chr' + target_chr
                if os.path.exists(hicat_assembly_out_dir):
                    cmd = 'python ' + hicat_mini_aggregate + ' -p ' + pattern_file + ' -c ' + target_chr + ' -i ' + \
                          hicat_assembly_out_dir + ' -sp ' + str(show_hor_number) + ' -sn ' + str(show_hor_min_repeat_number)
                    os.system(cmd)

def hicat_human_only_reads(args):
    reads_dir = args.input_reads_dir
    reads_sample_list = args.reads_sample_file
    reads_workdir = args.reads_output_dir
    thread = args.thread

    hicat_reads_classification = script_path + '/HiCATReads_classification.py'
    hicat_reads_annotation = script_path + '/HiCATReads_annotation.py'
    ref_dir = script_path + '/HiCAT_ref'

    cmd = 'mkdir -p ' + reads_workdir
    os.system(cmd)

    with open(reads_workdir + '/read_count', 'w') as o:
        with open(reads_sample_list, 'r') as f:
            for line in f:
                items = line.strip().split('\t')
                ID = items[0]
                fai_file = reads_dir + '/' + ID + '.fasta.fai'
                read_number = 0
                total_base_number = 0
                with open(fai_file, 'r') as fai:
                    for _ in fai:
                        read_number += 1
                        length = int(_.split('\t')[1])
                        total_base_number += length
                o.writelines(ID + '\t' + str(read_number) + '\t' + str(total_base_number) + '\n')

    read_number_file = reads_workdir + '/read_count'
    with open(read_number_file, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            ID = items[0]
            count = items[1]

            cmd = 'mkdir -p ' + reads_workdir + '/' + ID
            print(cmd)
            os.system(cmd)
            outdir = reads_workdir + '/' + ID
            input_file = reads_dir + '/' + ID + '.fasta'
            cmd = 'python ' + hicat_reads_classification + ' -i ' + input_file + ' -rn ' + count + ' -o ' + outdir + ' -th ' + str(
                thread)
            print(cmd)
            os.system(cmd)
            cmd = 'python ' + hicat_reads_annotation + ' -i ' + outdir + '/reads_out' + ' -r ' + ref_dir + ' -th ' + str(
                thread)
            print(cmd)
            os.system(cmd)

            cmd = 'rm -rf ' + outdir + '/split_fa'
            os.system(cmd)
            cmd = 'rm -rf ' + outdir + '/reads_out/split_fasta_*'
            os.system(cmd)

    #aggregate
    reads_output_dir = reads_workdir
    ref_genome_size = args.ref_genome_size
    rare_ratio = args.rare_ratio
    plot_upper = args.plot_upper

    hicat_reads_aggregate = script_path + '/HiCATReads_sampleAggregate.py'
    pattern_summary = script_path + '/getPatternSummary.py'
    pattern_table = script_path + '/getPatternTable.py'
    normal_HOR = script_path + '/normalHORnumber.py'

    cmd = 'python ' + hicat_reads_aggregate + ' -i ' + reads_output_dir + ' -s ' + reads_sample_list
    os.system(cmd)

    # pattern summary and pattern table
    cmd = 'python ' + pattern_summary + ' -i ' + reads_output_dir + ' -s ' + reads_sample_list
    os.system(cmd)
    cmd = 'python ' + pattern_table + ' -i ' + reads_output_dir
    os.system(cmd)

    # normal HOR number and plot
    cmd = 'python ' + normal_HOR + ' -i ' + reads_output_dir + ' -n ' + read_number_file + ' -s ' + reads_sample_list + ' -m ' + reads_output_dir + '/all.sample.HOR.repeatnumber.matrix.xls' \
          + ' -o ' + reads_output_dir + ' -ref ' + str(ref_genome_size) + ' -r ' + str(rare_ratio) + ' -u ' + str(plot_upper)
    os.system(cmd)

def hicat_human_reads_with_assembly(args):
    #
    # reads
    reads_dir = args.input_reads_dir
    reads_sample_list = args.reads_sample_file
    reads_workdir = args.reads_output_dir
    thread = args.thread

    hicat_reads_classification = script_path + '/HiCATReads_classification.py'
    hicat_reads_annotation = script_path + '/HiCATReads_annotation.py'
    ref_dir = script_path + '/HiCAT_ref'

    cmd = 'mkdir -p ' + reads_workdir
    os.system(cmd)

    with open(reads_workdir + '/read_count', 'w') as o:
        with open(reads_sample_list, 'r') as f:
            for line in f:
                items = line.strip().split('\t')
                ID = items[0]
                fai_file = reads_dir + '/' + ID + '.fasta.fai'
                read_number = 0
                total_base_number = 0
                with open(fai_file, 'r') as fai:
                    for _ in fai:
                        read_number += 1
                        length = int(_.split('\t')[1])
                        total_base_number += length
                o.writelines(ID + '\t' + str(read_number) + '\t' + str(total_base_number) + '\n')

    read_number_file = reads_workdir + '/read_count'
    with open(read_number_file, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            ID = items[0]
            count = items[1]

            cmd = 'mkdir -p ' + reads_workdir + '/' + ID
            print(cmd)
            os.system(cmd)
            outdir = reads_workdir + '/' + ID
            input_file = reads_dir + '/' + ID + '.fasta'
            cmd = 'python ' + hicat_reads_classification + ' -i ' + input_file + ' -rn ' + count + ' -o ' + outdir + ' -th ' + str(
                thread)
            print(cmd)
            os.system(cmd)
            cmd = 'python ' + hicat_reads_annotation + ' -i ' + outdir + '/reads_out' + ' -r ' + ref_dir + ' -th ' + str(
                thread)
            print(cmd)
            os.system(cmd)

            cmd = 'rm -rf ' + outdir + '/split_fa'
            os.system(cmd)
            cmd = 'rm -rf ' + outdir + '/reads_out/split_fasta_*'
            os.system(cmd)

    #
    # assembly
    assembly_dir = args.input_assembly_dir
    assembly_sample_list = args.assembly_sample_file
    thread = args.thread
    min_similarity = args.min_similarity
    max_hor_len = args.max_hor_len
    show_hor_number = args.show_hor_number
    show_hor_min_repeat_number = args.show_hor_min_repeat_number

    hicat_mini = script_path + '/HiCAT_mini.py'

    alpha_satellite_template = script_path + '/AlphaSat.fa'
    lastz_script = script_path + '/processLastz.py'
    region_script = script_path + '/getCENRegion.py'
    cen_script = script_path + '/prepareAlphaSat.py'

    with open(assembly_sample_list, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            ID = items[0]
            assembly_file = assembly_dir + '/' + ID + '.fasta'
            fai_file = assembly_dir + '/' + ID + '.fasta.fai'

            cmd = 'mkdir -p ' + assembly_dir + '/' + ID
            os.system(cmd)
            workdir = assembly_dir + '/' + ID

            cmd = 'python ' + lastz_script + ' -r ' + assembly_file + ' -fi ' + fai_file + ' -w ' + workdir + ' -a ' + alpha_satellite_template
            os.system(cmd)

            cmd = 'python ' + region_script + ' -fi ' + fai_file + ' -w ' + workdir
            os.system(cmd)

            cmd = 'python ' + cen_script + ' -r ' + assembly_file + ' -f ' + workdir + '/lastz_regions.bed' + ' -w ' + workdir
            os.system(cmd)

            censeq_dir = workdir + '/censeq'
            chr_file_list = os.listdir(censeq_dir)
            for i in chr_file_list:
                if i.startswith('chr'):
                    cen_seq = censeq_dir + '/' + i + '/' + i + '.cen.fa'
                    template_file = script_path + '/HiCAT_ref/' + i[3:] + '/out/out_monomer_represent.fa'
                    outdir = censeq_dir + '/' + i
                    cmd = 'python ' + hicat_mini + ' -i ' + cen_seq + ' -t ' + alpha_satellite_template + ' -tm ' + template_file + ' -o ' + outdir + ' -th ' + str(thread) + \
                          ' -ms ' + str(min_similarity) + ' -mh ' + str(max_hor_len) + ' -sp ' + str(show_hor_number) + ' -sn ' + str(show_hor_min_repeat_number)
                    os.system(cmd)

    #
    #reads aggregate
    reads_output_dir = reads_workdir
    ref_genome_size = args.ref_genome_size
    rare_ratio = args.rare_ratio
    plot_upper = args.plot_upper

    hicat_reads_aggregate = script_path + '/HiCATReads_sampleAggregate.py'
    pattern_summary = script_path + '/getPatternSummary.py'
    pattern_table = script_path + '/getPatternTable.py'
    normal_HOR = script_path + '/normalHORnumber.py'

    cmd = 'python ' + hicat_reads_aggregate + ' -i ' + reads_output_dir + ' -s ' + reads_sample_list
    os.system(cmd)

    # pattern summary and pattern table
    cmd = 'python ' + pattern_summary + ' -i ' + reads_output_dir + ' -s ' + reads_sample_list
    os.system(cmd)
    cmd = 'python ' + pattern_table + ' -i ' + reads_output_dir
    os.system(cmd)

    # normal HOR number and plot
    cmd = 'python ' + normal_HOR + ' -i ' + reads_output_dir + ' -n ' + read_number_file + ' -s ' + reads_sample_list + ' -m ' + reads_output_dir + '/all.sample.HOR.repeatnumber.matrix.xls' \
          + ' -o ' + reads_output_dir + ' -ref ' + str(ref_genome_size) + ' -r ' + str(rare_ratio) + ' -u ' + str(plot_upper)
    os.system(cmd)

    #
    # assembly match to reads
    pattern_file = reads_output_dir + '/pattern.table.xls'

    hicat_mini_aggregate = script_path + '/HiCAT_aggregate.py'

    chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                '19', '20', '21', '22', 'X', 'Y']
    with open(assembly_sample_list, 'r') as f:
        for line in f:
            ID = line.strip()
            for target_chr in chr_list:
                hicat_assembly_out_dir = assembly_dir + '/' + ID + '/censeq/chr' + target_chr
                if os.path.exists(hicat_assembly_out_dir):
                    cmd = 'python ' + hicat_mini_aggregate + ' -p ' + pattern_file + ' -c ' + target_chr + ' -i ' + \
                          hicat_assembly_out_dir + ' -sp ' + str(show_hor_number) + ' -sn ' + str(show_hor_min_repeat_number)
                    os.system(cmd)

def main():
    parser = argparse.ArgumentParser(description="automatically annotate centromere HOR patterns from reads and assemblies")

    subparsers = parser.add_subparsers(help="Modes: reads,assembly,reads_aggregate,assembly_match,only_reads,reads_with_assembly")

    # python HiCAT_human.py reads
    parser_reads = subparsers.add_parser('reads', help='Annotate centromere HOR patterns from reads')
    parser_reads.add_argument("-r", "--input_reads_dir", help="Input reads directory containing all sample read file (.fasta) and corresponding .fai file, required", required=True)
    parser_reads.add_argument("-rs", "--reads_sample_file", help="File record all sample names (the prefix of .fasta file), one sample name per line, required", required=True)
    parser_reads.add_argument("-o", "--reads_output_dir", help="HiCAT human reads output path, required", required=True)
    parser_reads.add_argument("-th", "--thread", help="The number of threads, default is 1", type=int, default=1, required=False)
    parser_reads.set_defaults(func=hicat_human_reads)

    # python HiCAT_human.py assembly
    parser_assembly = subparsers.add_parser('assembly', help='Annotate centromere HOR patterns from assembly')
    parser_assembly.add_argument("-a", "--input_assembly_dir", help="Input assembly directory containing all sample assembly file (.fasta) corresponding .fai file, required", required=True)
    parser_assembly.add_argument("-as", "--assembly_sample_file", help="File record all sample names (the prefix of .fasta file), one sample name per line, required", required=True)
    parser_assembly.add_argument("-ms", "--min_similarity", help="Lower bound for similarity threshold between block and target monomer, default is 0.9", type=float, default=0.9, required=False)
    parser_assembly.add_argument("-mh", "--max_hor_len", help="Upper bound of the tandem repeat unit length for improving efficiency, default 40 monomers", type=int, default=40, required=False)
    parser_assembly.add_argument("-sp", "--show_hor_number", help="Visualized HOR number, default 5", type=int, default=5, required=False)
    parser_assembly.add_argument("-sn", "--show_hor_min_repeat_number", help="HORs with repeat number less than this threshold will not be visualized, default 10", type=int, default=10, required=False)
    parser_assembly.add_argument("-th", "--thread", help="The number of threads, default is 1", type=int, default=1, required=False)
    parser_assembly.set_defaults(func=hicat_human_assembly)

    # python HiCAT_human.py reads_aggregate
    parser_reads_aggregate = subparsers.add_parser('reads_aggregate', help='Merge the reads HOR annotation results for multi-samples, and plot the normalized HOR number change graph among samples')
    parser_reads_aggregate.add_argument("-i", "--all_sample_dir", help="All sample HiCAT human reads result path, required", required=True)
    parser_reads_aggregate.add_argument("-s", "--sample_file", help="Two-column file record all sample names and gender, \\t separator, required", required=True)
    parser_reads_aggregate.add_argument("-ref", "--ref_genome_size", help="Reference human genome size (base number), default is the base number of CHM13", default=3054815472, type=int, required=False)
    parser_reads_aggregate.add_argument("-rr", "--rare_ratio", help="Exclude HORs with less frequency than this ratio in all samples, default is 0.1", default=0.1, type=float, required=False)
    parser_reads_aggregate.add_argument("-u", "--plot_upper", help="Mean fold-change upper bound in the output plot, default is 5", type=float, default=5, required=False)
    parser_reads_aggregate.set_defaults(func=hicat_human_reads_aggregate)

    # python HiCAT_human.py assembly_match
    parser_reads_aggregate_match = subparsers.add_parser('assembly_match', help='Match the assembly HOR annotation results for multi-samples based on reads annotation results')
    parser_reads_aggregate_match.add_argument("-p", "--pattern_file", help="pattern.table.xls, output file of reads_aggregate pipeline, record all patterns in HiCAT human reads result, required", required=True)
    parser_reads_aggregate_match.add_argument("-s", "--sample_file", help="File record all sample names, required", required=True)
    parser_reads_aggregate_match.add_argument("-a", "--hicat_human_assembly_dir", help="HiCAT human assembly result path", required=True)
    parser_reads_aggregate_match.add_argument("-sp", "--show_hor_number", help="Visualized HOR number, default 5", type=int, default=5, required=False)
    parser_reads_aggregate_match.add_argument("-sn", "--show_hor_min_repeat_number", help="HORs with repeat number less than this threshold will not be visualized, default 10", type=int, default=10, required=False)
    parser_reads_aggregate_match.set_defaults(func=hicat_human_assembly_match_to_reads)

    # python HiCAT_human.py only_reads
    parser_only_reads = subparsers.add_parser('only_reads', help='Reads annotation and aggregate')
    parser_only_reads.add_argument("-r", "--input_reads_dir", help="Input reads directory containing all sample read file (.fasta) and corresponding .fai file, required", required=True)
    parser_only_reads.add_argument("-rs", "--reads_sample_file", help="File record all sample names (the prefix of .fasta file), one sample name per line, required", required=True)
    parser_only_reads.add_argument("-o", "--reads_output_dir", help="HiCAT human reads output path, required", required=True)
    parser_only_reads.add_argument("-th", "--thread", help="The number of threads, default is 1", type=int, default=1, required=False)

    parser_only_reads.add_argument("-ref", "--ref_genome_size", help="Reference human genome size (base number), default is the base number of CHM13", default=3054815472, type=int, required=False)
    parser_only_reads.add_argument("-rr", "--rare_ratio", help="Exclude HORs with less frequency than this ratio in each sample, default is 0.1", default=0.1, type=float, required=False)
    parser_only_reads.add_argument("-u", "--plot_upper", help="Mean fold-change upper bound in the output plot, default is 5", type=float, default=5, required=False)

    parser_only_reads.set_defaults(func=hicat_human_only_reads)

    # python HiCAT_human.py reads_with_assembly
    parser_reads_with_assembly = subparsers.add_parser('reads_with_assembly', help='Annotation for reads and assemblies, aggregate the results')
    parser_reads_with_assembly.add_argument("-r", "--input_reads_dir", help="Input reads directory containing all sample read file (.fasta) and corresponding .fai file, required", required=True)
    parser_reads_with_assembly.add_argument("-rs", "--reads_sample_file", help="File record all sample names (the prefix of .fasta file), one sample name per line, required", required=True)
    parser_reads_with_assembly.add_argument("-o", "--reads_output_dir", help="HiCAT human reads output path, required", required=True)
    parser_reads_with_assembly.add_argument("-th", "--thread", help="The number of threads, default is 1", type=int, default=1, required=False)

    parser_reads_with_assembly.add_argument("-a", "--input_assembly_dir", help="Input assembly directory containing all sample assembly file (.fasta) corresponding .fai file, required", required=True)
    parser_reads_with_assembly.add_argument("-as", "--assembly_sample_file", help="File record all sample names (the prefix of .fasta file), one sample name per line, required", required=True)
    parser_reads_with_assembly.add_argument("-ms", "--min_similarity", help="Lower bound for similarity threshold between block and target monomer, default is 0.9", type=float, default=0.9, required=False)
    parser_reads_with_assembly.add_argument("-mh", "--max_hor_len", help="Upper bound of the tandem repeat unit length for improving efficiency, default 40 monomers", type=int, default=40, required=False)
    parser_reads_with_assembly.add_argument("-sp", "--show_hor_number", help="Visualized HOR number, default 5", type=int, default=5, required=False)
    parser_reads_with_assembly.add_argument("-sn", "--show_hor_min_repeat_number", help="HORs with repeat number less than this threshold will not be visualized, default 10", type=int, default=10, required=False)

    parser_reads_with_assembly.add_argument("-ref", "--ref_genome_size", help="Reference human genome size (base number), default is the base number of CHM13", default=3054815472, type=int, required=False)
    parser_reads_with_assembly.add_argument("-rr", "--rare_ratio", help="Exclude HORs with less frequency than this ratio in each sample, default is 0.1", default=0.1, type=float, required=False)
    parser_reads_with_assembly.add_argument("-u", "--plot_upper", help="Mean fold-change upper bound in the output plot, default is 5", type=float, default=5, required=False)

    parser_reads_with_assembly.set_defaults(func=hicat_human_reads_with_assembly)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
