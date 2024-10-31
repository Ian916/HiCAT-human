import argparse
import os
from joblib import Parallel, delayed


def getReads(cen_reads_output_dir, file_list):
    reads = {}
    header = ''
    for i in file_list:
        read_file = cen_reads_output_dir + '/' + i
        with open(read_file, 'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                if line.startswith('>'):
                    header = line[1:].split(' ')[0]
                else:
                    reads[header] = line
    return reads


def getLabels(cen_reads_output_dir, cen_label_files):
    reads_label = {}
    chr_reads = {
        'chr1': [], 'chr2': [], 'chr3': [], 'chr4': [], 'chr5': [],
        'chr6': [], 'chr7': [], 'chr8': [], 'chr9': [], 'chr10': [],
        'chr11': [], 'chr12': [], 'chr13': [], 'chr14': [], 'chr15': [],
        'chr16': [], 'chr17': [], 'chr18': [], 'chr19': [], 'chr20': [],
        'chr21': [], 'chr22': [], 'chrX': [], 'chrY': []
    }
    for i in cen_label_files:
        label_file = cen_reads_output_dir + '/' + i
        with open(label_file, 'r') as f:
            f.readline()
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                readname = items[0].split(' ')[0]
                label = items[1]
                reads_label[readname] = label
                chr_reads[label].append(readname)
    return reads_label, chr_reads


def classifation(file, split_reads_dir, predict_script, model1_file, model2_file, featurename_file,
                 cen_reads_output_dir):
    prefix = str(file.split('.')[0])
    print(file)
    print(prefix)
    reads_file = split_reads_dir + '/' + file
    cmd = 'python ' + \
          predict_script + ' ' + \
          '-r' + ' ' + reads_file + ' ' + \
          '-m1' + ' ' + model1_file + ' ' + \
          '-m2' + ' ' + model2_file + ' ' + \
          '-f' + ' ' + featurename_file + ' ' + \
          '-p' + ' ' + prefix + ' ' + \
          '-o' + ' ' + cen_reads_output_dir
    os.system(cmd)
    print(cmd)


def main():
    parser = argparse.ArgumentParser(description="Finding centromere reads and classifying into each chromosome")
    parser.add_argument("-i", "--input_fasta", help="input fasta, required", required=True)
    parser.add_argument("-rn", "--number_of_reads", help="the number of reads, required", required=True)
    parser.add_argument("-o", "--output_dir", help="HiCAT reads output path default is ./HiCAT_out",
                        default='./HiCAT_out',
                        required=False)

    parser.add_argument("-th", "--thread", help="The number of threads, default is 1", type=int, default=1,
                        required=False)

    args = parser.parse_args()

    input_fasta = args.input_fasta
    number_of_reads = args.number_of_reads
    output_dir = args.output_dir
    thread = args.thread

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    script_path = os.path.split(os.path.realpath(__file__))[0]
    split_reads_script = script_path + '/' + 'splitReadsFa.sh'
    if thread == 1:
        cmd = 'sh ' + \
              split_reads_script + ' ' + \
              input_fasta + ' ' + output_dir + ' ' + str(thread) + ' ' + number_of_reads
    else:
        cmd = 'sh ' + \
              split_reads_script + ' ' + \
              input_fasta + ' ' + output_dir + ' ' + str(thread - 1) + ' ' + number_of_reads
    print('split reads file')
    print(cmd)
    os.system(cmd)

    split_reads_dir = output_dir + '/split_fa'
    if not os.path.exists(split_reads_dir):
        os.mkdir(split_reads_dir)
    script_path = os.path.split(os.path.realpath(__file__))[0]
    model1_file = script_path + '/model.pkl'
    model2_file = script_path + '/model_cen.pkl'
    predict_script = script_path + '/predictReads.py'
    featurename_file = script_path + '/kmer_feature.txt'
    file_list = os.listdir(split_reads_dir)
    read_file_list = []
    for i in file_list:
        if i.startswith('split_fasta') and i.endswith('.fasta'):
            read_file_list.append(i)
    print('classify reads')
    cen_reads_output_dir = output_dir + '/reads_out'
    if not os.path.exists(cen_reads_output_dir):
        os.mkdir(cen_reads_output_dir)

    Parallel(n_jobs=thread)(delayed(classifation)(i, split_reads_dir,
                                                  predict_script,
                                                  model1_file,
                                                  model2_file,
                                                  featurename_file,
                                                  cen_reads_output_dir) for i in read_file_list)

    print('out results')
    file_list = os.listdir(cen_reads_output_dir)
    cen_fa_files = []
    cen_label_files = []
    for i in file_list:
        if i.endswith('cen.fa'):
            cen_fa_files.append(i)
        if i.endswith('label.xls'):
            cen_label_files.append(i)

    cen_reads = getReads(cen_reads_output_dir, cen_fa_files)
    reads_label, chr_reads = getLabels(cen_reads_output_dir, cen_label_files)

    out_reads_file = cen_reads_output_dir + '/merge.fa'
    out_reads_file = open(out_reads_file, 'w')
    for i in cen_reads:
        out_reads_file.write('>' + i + '\n')
        out_reads_file.write(cen_reads[i] + '\n')
    out_reads_file.close()
    out_label_file = cen_reads_output_dir + '/merge.label.xls'
    out_label_file = open(out_label_file, 'w')
    out_label_file.write('readname\tlabel\n')
    for i in reads_label:
        out_label_file.write(i + '\t' + reads_label[i] + '\n')
    out_label_file.close()

    out_result_dir = cen_reads_output_dir + '/chr_data'
    if not os.path.exists(out_result_dir):
        os.mkdir(out_result_dir)
    for i in chr_reads.keys():
        out_chr_reads_file = out_result_dir + '/' + i + '.fa'
        out_chr_reads_file = open(out_chr_reads_file, 'w')
        for j in chr_reads[i]:
            out_chr_reads_file.write('>' + j + '\n')
            out_chr_reads_file.write(cen_reads[j] + '\n')
        out_chr_reads_file.close()


if __name__ == '__main__':
    main()
