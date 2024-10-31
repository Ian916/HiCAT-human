import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse


def main():
    parser = argparse.ArgumentParser(description="Normal sample HOR numbers")
    parser.add_argument("-i", "--all_sample_dir", help="All sample result path, required", required=True)
    parser.add_argument("-n", "--base_number_file", help="Reads total base number", required=True)
    parser.add_argument("-s", "--sample_file", help="The file record all sample names, required", required=True)
    parser.add_argument("-m", "--HOR_matrix_file", help="matrix record all HOR numbers in each sample, from aggregate output", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-ref", "--ref_genome_size", help="reference human genome size (base number), default is the base number of chm13v2.0", default=3054815472, type=int, required=False)
    parser.add_argument("-r", "--rare_ratio", help="exclude HORs with less frequency than this ratio in each sample, default is 0.1", default=0.1, type=float, required=False)
    parser.add_argument("-u", "--plot_upper", help="mean fold-change upper bound in the output plot, default is 5", type=float, default=5, required=False)

    args = parser.parse_args()
    workdir = args.all_sample_dir
    base_number_file = args.base_number_file
    sample_file = args.sample_file
    HOR_matrix_file = args.HOR_matrix_file
    outdir = args.outdir

    chm13_genome_size = 3054815472
    rare_ratio = 0.1
    plot_upper = 5

    sample_cov = {}
    with open(base_number_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            sample_cov[items[0]] = int(items[2]) / chm13_genome_size
    sample_gender = {}
    with open(sample_file, 'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            sample_gender[items[0]] = items[1]

    HOR_matrix = pd.read_csv(HOR_matrix_file,sep='\t',index_col=0)
    patterns = HOR_matrix.index.tolist()
    chr_patterns_list = {}
    for i in patterns:
        items = i.split('_')
        chr = items[0]
        pattern = items[1]
        index = int(pattern[1:].split('L')[0])
        if chr not in chr_patterns_list.keys():
            chr_patterns_list[chr] = []
            chr_patterns_list[chr].append([index,pattern])
        else:
            chr_patterns_list[chr].append([index, pattern])


    data_dict = HOR_matrix.to_dict()


    new_data_dict = {}
    for i in data_dict.keys():
        new_data_dict[i] = data_dict[i]
        for j in new_data_dict[i].keys():
            chr = j.split('_')[0]
            if chr == 'X':
                if sample_gender[i] == 'female':
                    new_data_dict[i][j] = new_data_dict[i][j] / sample_cov[i]
                else:
                    new_data_dict[i][j] = (new_data_dict[i][j] * 2) / sample_cov[i]
            elif chr == 'Y':
                if sample_gender[i] == 'female':
                    new_data_dict[i][j] = 0
                else:
                    new_data_dict[i][j] = (new_data_dict[i][j] * 2) / sample_cov[i]
            else:
                new_data_dict[i][j] = new_data_dict[i][j] / sample_cov[i]

    normal_data = pd.DataFrame(new_data_dict)
    HOR_nnumber_matrix_file = outdir + '/HOR.nnumber.xls'
    normal_data.to_csv(HOR_nnumber_matrix_file,sep='\t')

    # filter rare HOR
    pattern = normal_data.index.tolist()
    chr_patterns = {}
    for i in pattern:
        items = i.split('_')
        if items[0] not in chr_patterns.keys():
            chr_patterns[items[0]] = [i]
        else:
            chr_patterns[items[0]].append(i)
    normal_data = normal_data.T
    chr_save_patterns = {}
    for i in chr_patterns.keys():
        chr_save_patterns[i] = []
        chr_data = normal_data[chr_patterns[i]]
        chr_data = chr_data.T
        samples = chr_data.columns.tolist()
        sample_sum = []
        for j in samples:
            sumnumber = np.sum(chr_data[j])
            sample_sum.append(sumnumber)
        patterns = chr_data.index.tolist()
        np_chr_data = np.asarray(chr_data)
        for j in range(len(patterns)):
            ok = 0
            for k in range(len(sample_sum)):
                if np_chr_data[j][k] == 0:
                    ratio = 0
                else:
                    ratio = np_chr_data[j][k] / sample_sum[k]
                if ratio > rare_ratio:
                    ok = 1
            if ok == 1:
                index = int(patterns[j].split('_')[1][1:].split('L')[0])
                chr_save_patterns[i].append([index,patterns[j]])
        chr_save_patterns[i] = sorted(chr_save_patterns[i],key=lambda x:x[0])

    save_patterns = []
    for i in chr_save_patterns.keys():
        for j in chr_save_patterns[i]:
            save_patterns.append(j[1])

    non_rare_normal_data = normal_data[save_patterns]
    non_rare_HOR_nnumber_matrix_file = outdir + '/non-rare.HOR.nnumber.xls'
    non_rare_normal_data.to_csv(non_rare_HOR_nnumber_matrix_file, sep='\t')

    # mean FC
    # 额外算Y
    non_rare_normal_data_table = non_rare_normal_data.to_dict()

    HOR_mean = {}

    for i in non_rare_normal_data_table.keys():
        HOR_mean[i] = 0
        sample_number = 0
        for j in non_rare_normal_data_table[i].keys():
            if i.startswith('Y_M'):
                if sample_gender[j] == 'male':
                    sample_number += 1
                    HOR_mean[i] += non_rare_normal_data_table[i][j]
            else:
                sample_number += 1
                HOR_mean[i] += non_rare_normal_data_table[i][j]
        if sample_number == 0:
            HOR_mean[i] = 0
        else:
            HOR_mean[i] = HOR_mean[i] / sample_number

    non_rare_normal_meanFC_data_table = {}
    for i in non_rare_normal_data_table.keys():
        non_rare_normal_meanFC_data_table[i] = {}
        for j in non_rare_normal_data_table[i].keys():
            non_rare_normal_meanFC_data_table[i][j] = non_rare_normal_data_table[i][j] / HOR_mean[i]

    non_rare_normal_meanFC_data = pd.DataFrame(non_rare_normal_meanFC_data_table).T
    non_rare_HOR_nnumber_meanFC_matrix_file = outdir + '/non-rare.HOR.nnumber.meanFC.xls'
    non_rare_normal_meanFC_data.to_csv(non_rare_HOR_nnumber_meanFC_matrix_file, sep='\t')

    samples = non_rare_normal_meanFC_data.columns.tolist()
    patterns = non_rare_normal_meanFC_data.index.tolist()
    non_rare_normal_meanFC_data = np.asarray(non_rare_normal_meanFC_data)

    long_data = []
    long_data_filter = []
    for i in range(len(patterns)):
        for j in range(len(samples)):
            long_data.append([patterns[i], samples[j], non_rare_normal_meanFC_data[i][j]])
            if non_rare_normal_meanFC_data[i][j] > plot_upper:
                long_data_filter.append([patterns[i], samples[j], plot_upper])
            else:
                long_data_filter.append([patterns[i], samples[j], non_rare_normal_meanFC_data[i][j]])
    long_data = pd.DataFrame(long_data, columns=[['HOR', 'sample', 'meanFC']])
    plot_file = outdir + '/non-rare.HOR.nnumber.meanFC.long.xls'
    long_data.to_csv(plot_file, sep='\t', index=None)
    plot_filter_file = outdir + '/non-rare.HOR.nnumber.meanFC.long.filter.xls'
    long_data_filter = pd.DataFrame(long_data_filter, columns=[['HOR', 'sample', 'meanFC']])
    long_data_filter.to_csv(plot_filter_file, sep='\t', index=None)

    plt.figure(figsize=(len(patterns),len(patterns) / 5))

    for i in range(len(samples)):
        sample = samples[i]
        x = []
        y = []
        for j in range(len(patterns)):
            pattern = patterns[j]
            if pattern.startswith('Y_M'):
                if sample_gender[sample] == 'male':
                    x.append(pattern)
                    if non_rare_normal_meanFC_data[j][i] > plot_upper:
                        y.append(plot_upper)
                    else:
                        y.append(non_rare_normal_meanFC_data[j][i])
            else:
                x.append(pattern)
                if non_rare_normal_meanFC_data[j][i] > plot_upper:
                    y.append(plot_upper)
                else:
                    y.append(non_rare_normal_meanFC_data[j][i])
        plt.plot(x, y,color='grey')

    plt.title('HOR variation (upper = '+str(plot_upper)+')',fontsize = len(patterns) / 2)
    plt.xlabel('HOR',fontsize = len(patterns) / 2)
    plt.ylabel('mean FC',fontsize = len(patterns) / 2)
    plt.xticks(rotation=45,fontsize = len(patterns) / 5)
    plt.yticks(fontsize=len(patterns) / 5)
    plt.tight_layout()
    plt.savefig(workdir + '/HOR_variation.pdf')
    plt.close()

    # 计算SD
    HOR_sd = []
    for i in non_rare_normal_meanFC_data_table.keys():
        values = []
        for j in non_rare_normal_meanFC_data_table[i].keys():
            if i.startswith('Y_M'):
                if sample_gender[j] == 'male':
                    values.append(non_rare_normal_meanFC_data_table[i][j])
            else:
                values.append(non_rare_normal_meanFC_data_table[i][j])
        values = np.asarray(values)
        HOR_sd.append([i,values.std()])
    HOR_sd = sorted(HOR_sd,key=lambda x:x[1],reverse=True)
    HOR_sd_file = outdir + '/non-rare.HOR.nnumber.meanFC.std.xls'
    HOR_sd = pd.DataFrame(HOR_sd, columns=[['HOR', 'std']])
    HOR_sd.to_csv(HOR_sd_file, sep='\t', index=None)



if __name__ == '__main__':
    main()

