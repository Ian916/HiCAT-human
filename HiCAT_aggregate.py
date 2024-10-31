# hor.repeatnumber.xls
# out_all_layer.xls
# out_top_layer.xls
# out_hor.normal.fa
# out_hor.raw.fa
# out_final_hor.xls

import matplotlib.pyplot as plt
import matplotlib.patches as mpathes
from matplotlib.lines import Line2D
import numpy as np
import argparse

def Plot(monomer_sequence, patterns,pattern_static,
         block_seuqence, outdir, show_number = 5, show_min_repeat_number = 10):
    fig, ax = plt.subplots(figsize=(10, 10))
    monomer_len = len(monomer_sequence)
    color = '#D14524'
    custom_lines = []
    legend_text = []

    filter_patterns = {}
    pattern_count = 0
    for i in patterns.keys():
        if pattern_count >= show_number:
            break
        pattern_name = pattern_static[i][0]
        pattern_repeat_number = pattern_static[i][1]
        if pattern_repeat_number < show_min_repeat_number:
            continue
        filter_patterns[i] = patterns[i]
        pattern_count += 1

    re_patterns = list(filter_patterns.keys())[::-1]
    pattern_count = 0

    for i in re_patterns:
        # print(pattern_static[i])
        pattern_name = pattern_static[i][0]
        xy = np.array([0, pattern_count * monomer_len / 25])
        rect = mpathes.Rectangle(xy, monomer_len, monomer_len / 50, color='#D0CECE')
        ax.add_patch(rect)
        custom_lines.append(Line2D([0], [0], color=color, lw=2))
        legend_text.append(i)
        for j in patterns[i]:
            start = j[0]
            end = j[1]
            xy2 = np.asarray([start, pattern_count * monomer_len / 25])
            rect = mpathes.Rectangle(xy2, end + 1 - start, monomer_len / 50, color=color,lw=0)
            ax.add_patch(rect)
        plt.text(monomer_len + monomer_len / 50, pattern_count * monomer_len / 25, pattern_name, fontsize=10)
        pattern_count += 1

    xy3 = np.asarray([0, -monomer_len / 50])
    rect = mpathes.Rectangle(xy3, monomer_len, monomer_len / 1000, color='black')
    ax.add_patch(rect)
    point_bar = int(monomer_len / 10)
    for i in range(10):
        xy3 = np.asarray([0 + i * point_bar, -monomer_len / 50])
        rect = mpathes.Rectangle(xy3, monomer_len / 1000, -monomer_len / 100, color='black')
        ax.add_patch(rect)
        plt.text(0 + i * point_bar, -monomer_len / 50 - monomer_len / 50, str(block_seuqence[0 + i * point_bar][0]), fontsize=5)

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # ax.legend(custom_lines,legend_text)
    plt.xticks([])
    plt.yticks([])
    plt.axis('equal')
    plt.savefig(outdir + '/plot_pattern.pdf')
    plt.close()

    x = np.arange(min(len(filter_patterns.keys()),show_number))
    y = []
    y1 = []

    bar_width = 0.35

    tick_label = []
    for i in filter_patterns.keys():
        database = patterns[i]
        pattern = i.split('_')
        canonical = 0
        nested = 0
        for j in database:
            item_len = int(j[1]) + 1 - int(j[0])
            if item_len == len(pattern) * j[4]:
                canonical += j[4]
            else:
                nested += j[4]
        y.append(canonical)
        y1.append(nested)
        tick_label.append(pattern_static[i][0])

    pattern_static_file = outdir + '/pattern_static.xls'
    pattern_static_file = open(pattern_static_file,'w')
    pattern_static_file.write('HORs\tCanonical\tNested\n')

    for i in range(len(tick_label)):
        pattern_static_file.write(tick_label[i]+'\t'+str(y[i]) +'\t' +str(y1[i]) + '\n')
    pattern_static_file.close()

    plt.figure(figsize=(10, 10))
    plt.bar(x, y, bar_width, align="center", color="c", label="canonical", alpha=0.5)
    plt.bar(x + bar_width, y1, bar_width, color="b", align="center", label="nested", alpha=0.5)

    plt.xlabel("HORs")
    plt.ylabel("Repeat Number")

    plt.xticks(x + bar_width / 2, tick_label)

    plt.legend()

    plt.savefig(outdir + '/pattern_static.pdf')
    plt.close()

def readPattern(pattern_file):
    patterns = {}
    key = ''
    patterns_name = {}
    index = 1
    with open(pattern_file,'r') as pf:
        while True:
            line = pf.readline()[:-1]
            if not line:
                break
            key = line
            pattern_len = len(key.split('_'))
            pattern_name = 'R' + str(index) + 'L' + str(pattern_len)
            patterns_name[key] = pattern_name
            patterns[key] = []
            line = pf.readline()[:-1]
            itemsets = line.split('\t')[1:]

            # 增加strand
            for i in itemsets:
                r = i.split(',')
                start = int(r[0])
                end = int(r[1])
                strand = r[2]
                pattern = r[3]
                repeat_number = int(r[4])
                patterns[key].append([start,end,strand,pattern,repeat_number])
            index += 1
    return patterns,patterns_name

def readMonomerSequence(monomer_sequence_file):
    monomer_sequence = []
    with open(monomer_sequence_file,'r') as msf:
        while True:
            line = msf.readline()[:-2]
            if not line:
                break
            items = line.split('\t')
            monomer_sequence = items[1].split(' ')
    return monomer_sequence

def buildMonomerBlockSequence(decomposition_path):
    block_sequence = []
    with open(decomposition_path, 'r') as bf:
        while True:
            line = bf.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            if items[1][-1] == "'":
                block_sequence.append([int(items[2]), int(items[3]) ,'-'])
            else:
                block_sequence.append([int(items[2]), int(items[3]) ,'+'])
    return block_sequence

def reverse(sequence):
    base_map = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    new_sequence = ''
    for i in sequence[::-1]:
        new_sequence += base_map[i]
    return new_sequence

def buildHORFile(patterns, patterns_name, pattern_name_table,pattern_mono_table,base_sequence,
                 monomer_sequence, block_sequence, outdir):
    out_hor_raw_file = outdir + '/out_hor.raw.merge.fa'
    out_hor_raw_file = open(out_hor_raw_file,'w')
    out_hor_normal_file = outdir + '/out_hor.normal.merge.fa'
    out_hor_normal_file = open(out_hor_normal_file,'w')
    for i in patterns.keys():
        pattern_name = pattern_name_table[patterns_name[i]]
        pattern = pattern_mono_table[pattern_name].split('_')
        database = patterns[i]
        # ([start,end,strand,pattern,repeat_number])
        for j in database:
            start = j[0]
            end = j[1]
            strand = j[2] # 更新增加strand
            monomer_sequence_item = monomer_sequence[start:end+1]
            repeat_number = j[4]
            # patternname.index start end pattern repeatnumber rawpattern
            monomer_sequence_item_str = ''
            for k in monomer_sequence_item:
                monomer_sequence_item_str += str(k) + '_'
            monomer_sequence_item_str = str(monomer_sequence_item_str[:-1])
            out_hor_raw_file.write('>' + pattern_name  + '::' +
                                       str(block_sequence[start][0]) + '-' + str(block_sequence[end][1] + 1) +
                                   '::' + strand +
                                       ' nHOR-' + pattern_mono_table[pattern_name] + '::rHOR-' + monomer_sequence_item_str +
                                   ' repeat_number-' +str(repeat_number)+ '\n')
            out_hor_raw_file.write(base_sequence[block_sequence[start][0]:block_sequence[end][1] + 1] + '\n')


            out_hor_normal_file.write('>' + pattern_name + '::' +
                                   str(block_sequence[start][0]) + '-' + str(block_sequence[end][1] + 1) +
                                      '::' + strand +
                                   ' nHOR-' + pattern_mono_table[pattern_name] + '::rHOR-' + monomer_sequence_item_str +
                                      ' repeat_number-' +str(repeat_number)+ '\n')

            if len(pattern) == 1:
                # 考虑反链 '-' 链标准化变正
                normal_sequence = base_sequence[block_sequence[start][0]:block_sequence[end][1] + 1]
                if strand == '-':
                    normal_sequence = reverse(normal_sequence)
                out_hor_normal_file.write(normal_sequence + '\n')
            else:
                if strand == '-':
                    monomer_sequence_item = monomer_sequence_item[::-1] # 反链序列翻转
                    double_sequence = monomer_sequence_item + monomer_sequence_item
                    double_index = list(range(len(monomer_sequence_item)))[::-1] + \
                                   list(range(len(monomer_sequence_item)))[::-1] # 反链index翻转
                    count = 0
                    prefix = []
                    pattern_index = 0
                    for k in range(len(double_sequence)):
                        if pattern[pattern_index] == double_sequence[k]:
                            prefix.append([k, double_sequence[k], double_index[k], pattern_index])
                    normal_pattern = []
                    for k in prefix:
                        record = [k]
                        pattern_index = k[3] + 1
                        not_find = 0
                        for l in range(k[0] + 1, len(double_sequence)):
                            if double_sequence[l] == pattern[pattern_index]:
                                record.append([l, double_sequence[l], double_index[l], pattern_index])
                                pattern_index += 1
                                if pattern_index == len(pattern):
                                    break
                            else:
                                continue_flag = 0
                                for m in record:
                                    if double_sequence[l] == m[1]:
                                        continue_flag = 1
                                if continue_flag == 1:
                                    continue
                                else:
                                    not_find = 1
                                    break
                        if not_find == 1:
                            continue
                        if len(record) != len(pattern):
                            continue
                        normal_pattern = record
                    normal_sequence = ''
                    for k in normal_pattern:
                        block_start = block_sequence[start + k[2]][0]
                        block_end = block_sequence[start + k[2]][1] + 1
                        normal_sequence += reverse(base_sequence[block_start:block_end]) # 每个block变反
                    out_hor_normal_file.write(normal_sequence + '\n')
                else:
                    # +
                    double_sequence = monomer_sequence_item + monomer_sequence_item
                    double_index = list(range(len(monomer_sequence_item))) + list(range(len(monomer_sequence_item)))
                    count = 0
                    prefix = []
                    pattern_index = 0
                    for k in range(len(double_sequence)):
                        if pattern[pattern_index] == double_sequence[k]:
                            prefix.append([k, double_sequence[k], double_index[k], pattern_index])
                    normal_pattern = []
                    for k in prefix:
                        record = [k]
                        pattern_index = k[3] + 1
                        not_find = 0
                        for l in range(k[0] + 1, len(double_sequence)):
                            if double_sequence[l] == pattern[pattern_index]:
                                record.append([l, double_sequence[l], double_index[l], pattern_index])
                                pattern_index += 1
                                if pattern_index == len(pattern):
                                    break
                            else:
                                continue_flag = 0
                                for m in record:
                                    if double_sequence[l] == m[1]:
                                        continue_flag = 1
                                if continue_flag == 1:
                                    continue
                                else:
                                    not_find = 1
                                    break

                        if not_find == 1:
                            continue
                        if len(record) != len(pattern):
                            continue
                        normal_pattern = record

                    normal_sequence = ''
                    for k in normal_pattern:
                        block_start = block_sequence[start + k[2]][0]
                        block_end = block_sequence[start + k[2]][1]+1
                        normal_sequence += base_sequence[block_start:block_end]

                    out_hor_normal_file.write(normal_sequence + '\n')

    out_hor_raw_file.close()
    out_hor_normal_file.close()



def main():
    parser = argparse.ArgumentParser(description="Aggregate HOR between ref and reads")
    parser.add_argument("-p", "--pattern_file", help="All patterns in HiCATReads with name, required", required=True)
    parser.add_argument("-c", "--target_chr", help="Target chromosome", required=True)
    parser.add_argument("-i", "--hicat_outdir", help="HiCAT result path", required=True)
    
    parser.add_argument("-sp", "--show_hor_number", help="Default visualized the top five HORs", type=int, default=5,
                        required=False)
    parser.add_argument("-sn", "--show_hor_min_repeat_number",
                        help="Default visualized the HORs with repeat numbers greater than 10", type=int, default=10,
                        required=False)
    
    
    args = parser.parse_args()
    
    pattern_file = args.pattern_file
    target_chr = args.target_chr
    hicat_outdir = args.hicat_outdir
    
    show_hor_number = args.show_hor_number
    show_hor_min_repeat_number = args.show_hor_min_repeat_number

    pattern_table = {}
    pattern_mono_table = {}
    with open(pattern_file,'r') as pf:
        while True:
            line = pf.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            info = items[0].split('_')
            chr = info[0]
            pattern = info[1]
            mon_pattern = items[1]
            if chr != target_chr:
                continue
            pattern_table[mon_pattern] = pattern
            pattern_mono_table[pattern] = mon_pattern

    final_hor_file = hicat_outdir + '/out_final_hor.xls'
    pattern_name_table = {}
    index = 1
    with open(final_hor_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            f.readline()
            # 在table中查找
            ori_monomer_pattern = line
            mon_pattern = ori_monomer_pattern.split('_')
            in_flag = 0
            pattern_name = 'R' + str(index) + 'L' + str(len(mon_pattern))
            index += 1
            for i in range(len(mon_pattern)):  # 循环pattern
                prefix_pattern = mon_pattern[i:]
                suffix_pattern = mon_pattern[:i]
                loop_pattern = prefix_pattern + suffix_pattern
                s_loop_pattern = ''
                for j in loop_pattern:
                    s_loop_pattern += str(j) + '_'
                s_loop_pattern = s_loop_pattern[:-1]
                if s_loop_pattern in pattern_table.keys():
                    in_flag = 1
                    pattern_name_table[pattern_name] = pattern_table[s_loop_pattern]
                    break
            # 没找到，还是原始命名，还是原始标准化
            if in_flag == 0:
                pattern_name_table[pattern_name] = pattern_name
                pattern_mono_table[pattern_name] = ori_monomer_pattern


    hor_repeatnumber_file = hicat_outdir + '/hor.repeatnumber.xls'
    out_hor_repeatnumber_file = hicat_outdir + '/hor.repeatnumber.merge.xls'
    out_hor_repeatnumber_file = open(out_hor_repeatnumber_file,'w')
    out_hor_repeatnumber_file.write('HORs\tRepeatNumber\n')
    with open(hor_repeatnumber_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if line.startswith('HOR'):
                continue
            items = line.split('\t')
            out_hor_repeatnumber_file.write(pattern_name_table[items[0]]+'\t'+items[1]+'\n')
    out_hor_repeatnumber_file.close()

    all_layer_file = hicat_outdir + '/out_all_layer.xls'
    out_all_layer_file = hicat_outdir + '/out_all_layer.merge.xls'
    out_all_layer_file = open(out_all_layer_file,'w')
    with open(all_layer_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            out_all_layer_file.write(pattern_name_table[items[0]] + '\t' +
                                     items[1] + '\t' +
                                     items[2] + '\t' +
                                     items[3] + '\t' +
                                     items[4] + '\t' +
                                     items[5] + '\n')
    out_all_layer_file.close()


    top_layer_file = hicat_outdir + '/out_top_layer.xls'
    out_top_layer_file = hicat_outdir + '/out_top_layer.merge.xls'
    out_top_layer_file = open(out_top_layer_file, 'w')
    with open(top_layer_file, 'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            out_top_layer_file.write(pattern_name_table[items[0]] + '\t' +
                                     items[1] + '\t' +
                                     items[2] + '\t' +
                                     items[3] + '\t' +
                                     items[4] + '\n')
    out_top_layer_file.close()


    # 读取mono seq
    mon_file = hicat_outdir + '/out_monomer_seq.xls'
    monomer_sequence = ''
    with open(mon_file,'r') as f:
        items = f.readline()[:-2].split('\t')
        monomer_sequence = items[1].split(' ')
    base_sequence = ''
    base_file = hicat_outdir + '/input_fasta.1.fa'
    with open(base_file,'r') as f:
        f.readline()
        base_sequence = f.readline()

    # 读取final sd
    # base_sequence[block_sequence[start][0]:block_sequence[end][1] + 1]
    # monomer_sequence_item = monomer_sequence[start:end+1]
    sd_file = hicat_outdir + '/' + '/final_decomposition.tsv'
    block_sequence = []
    with open(sd_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            start = int(items[2])
            end = int(items[3])
            strand = '+'
            if items[1].endswith("'"):
                strand = '-'
            block_sequence.append([start,end,strand])

    patterns,pattern_name = readPattern(final_hor_file)  # 更新增加strand [start,end,strand,pattern,repeat_number]
    buildHORFile(patterns, pattern_name, pattern_name_table,pattern_mono_table,base_sequence,
                 monomer_sequence, block_sequence, hicat_outdir)


    # 作图

    monomer_sequence_file = hicat_outdir + '/out_monomer_seq.xls'
    monomer_sequence = readMonomerSequence(monomer_sequence_file)
    decomposition_path = hicat_outdir + '/final_decomposition.tsv'
    block_sequence = buildMonomerBlockSequence(decomposition_path)

    pattern_static = {}  #
    pattern_index = 1
    for i in patterns.keys():
        pattern = i.split('_')
        database = patterns[i]  # 增加了strand
        repeat_number = 0
        for j in database:
            repeat_number += j[4]
        pattern_name = 'R' + str(pattern_index) + 'L' + str(len(pattern))
        pattern_name = pattern_name_table[pattern_name]
        pattern_static[i] = [pattern_name, repeat_number]
        pattern_index += 1

    Plot(monomer_sequence, patterns, pattern_static, block_sequence, hicat_outdir,
         show_number=show_hor_number, show_min_repeat_number=show_hor_min_repeat_number)



if __name__ == '__main__':
    main()