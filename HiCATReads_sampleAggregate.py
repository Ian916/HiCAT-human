import argparse
import pandas as pd

def readPattern(pattern_file):
    patterns = {}
    patterns_info = {}
    with open(pattern_file,'r') as pf:
        while True:
            line = pf.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            pattern_name = items[0]
            pattern = items[1]
            all_repeat_number = int(items[2])
            patterns[pattern] = []
            line = pf.readline()[:-1]
            itemsets = line.split('\t')[1:]
            patterns_info[pattern] = [pattern_name,all_repeat_number]
            # 增加strand
            for i in itemsets:
                r = i.split(',')
                read_name = r[0]
                start = int(r[1])
                end = int(r[2])
                strand = r[3]
                r_pattern = r[4]
                repeat_number = int(r[5])
                patterns[pattern].append([read_name,start,end,strand,r_pattern,repeat_number])
    return patterns,patterns_info

def outPattern(patterns_info,patterns,patterns_name,out_pattern_file):
    out_pattern_file = open(out_pattern_file,'w')
    for i in patterns.keys():
        out_pattern_file.write(patterns_name[i][0]+'\t'+patterns_name[i][1]+'\t'+str(patterns_info[i][1])+'\n')
        out_pattern_file.write('P(read_name,start,end,strand,pattern,repeat_number):')
        for j in patterns[i]:
            out_pattern_file.write('\t' + j[0] + ',' + str(j[1]) + ',' + str(j[2]) + ',' +
                                   j[3] + ',' + j[4] + ',' + str(j[5]))
        out_pattern_file.write('\n')

def outAllLayer(all_layer_file,patterns_name,out_all_layer_file):
    out_all_layer_file = open(out_all_layer_file,'w')
    with open(all_layer_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            pattern = items[-3].split('_')
            in_flag = 0
            pattern_name = ''
            for i in range(len(pattern)):  # 循环pattern
                prefix_pattern = pattern[i:]
                suffix_pattern = pattern[:i]
                loop_pattern = prefix_pattern + suffix_pattern
                s_loop_pattern = ''
                for j in loop_pattern:
                    s_loop_pattern += str(j) + '_'
                s_loop_pattern = s_loop_pattern[:-1]
                if s_loop_pattern in patterns_name.keys():
                    in_flag = 1
                    pattern_name = patterns_name[s_loop_pattern][0]
                    break
            if in_flag == 0:
                r_pattern = pattern[::-1]  # 反向
                for i in range(len(r_pattern)):  # 循环r_pattern
                    prefix_pattern = r_pattern[i:]
                    suffix_pattern = r_pattern[:i]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for j in loop_pattern:
                        s_loop_pattern += str(j) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in patterns_name.keys():
                        in_flag = 1
                        pattern_name = patterns_name[s_loop_pattern][0]
                        break
            out_all_layer_file.write(items[0] + '\t' +
                                     items[1] +'\t' +
                                     items[2] + '\t' +
                                     items[3] + '\t' +
                                     items[4] + '\t' +
                                     pattern_name + '\t' + items[6]+'\n')
    out_all_layer_file.close()

def outTopLayer(top_layer_file,patterns_name,out_top_layer_file):
    out_top_layer_file = open(out_top_layer_file, 'w')
    with open(top_layer_file, 'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')

            pattern = items[-2].split('_')
            in_flag = 0
            pattern_name = ''
            for i in range(len(pattern)):  # 循环pattern
                prefix_pattern = pattern[i:]
                suffix_pattern = pattern[:i]
                loop_pattern = prefix_pattern + suffix_pattern
                s_loop_pattern = ''
                for j in loop_pattern:
                    s_loop_pattern += str(j) + '_'
                s_loop_pattern = s_loop_pattern[:-1]
                if s_loop_pattern in patterns_name.keys():
                    in_flag = 1
                    pattern_name = patterns_name[s_loop_pattern][0]
                    break
            if in_flag == 0:
                r_pattern = pattern[::-1]  # 反向
                for i in range(len(r_pattern)):  # 循环r_pattern
                    prefix_pattern = r_pattern[i:]
                    suffix_pattern = r_pattern[:i]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for j in loop_pattern:
                        s_loop_pattern += str(j) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in patterns_name.keys():
                        in_flag = 1
                        pattern_name = patterns_name[s_loop_pattern][0]
                        break

            out_top_layer_file.write(items[0] + '\t' +
                                     items[1] + '\t' +
                                     items[2] + '\t' +
                                     items[3] + '\t' +
                                     items[4] + '\t' +
                                     pattern_name + '\n')
    out_top_layer_file.close()

def reverse(sequence):
    base_map = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    new_sequence = ''
    for i in sequence[::-1]:
        new_sequence += base_map[i]
    return new_sequence

def readSDFile(file):
    read_blocks = {}
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            readname = items[0]
            start = int(items[1])
            end = int(items[2])
            strand = items[-1]
            if readname not in read_blocks.keys():
                read_blocks[readname] = [[start,end,strand]]
            else:
                read_blocks[readname].append([start,end,strand])
    return read_blocks

def readSequences(file):
    seqs= {}
    header = ''
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if line.startswith('>'):
                header = line[1:]
                seqs[header] = ''
            else:
                seqs[header] = line
    return seqs

def readMonSeqs(file):
    reads_mono = {}
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-2]
            if not line:
                break
            items = line.split('\t')
            readname = items[0]
            monomer_seq = items[1].split(' ')
            reads_mono[readname] = monomer_seq
    return reads_mono


def buildHORFile(patterns, patterns_info, base_sequences, reads_mono, read_blocks, outdir):
    out_hor_raw_file = outdir + '/out_hor.raw.merge.fa'
    out_hor_raw_file = open(out_hor_raw_file,'w')
    out_hor_normal_file = outdir + '/out_hor.normal.merge.fa'
    out_hor_normal_file = open(out_hor_normal_file,'w')
    for i in patterns.keys():
        pattern_name = patterns_info[i][0]
        pattern = i.split('_')
        database = patterns[i]
        # ([read_name,start,end,strand,pattern,repeat_number])
        for j in database:
            read_name = j[0]
            start = j[1]
            end = j[2]
            strand = j[3] # 更新增加strand
            repeat_number = j[5]
            monomer_sequence_item = reads_mono[read_name][start:end+1]
            # patternname.index start end pattern repeatnumber rawpattern
            monomer_sequence_item_str = ''
            for k in monomer_sequence_item:
                monomer_sequence_item_str += k + '_'
            monomer_sequence_item_str = monomer_sequence_item_str[:-1]
            out_hor_raw_file.write('>' + pattern_name  + '::' + read_name +'-' +
                                       str(read_blocks[read_name][start][0]) + '-' + str(read_blocks[read_name][end][1] + 1) +
                                   '::' + strand +
                                       ' nHOR-' + i + '::rHOR-' + monomer_sequence_item_str +
                                      ' repeat_number-' +  str(repeat_number) + '\n')
            out_hor_raw_file.write(base_sequences[read_name][read_blocks[read_name][start][0]:read_blocks[read_name][end][1] + 1] + '\n')

            out_hor_normal_file.write('>' + pattern_name + '::' + read_name +'-' +
                                   str(read_blocks[read_name][start][0]) + '-' + str(read_blocks[read_name][end][1] + 1) +
                                      '::' + strand +
                                   ' nHOR-' + i + '::rHOR-' + monomer_sequence_item_str +
                                      ' repeat_number-' +  str(repeat_number) + '\n')

            if len(pattern) == 1:
                # 考虑反链 '-' 链标准化变正
                normal_sequence = base_sequences[read_name][read_blocks[read_name][start][0]:read_blocks[read_name][end][1] + 1]
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
                        block_start = read_blocks[read_name][start + k[2]][0]
                        block_end = read_blocks[read_name][start + k[2]][1] + 1
                        normal_sequence += reverse(base_sequences[read_name][block_start:block_end]) # 每个block变反
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
                        block_start = read_blocks[read_name][start + k[2]][0]
                        block_end = read_blocks[read_name][start + k[2]][1]+1
                        normal_sequence += base_sequences[read_name][block_start:block_end]
                    out_hor_normal_file.write(normal_sequence + '\n')
    out_hor_raw_file.close()
    out_hor_normal_file.close()


def main():
    parser = argparse.ArgumentParser(description="Merge sample HOR results")
    parser.add_argument("-i", "--all_sample_dir", help="All sample result path, required", required=True)
    parser.add_argument("-s", "--sample_file", help="The file record all sample names, required", required=True)

    args = parser.parse_args()

    workdir = args.all_sample_dir
    sample_file = args.sample_file

    samples = []
    with open(sample_file, 'r') as sf:
        while True:
            line = sf.readline()[:-1]
            if not line:
                break
            ID = line.split('\t')[0]
            samples.append(ID)
    chr_sample_patterns = {}
    chr_sample_patterns_info = {}
    for i in samples:
        print(i)
        hicat_result_dir = workdir + '/' + i + '/reads_out/hicat_reads'
        chr_file = ['1','2','3','4','5','6','7','8','9','10',
                    '11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
        for j in chr_file:
            pattern_file = hicat_result_dir + '/' + j + '/out_final_hor.xls'
            patterns,patterns_info = readPattern(pattern_file)
            if j not in chr_sample_patterns.keys():
                chr_sample_patterns[j] = {}
                chr_sample_patterns[j][i] = patterns
                chr_sample_patterns_info[j] = {}
                chr_sample_patterns_info[j][i] = patterns_info
            else:
                chr_sample_patterns[j][i] = patterns
                chr_sample_patterns_info[j][i] = patterns_info

    # 每个chr进行标准化
    print('merge patterns')
    chr_merge_pattern = {}
    for chr in chr_sample_patterns.keys():
        print(chr)
        merge_patterns = {}
        merge_patterns_info = {}
        for sample in chr_sample_patterns[chr]:
            print(sample)
            patterns = chr_sample_patterns[chr][sample]
            patterns_info = chr_sample_patterns_info[chr][sample]
            for i in patterns.keys():
                pattern = i.split('_')
                in_flag = 0
                for j in range(len(pattern)):  # 循环pattern
                    prefix_pattern = pattern[j:]
                    suffix_pattern = pattern[:j]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for k in loop_pattern:
                        s_loop_pattern += str(k) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in merge_patterns.keys():
                        in_flag = 1
                        merge_patterns[s_loop_pattern].append([sample,patterns[i]])
                        merge_patterns_info[s_loop_pattern] += patterns_info[i][1]
                        break
                # 如果pattern 没有出现，则建立新的
                if in_flag == 0:
                    s_pattern = ''
                    for j in pattern:
                        s_pattern += str(j) + '_'
                    s_pattern = s_pattern[:-1]
                    merge_patterns[s_pattern] = [[sample,patterns[i]]]
                    merge_patterns_info[s_pattern] = patterns_info[i][1]
        chr_merge_pattern[chr] = [merge_patterns,merge_patterns_info]

    print('naming patterns')
    chr_pattern_name = {}
    for chr in chr_merge_pattern.keys():
        # 按照info排序pattern，命名MXLXX
        print(chr)
        merge_patterns = chr_merge_pattern[chr][0]
        merge_patterns_info = chr_merge_pattern[chr][1]
        sorted_patterns = sorted(merge_patterns_info.items(),key=lambda x:x[1],reverse=True)
        index = 1
        patterns_name_table = {}
        for i in sorted_patterns:
            pattern = i[0].split('_')
            patterns_name_table[i[0]] = 'M' + str(index) + 'L' + str(len(pattern))
            index += 1
        chr_pattern_name[chr] = patterns_name_table

        # 重新读文件替换
        for i in samples:
            hicat_result_dir = workdir + '/' + i + '/reads_out/hicat_reads'
            pattern_file = hicat_result_dir + '/' + chr + '/out_final_hor.xls'
            patterns,patterns_info = readPattern(pattern_file)
            patterns_name = {}
            for p in patterns.keys():
                patterns_name[p] = ''
                pattern = p.split('_')
                for j in range(len(pattern)):  # 循环pattern
                    prefix_pattern = pattern[j:]
                    suffix_pattern = pattern[:j]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for k in loop_pattern:
                        s_loop_pattern += str(k) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in patterns_name_table.keys():
                        patterns_name[p] = [patterns_name_table[s_loop_pattern],s_loop_pattern]
                        break
            out_pattern_file = hicat_result_dir + '/' + chr + '/out_final_hor.merge.xls'
            outPattern(patterns_info, patterns, patterns_name,out_pattern_file)

            # pattern fa需要重新矫正
            sd_path = workdir +'/'+i+'/reads_out/chr_data'+ '/chr'+chr+'.sd.xls'
            read_blocks = readSDFile(sd_path)
            reads_path = workdir +'/'+i+'/reads_out/chr_data'+ '/chr'+chr+'.fa'
            base_sequences = readSequences(reads_path)

            monomer_seq_file = workdir + '/' + i + '/reads_out/hicat_reads/'+chr+'/monseq.xls'
            reads_mono = readMonSeqs(monomer_seq_file)

            patterns, patterns_info = readPattern(out_pattern_file)
            buildHORFile(patterns, patterns_info,
                         base_sequences, reads_mono, read_blocks, hicat_result_dir + '/' + chr)

            all_layer_file = hicat_result_dir + '/' + chr + '/out_all_layer.xls'
            out_all_layer_file = hicat_result_dir + '/' + chr + '/out_all_layer.merge.xls'
            outAllLayer(all_layer_file,patterns_name,out_all_layer_file)
            top_layer_file = hicat_result_dir + '/' + chr + '/out_top_layer.xls'
            out_top_layer_file = hicat_result_dir + '/' + chr + '/out_top_layer.merge.xls'
            outTopLayer(top_layer_file,patterns_name,out_top_layer_file)

    # 输出总的数量表
    outfile = workdir + '/' + 'all.sample.HOR.repeatnumber.matrix.xls'
    sample_pattern_table = {}
    all_sample_pattern = set()
    for i in samples:
        print(i)
        sample_pattern_table[i] = {}
        hicat_result_dir = workdir + '/' + i + '/reads_out/hicat_reads'
        chr_file = ['1','2','3','4','5','6','7','8','9','10',
                    '11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
        for j in chr_file:
            pattern_file = hicat_result_dir + '/' + j + '/out_final_hor.merge.xls'
            patterns_info = {}
            with open(pattern_file, 'r') as pf:
                while True:
                    line = pf.readline()[:-1]
                    if not line:
                        break
                    key = line.split('\t')[0]
                    patterns_info[key] = int(line.split('\t')[-1])
                    line = pf.readline()[:-1]
            patterns, patterns_info = readPattern(pattern_file)
            for k in patterns_info.keys():
                all_sample_pattern.add(j+'_'+chr_pattern_name[j][k])
                sample_pattern_table[i][j+'_'+chr_pattern_name[j][k]] = patterns_info[k][1]
    for i in samples:
        for j in all_sample_pattern:
            if j not in sample_pattern_table[i].keys():
                sample_pattern_table[i][j] = 0
    sample_pattern = pd.DataFrame(sample_pattern_table)
    sample_pattern.to_csv(outfile,sep='\t')


if __name__ == '__main__':
    main()