import os
import edlib
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import argparse
os.environ['OPENBLAS_NUM_THREADS'] = '1'


def ed_distance(sequenceA, sequenceB):
    return (edlib.align(sequenceA,sequenceB))['editDistance'] / max(len(sequenceA), len(sequenceB))

def ed_distance_apply(target, col):
    return col.apply(ed_distance, args=(target,))

def ed_distance_apply_apply(monomers,data, region):
    return data[region[0]:region[1]]['seq'].apply(ed_distance_apply, args=(monomers['seq'],))

def reverse(sequence):
    base_map = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    new_sequence = ''
    for i in sequence[::-1]:
        new_sequence += base_map[i]
    return new_sequence

def calculateED(All_blocks,monomers, base_sequences,thread):
    names = []
    seq = []
    for i in All_blocks:
        names.append(i)
        items = i.split('@')
        readname = items[0]
        strand = items[3]
        if strand == '-':
            split_base_sequence = reverse(base_sequences[readname][int(items[1]):int(items[2])])
        else:
            split_base_sequence = base_sequences[readname][int(items[1]):int(items[2])]
        seq.append(split_base_sequence)
    data = pd.DataFrame({'name': names, 'seq': seq})
    names = []
    seq = []
    for i in monomers.keys():
        names.append(i)
        seq.append(monomers[i])
    monomers = pd.DataFrame({'name': names, 'seq': seq})

    binsize = int(len(data) / thread) + 1
    split_in = [[binsize * i, binsize * i + binsize] for i in range(thread)]
    res = Parallel(n_jobs=thread)(delayed(ed_distance_apply_apply)(monomers,data, i) for i in split_in)
    res = pd.concat(res)
    edit_distance_matrix = np.array(res)
    block_name_index = list(data['name'])
    monomer_name_index = list(monomers['name'])
    return edit_distance_matrix, block_name_index,monomer_name_index

def readLabel(label_file):
    chr_reads = {
        'chr1': [], 'chr2': [], 'chr3': [], 'chr4': [], 'chr5': [],
        'chr6': [], 'chr7': [], 'chr8': [], 'chr9': [], 'chr10': [],
        'chr11': [], 'chr12': [], 'chr13': [], 'chr14': [], 'chr15': [],
        'chr16': [], 'chr17': [], 'chr18': [], 'chr19': [], 'chr20': [],
        'chr21': [], 'chr22': [], 'chrX': [], 'chrY': []
    }
    with open(label_file,'r') as lf:
        while True:
            line = lf.readline()[:-1]
            if not line:
                break
            if line.startswith('readname'):
                continue
            items = line.split('\t')
            readname = items[0]
            chr = items[1]
            chr_reads[chr].append(readname)
    return chr_reads

def readSDBlock(file):
    read_blocks = {}
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            readname = items[0]
            strand = '+'
            if items[1][-1] == "'":
                strand = '-'
            start = int(items[2])
            end = int(items[3])
            if readname not in read_blocks.keys():
                read_blocks[readname] = [[start,end,strand]]
            else:
                read_blocks[readname].append([start,end,strand])
    return read_blocks

def readHiCATMonomers(file):
    monomer_seq = {} # 与出现次数超过2个的代表monomer比较
    monomer_number = {}
    monomer = ''
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if line.startswith('>'):
                monomer_info = line[1:].split(' ')
                monomer = monomer_info[0]
                number = int(monomer_info[1])
                monomer_number[monomer] = number
            else:
                monomer_seq[monomer] = line

    select_monomers = set() # 选择至少出现2次的monomer
    for i in monomer_number.keys():
        if monomer_number[i] > 1:
            select_monomers.add(i)
    monomers = {}
    for i in monomer_seq.keys():
        monomer = i
        if monomer in select_monomers:
            monomers[monomer] = monomer_seq[i]
    return monomers

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

def miningMonomerTDPattern(new_monomer_sequence, max_hor_len):
    new_monomer_sequence_index_left = []
    new_monomer_sequence_index_right = []
    for i in range(len(new_monomer_sequence)):
        new_monomer_sequence_index_left.append(i)
        new_monomer_sequence_index_right.append(i)
    ori_monomer_sequence = new_monomer_sequence
    top_layer = []
    all_layer = []
    all_layer_marker = set()
    start_d = 1

    while start_d < min(max_hor_len, len(new_monomer_sequence)):
        monomer_sequence = new_monomer_sequence
        # print('--------------------------')
        monomer_sequence_index_left = new_monomer_sequence_index_left
        monomer_sequence_index_right = new_monomer_sequence_index_right
        candidate_pattern = {}
        for i in range(len(monomer_sequence)):
            if str(monomer_sequence[i]) not in candidate_pattern.keys():
                candidate_pattern[str(monomer_sequence[i])] = [[i, i + 1, str(monomer_sequence[i])]]
            else:
                candidate_pattern[str(monomer_sequence[i])].append([i, i + 1, str(monomer_sequence[i])])

        d_database = []
        for i in candidate_pattern.keys():
            pattern_database = candidate_pattern[i]
            for j in range(len(pattern_database)):
                current = pattern_database[j]
                current_start = current[0]
                if j == 0:
                    continue
                candidate_database = pattern_database[:j]
                for k in candidate_database[::-1]:
                    pre_database = k
                    pre_start = pre_database[0]
                    d = current_start - pre_start
                    if d == start_d:
                        d_database.append([pre_database, current])
                    if d > start_d:
                        break
        if len(d_database) == 0:
            start_d += 1
            continue
        else:
            sorted_d_database = sorted(d_database, key=lambda x: x[0])
            # print(sorted_d_database)
            chain_list = []
            chain = []
            chain.append(sorted_d_database[0][0])
            index = 1
            while index < len(sorted_d_database):
                if sorted_d_database[index][0][0] - sorted_d_database[index - 1][0][0] != 1:
                    final_start = index - start_d
                    if final_start >= 0:
                        for i in range(start_d):
                            chain.append(sorted_d_database[final_start + i][1])

                    if int(len(chain) / start_d) > 1:
                        chain_list.append(chain[:int(len(chain) / start_d) * start_d])
                        # chain_list.append(chain)
                    chain = [sorted_d_database[index][0]]
                    index += 1
                else:
                    chain.append(sorted_d_database[index][0])
                    index += 1
                    if index >= len(sorted_d_database):
                        final_start = index - start_d
                        if final_start >= 0:
                            for i in range(start_d):
                                chain.append(sorted_d_database[final_start + i][1])
                        if int(len(chain) / start_d) > 1:
                            chain_list.append(chain[:int(len(chain) / start_d) * start_d])
                            # chain_list.append(chain)
                        chain = []
                        break

                    while sorted_d_database[index][0][0] - sorted_d_database[index - 1][0][0] == 1:
                        chain.append(sorted_d_database[index][0])
                        index += 1
                        if index >= len(sorted_d_database):
                            break
                    final_start = index - start_d
                    if final_start >= 0:
                        for i in range(start_d):
                            chain.append(sorted_d_database[final_start + i][1])
                    if int(len(chain) / start_d) > 1:
                        chain_list.append(chain[:int(len(chain) / start_d) * start_d])
                        # chain_list.append(chain)
                    chain = []
            if len(chain) != 0:
                final_start = index - start_d
                if final_start >= 0:
                    for i in range(start_d):
                        chain.append(sorted_d_database[final_start + i][1])
                if int(len(chain) / start_d) > 1:
                    chain_list.append(chain[:int(len(chain) / start_d) * start_d])
                    # chain_list.append(chain)
            # print(start_d)
            # print(chain_list)
            tmp_top_layer = []
            for i in range(len(chain_list)):
                relative_start = chain_list[i][0][0]
                relative_end = chain_list[i][-1][0]
                absolute_start = monomer_sequence_index_left[relative_start]
                absolute_end = monomer_sequence_index_right[relative_end]
                if len(top_layer) == 0:
                    relative_end = chain_list[i][int(len(chain_list[i]) / start_d) * start_d - 1][0]
                    absolute_end = monomer_sequence_index_right[relative_end]
                    TD_item = monomer_sequence[relative_start:relative_start + start_d]
                    TD_count = int((relative_end - relative_start + 1) / start_d)
                    TD_range_index = []
                    start_td_item = -1
                    end_td_item = -1
                    new_start_td_item = relative_start
                    new_end_td_item = relative_start + start_d - 1
                    pattern_number = 1
                    while True:
                        start_td_item = new_start_td_item
                        end_td_item = new_end_td_item
                        if end_td_item > relative_end:
                            break
                        adding_flag = 1
                        if monomer_sequence_index_right[end_td_item] in all_layer_marker:
                            adding_flag = 0
                        if adding_flag == 1:
                            TD_range_index.append(
                                [monomer_sequence_index_left[start_td_item], monomer_sequence_index_right[end_td_item],
                                 TD_item, pattern_number])
                            new_start_td_item = end_td_item + 1
                            new_end_td_item = end_td_item + start_d
                        else:
                            pattern_number += 1
                            new_start_td_item = start_td_item
                            new_end_td_item = end_td_item + start_d
                    td_pattern_item = [absolute_start, absolute_end, TD_item, TD_count, TD_range_index, relative_start,
                                       relative_end]
                    top_layer.append(td_pattern_item)
                    tmp_top_layer.append(td_pattern_item)
                    all_layer.append(td_pattern_item)

                    for j in range(absolute_start, absolute_end):
                        all_layer_marker.add(j)
                        # print('aaaaaaa')
                        # print(all_layer_marker)
                    continue

                non_overlap_flag = 0
                top_layer_state = []
                top_layer_ins = []
                complete_cover_flag = 0
                for j in top_layer:
                    top_layer_state.append(0)

                for j in range(len(top_layer)):
                    j_absolute_start = top_layer[j][0]
                    j_absolute_end = top_layer[j][1]
                    #            -------
                    # ---------           ------------
                    if absolute_end < j_absolute_start or j_absolute_end < absolute_start:
                        non_overlap_flag += 1
                        continue
                    #          j_ab_s-----------j_ab_e
                    #      ab_s------------------------ab_e
                    if absolute_start <= j_absolute_start and j_absolute_end <= absolute_end:
                        top_layer_state[j] = 1
                    #            ---------------                      ------------------
                    #                   -----------------                              -----------
                    if absolute_start <= j_absolute_end and j_absolute_start < absolute_start and absolute_end > j_absolute_end:
                        top_layer_state[j] = 2
                    #                     ----------------                  -----------------
                    #        ------------------                  ------------
                    if absolute_start < j_absolute_start and absolute_end >= j_absolute_start and absolute_end < j_absolute_end:
                        top_layer_state[j] = 3
                    #           ---------------------
                    #                 ---------
                    if absolute_start >= j_absolute_start and absolute_end <= j_absolute_end:
                        complete_cover_flag = 1

                if complete_cover_flag == 1:
                    continue

                if non_overlap_flag == len(top_layer):
                    TD_item = monomer_sequence[relative_start:relative_start + start_d]
                    TD_count = int((relative_end - relative_start + 1) / start_d)
                    TD_range_index = []
                    relative_end = chain_list[i][int(len(chain_list[i]) / start_d) * start_d - 1][0]
                    absolute_end = monomer_sequence_index_right[relative_end]
                    start_td_item = -1
                    end_td_item = -1
                    new_start_td_item = relative_start
                    new_end_td_item = relative_start + start_d - 1
                    pattern_number = 1
                    while True:
                        start_td_item = new_start_td_item
                        end_td_item = new_end_td_item
                        if end_td_item > relative_end:
                            break
                        adding_flag = 1

                        if monomer_sequence_index_right[end_td_item] in all_layer_marker:
                            adding_flag = 0
                        if adding_flag == 1:
                            TD_range_index.append(
                                [monomer_sequence_index_left[start_td_item], monomer_sequence_index_right[end_td_item],
                                 TD_item, pattern_number])
                            pattern_number = 1
                            new_start_td_item = end_td_item + 1
                            new_end_td_item = end_td_item + start_d
                        else:
                            pattern_number += 1
                            new_start_td_item = start_td_item
                            new_end_td_item = end_td_item + start_d

                    top_layer_ins = [absolute_start, absolute_end, TD_item, TD_count, TD_range_index, relative_start,
                                     relative_end]

                    top_layer.append(top_layer_ins)
                    tmp_top_layer.append(top_layer_ins)
                    all_layer.append(top_layer_ins)
                    # print(absolute_start)
                    # print(absolute_end)
                    for j in range(absolute_start, absolute_end):
                        all_layer_marker.add(j)
                        # print('kkkk')
                        # print(all_layer_marker)
                else:
                    processed_absolute_start = absolute_start
                    processed_absolute_end = absolute_end
                    processed_relative_start = relative_start
                    processed_relative_end = relative_end
                    abandon_flag = 0
                    while True:
                        top_layer_ins = []
                        for j in range(len(top_layer_state)):
                            #           -------------                  --------------
                            #                 -------------                          ------------
                            if top_layer_state[j] == 2:
                                j_absolute_end = top_layer[j][1]
                                if j_absolute_end + 1 > processed_absolute_start:
                                    processed_absolute_start = j_absolute_end + 1
                            #           -------------                  --------------
                            #     -------------                  ------
                            if top_layer_state[j] == 3:
                                j_absolute_start = top_layer[j][0]
                                if j_absolute_start - 1 < processed_absolute_end:
                                    processed_absolute_end = j_absolute_start - 1

                        for k in range(len(monomer_sequence)):
                            if monomer_sequence_index_left[k] == processed_absolute_start:
                                processed_relative_start = k
                            if monomer_sequence_index_right[k] == processed_absolute_end:
                                processed_relative_end = k

                        TD_count = int((processed_relative_end - processed_relative_start + 1) / start_d)
                        if TD_count > 1:
                            processed_relative_end = processed_relative_start + TD_count * start_d - 1
                            processed_absolute_end = monomer_sequence_index_right[processed_relative_end]
                            TD_range_index = []
                            TD_item = monomer_sequence[processed_relative_start:processed_relative_start + start_d]

                            start_td_item = -1
                            end_td_item = -1
                            new_start_td_item = processed_relative_start
                            new_end_td_item = processed_relative_start + start_d - 1
                            pattern_number = 1

                            while True:
                                start_td_item = new_start_td_item
                                end_td_item = new_end_td_item
                                if end_td_item > processed_relative_end:
                                    break
                                adding_flag = 1
                                if monomer_sequence_index_right[end_td_item] in all_layer_marker:
                                    adding_flag = 0
                                if adding_flag == 1:
                                    # print([monomer_sequence_index[start_td_item], monomer_sequence_index[end_td_item],
                                    #      TD_item,pattern_number])
                                    TD_range_index.append(
                                        [monomer_sequence_index_left[start_td_item],
                                         monomer_sequence_index_right[end_td_item],
                                         TD_item, pattern_number])

                                    pattern_number = 1
                                    new_start_td_item = end_td_item + 1
                                    new_end_td_item = end_td_item + start_d
                                else:
                                    pattern_number += 1
                                    new_start_td_item = start_td_item
                                    new_end_td_item = end_td_item + start_d

                            top_layer_ins = [processed_absolute_start, processed_absolute_end,
                                             TD_item, TD_count, TD_range_index,
                                             processed_relative_start, processed_relative_end]

                        if len(top_layer_ins) == 0:
                            abandon_flag = 1
                            break
                        if len(top_layer_ins[4]) == 0:
                            abandon_flag = 1
                            break
                        non_overlap_flag = 0
                        top_layer_state = []
                        complete_cover_flag = 0
                        for j in top_layer:
                            top_layer_state.append(0)

                        for j in range(len(top_layer)):
                            j_absolute_start = top_layer[j][0]
                            j_absolute_end = top_layer[j][1]
                            #            -------
                            # ---------           ------------
                            if processed_absolute_end < j_absolute_start or j_absolute_end < processed_absolute_start:
                                non_overlap_flag += 1
                                continue
                            #          j_ab_s-----------j_ab_e
                            #      ab_s------------------------ab_e
                            if processed_absolute_start <= j_absolute_start and j_absolute_end <= processed_absolute_end:
                                top_layer_state[j] = 1
                            #            ---------------                      ------------------
                            #                   -----------------                              -----------
                            if processed_absolute_start <= j_absolute_end and j_absolute_start < processed_absolute_start and processed_absolute_end > j_absolute_end:
                                top_layer_state[j] = 2
                            #                     ----------------                  -----------------
                            #        ------------------                  ------------
                            if processed_absolute_start < j_absolute_start and processed_absolute_end >= j_absolute_start and processed_absolute_end < j_absolute_end:
                                top_layer_state[j] = 3
                            #           ---------------------
                            #                 ---------
                            if processed_absolute_start >= j_absolute_start and processed_absolute_end <= j_absolute_end:
                                complete_cover_flag = 1
                        processed_flag = 0
                        for j in top_layer_state:
                            if j == 3:
                                processed_flag = 1
                            if j == 2:
                                processed_flag = 1
                        if processed_flag == 0:
                            break

                    if abandon_flag == 1:
                        continue

                    new_top_layer = []

                    for j in range(len(top_layer_state)):
                        if top_layer_state[j] == 1:
                            continue
                        else:
                            new_top_layer.append(top_layer[j])
                    new_top_layer.append(top_layer_ins)
                    tmp_top_layer.append(top_layer_ins)
                    all_layer.append(top_layer_ins)
                    for j in range(processed_absolute_start, processed_absolute_end):
                        all_layer_marker.add(j)
                        # print('llll')
                        # print(all_layer_marker)
                    top_layer = new_top_layer

            if len(tmp_top_layer) == 0:
                start_d += 1
                continue

            new_monomer_sequence_index_left = []
            new_monomer_sequence_index_right = []
            tmp_monomer_sequence_left = []
            tmp_monomer_sequence_right = []
            new_monomer_sequence = []
            for i in monomer_sequence:
                tmp_monomer_sequence_left.append(i)
                tmp_monomer_sequence_right.append(i)

            for i in tmp_top_layer:
                relative_start = i[5]
                relative_end = i[6]
                for j in range(relative_end - relative_start + 1):
                    tmp_monomer_sequence_left[relative_start + j] = '*'
                    tmp_monomer_sequence_right[relative_end - j] = '*'

                for j in range(start_d):
                    tmp_monomer_sequence_left[relative_start + j] = monomer_sequence[relative_start + j]
                    tmp_monomer_sequence_right[relative_end - j] = monomer_sequence[relative_end - j]

            for i in range(len(tmp_monomer_sequence_left)):
                if tmp_monomer_sequence_left[i] == '*':
                    continue
                else:
                    new_monomer_sequence.append(monomer_sequence[i])
                    new_monomer_sequence_index_left.append(monomer_sequence_index_left[i])

            for i in range(len(tmp_monomer_sequence_right)):
                if tmp_monomer_sequence_right[i] == '*':
                    continue
                else:
                    new_monomer_sequence_index_right.append(monomer_sequence_index_right[i])

            start_d = 1
    return new_monomer_sequence, top_layer, all_layer

def buildingHor(block_sequence,all_layer):
    final_HOR = {}
    for i in all_layer:
        start = i[0]
        end = i[1]
        start_block = block_sequence[start]
        strand = start_block[-1]
        pattern = i[2]
        pattern_range = i[4]
        # print(pattern)
        in_flag = 0
        # 循环前后缀，拼接loop pattern，判断是否存在，
        # 如果存在in_flag为1：
        # final_HOR[s_loop_pattern][0]添加新的start end
        # final_HOR[s_loop_pattern][1] += pattern_range pattern range叠加
        # 如果不存在in_flag为0：
        # 创建新的pattern
        # 修改：以首位block正反表示HOR正反
        if strand == '-':
            pattern = pattern[::-1]

        for j in range(len(pattern)): # 循环pattern
            prefix_pattern = pattern[j:]
            suffix_pattern = pattern[:j]
            loop_pattern = prefix_pattern + suffix_pattern
            s_loop_pattern = ''
            for k in loop_pattern:
                s_loop_pattern += str(k) + '_'
            s_loop_pattern = s_loop_pattern[:-1]
            if s_loop_pattern in final_HOR.keys():
                in_flag = 1
                final_HOR[s_loop_pattern][0].append([start, end])
                # pattern_range 添加+ / -，然后拼接
                new_pattern_range = []
                for k in pattern_range:
                    new_pattern_range.append([k[0],k[1], strand,k[2],k[3]])
                final_HOR[s_loop_pattern][1] += new_pattern_range
                break

        # 如果pattern没有出现，则建立新的
        if in_flag == 0:
            s_pattern = ''
            for j in pattern:
                s_pattern += str(j) + '_'
            s_pattern = s_pattern[:-1]
            new_pattern_range = []
            for j in pattern_range:
                new_pattern_range.append([j[0], j[1], strand, j[2], j[3]])
            final_HOR[s_pattern] = [[[start, end]], new_pattern_range]
    return final_HOR

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

def outLayer(read_blocks,reads_HOR,pattern_names, outdir):
    out_top_layer = outdir + '/out_top_layer.xls'
    out_top_layer = open(out_top_layer, 'w')
    out_all_layer = outdir + '/out_all_layer.xls'
    out_all_layer = open(out_all_layer, 'w')
    for r in reads_HOR.keys():
        top_layer = reads_HOR[r][1]
        all_layer = reads_HOR[r][2]
        # sort output
        top_layer_set = set()
        sorted_top_layer = sorted(top_layer,key=lambda x:x[0])
        layers = {}
        layers_info = {}
        sd_blocks = read_blocks[r]
        for j in sorted_top_layer:
            start = sd_blocks[j[0]][0]
            end = sd_blocks[j[1]][1]
            pattern = ''
            repeat_number = j[3]
            for p in j[2]:
                pattern += str(p) + '_'
            pattern = pattern[:-1]
            top_layer_set.add(pattern+'_'+str(start)+'_'+str(end))
            # 搜索名字
            in_flag = 0
            pattern_name = ''
            for k in range(len(j[2])):  # 循环pattern
                prefix_pattern = j[2][k:]
                suffix_pattern = j[2][:k]
                loop_pattern = prefix_pattern + suffix_pattern
                s_loop_pattern = ''
                for l in loop_pattern:
                    s_loop_pattern += str(l) + '_'
                s_loop_pattern = s_loop_pattern[:-1]
                if s_loop_pattern in pattern_names.keys():
                    in_flag = 1
                    pattern_name = pattern_names[s_loop_pattern]
                    break
            if in_flag == 0:
                r_pattern = j[2][::-1]  # 反向
                for k in range(len(r_pattern)):  # 循环r_pattern
                    prefix_pattern = r_pattern[k:]
                    suffix_pattern = r_pattern[:k]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for l in loop_pattern:
                        s_loop_pattern += str(l) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in pattern_names.keys():
                        in_flag = 1
                        pattern_name = pattern_names[s_loop_pattern]
                        break

            out_top_layer.write(r+'\t'+str(start)+'\t'+str(end)+'\t'+str(repeat_number)+
                                '\t'+str(pattern) + '\t' + pattern_name+'\n')

            layers[pattern +'_'+str(start)+'_'+str(end)] = []
            layers_info[pattern +'_'+str(start)+'_'+str(end)] = [start,end,repeat_number,pattern]

        # fix same start region of top and cov
        sorted_all_layer = sorted(all_layer, key=lambda x: x[0])
        for j in sorted_all_layer:
            start = sd_blocks[j[0]][0]
            end = sd_blocks[j[1]][1]
            pattern = ''
            repeat_number = j[3]
            for p in j[2]:
                pattern += str(p) + '_'
            pattern = pattern[:-1]
            if pattern + '_' + str(start) + '_' + str(end) in layers.keys():
                pass
            else:
                for l in layers.keys():
                    top_start = int(l.split('_')[-2])
                    top_end = int(l.split('_')[-1])
                    if start >= top_start and end <= top_end:
                        layers[l].append([start,end,repeat_number,pattern])

        for j in layers_info.keys():
            start = layers_info[j][0]
            end = layers_info[j][1]
            repeat_number = layers_info[j][2]
            pattern = layers_info[j][3]
            pattern_list = pattern.split('_')
            # 搜索名字
            pattern_name = ''
            in_flag = 0
            for k in range(len(pattern_list)):  # 循环pattern
                prefix_pattern = pattern_list[k:]
                suffix_pattern = pattern_list[:k]
                loop_pattern = prefix_pattern + suffix_pattern
                s_loop_pattern = ''
                for l in loop_pattern:
                    s_loop_pattern += str(l) + '_'
                s_loop_pattern = s_loop_pattern[:-1]
                if s_loop_pattern in pattern_names.keys():
                    in_flag = 1
                    pattern_name = pattern_names[s_loop_pattern]
                    break
            if in_flag == 0:
                r_pattern = pattern_list[::-1]  # 反向
                for k in range(len(r_pattern)):  # 循环r_pattern
                    prefix_pattern = r_pattern[k:]
                    suffix_pattern = r_pattern[:k]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for l in loop_pattern:
                        s_loop_pattern += str(l) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in pattern_names.keys():
                        in_flag = 1
                        pattern_name = pattern_names[s_loop_pattern]
                        break
            out_all_layer.write(r+'\t' +
                str(start) + '\t' + str(end) + '\t' + str(repeat_number) + '\t' + str(pattern) + '\t' + pattern_name +'\t'+ 'top' + '\n')
            for k in layers[j]:
                sub_start = k[0]
                sub_end = k[1]
                sub_repeat_number = k[2]
                sub_pattern = k[3]

                pattern_list = sub_pattern.split('_')
                # 搜索名字
                in_flag = 0
                pattern_name = ''
                for l in range(len(pattern_list)):  # 循环pattern
                    prefix_pattern = pattern_list[l:]
                    suffix_pattern = pattern_list[:l]
                    loop_pattern = prefix_pattern + suffix_pattern
                    s_loop_pattern = ''
                    for m in loop_pattern:
                        s_loop_pattern += str(m) + '_'
                    s_loop_pattern = s_loop_pattern[:-1]
                    if s_loop_pattern in pattern_names.keys():
                        in_flag = 1
                        pattern_name = pattern_names[s_loop_pattern]
                        break
                if in_flag == 0:
                    r_pattern = pattern_list[::-1]  # 反向
                    for l in range(len(r_pattern)):  # 循环r_pattern
                        prefix_pattern = r_pattern[l:]
                        suffix_pattern = r_pattern[:l]
                        loop_pattern = prefix_pattern + suffix_pattern
                        s_loop_pattern = ''
                        for m in loop_pattern:
                            s_loop_pattern += str(m) + '_'
                        s_loop_pattern = s_loop_pattern[:-1]
                        if s_loop_pattern in pattern_names.keys():
                            in_flag = 1
                            pattern_name = pattern_names[s_loop_pattern]
                            break

                out_all_layer.write(r+'\t' +
                    str(sub_start) + '\t' + str(sub_end) + '\t' + str(sub_repeat_number) + '\t' + str(
                        sub_pattern) +'\t' + pattern_name + '\t' + 'cover' + '\n')
    out_all_layer.close()
    out_top_layer.close()

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

def buildHORFile(patterns, patterns_info, base_sequences, reads_mono, read_blocks, outdir):
    out_hor_raw_file = outdir + '/out_hor.raw.fa'
    out_hor_raw_file = open(out_hor_raw_file,'w')
    out_hor_normal_file = outdir + '/out_hor.normal.fa'
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

def buildHORForReads(read_blocks,base_sequences,monomer_seq_file,outdir,max_hor_len = 40):
    # HOR的正反性，尤其是CEN1上，原版HiCAT更新正负的处理，归类,、
    # Y染色体HOR太长，存在一个reads不够至少两个unit的情况
    # 未注释的block可以额外计算相互之间ed重新注释
    # -----
    # 需要将ref中的区间变换成绝对区间，方便与结果进行比较
    # 评估，
    # 1、评估主要pattern的情况，
    # 2、出错的区间在什么位置，是否在activate HOR 区间
    # 3、在第一阶段中错分的reads中的HOR注释情况
    # 考虑模拟不同长度的测序数据 30-45k 50-65k 70-85k
    # 其他可能尝试缩小工具应用范畴，放弃 非activate HOR区间，仅计算activate HOR
    print('start')
    reads_mono = readMonSeqs(monomer_seq_file)
    reads_HOR = {}
    for r in reads_mono.keys():
        monomer_sequence = reads_mono[r]
        block_sequence = read_blocks[r]
        new_monomer_sequence, top_layer, all_layer = miningMonomerTDPattern(monomer_sequence, max_hor_len)
        HORs = buildingHor(block_sequence,all_layer)
        filter_HORs = {}
        # 确保没有重叠的
        for i in HORs.keys():
            pattern = i
            database = HORs[i][1]
            filter_HOR = []
            init_sequences = []
            for j in monomer_sequence:
                init_sequences.append(0)
            sort_database = sorted(database, key=lambda x: x[1] - x[0])
            for j in sort_database:
                start = j[0]
                end = j[1]
                jump = 0
                for k in range(start, end + 1):
                    if init_sequences[k] == 1:
                        jump = 1
                        break
                if jump == 1:
                    continue
                else:
                    for k in range(start, end + 1):
                        init_sequences[k] = 1
                    filter_HOR.append(j)
            if len(filter_HOR) == 0:
                continue
            else:
                filter_HORs[pattern] = sorted(filter_HOR, key=lambda x: x[0])

        reads_HOR[r] = [filter_HORs, top_layer, all_layer]

    # 统计HOR，记录总数，记录每根reads上的数量
    # 1_2
    # P(read_name, start, end, strand, pattern, repeat_number):
    final_HORs = {}
    for r in reads_HOR.keys():
        for i in reads_HOR[r][0].keys():
            pattern = i.split('_')
            pattern_database = reads_HOR[r][0][i]
            in_flag = 0
            for j in range(len(pattern)):
                prefix_pattern = pattern[j:]
                suffix_pattern = pattern[:j]
                loop_pattern = prefix_pattern + suffix_pattern
                s_loop_pattern = ''
                for k in loop_pattern:
                    s_loop_pattern += str(k) + '_'
                s_loop_pattern = s_loop_pattern[:-1]
                if s_loop_pattern in final_HORs.keys():
                    in_flag = 1
                    for k in pattern_database:
                        final_HORs[s_loop_pattern].append([r, k[0], k[1], k[2], k[3], k[4]])
                    break
            if in_flag == 0:
                s_pattern = ''
                for j in pattern:
                    s_pattern += str(j) + '_'
                s_pattern = s_pattern[:-1]
                final_HORs[s_pattern] = []
                for j in pattern_database:
                    final_HORs[s_pattern].append([r, j[0], j[1], j[2], j[3], j[4]])

    out_pattern_file = outdir + '/' + 'out_final_hor.xls'
    out_pattern_file = open(out_pattern_file, 'w')
    sorted_HOR = sorted(final_HORs.items(), key=lambda x: len(x[1]), reverse=True)
    new_index = 1
    pattern_names = {}
    for h in sorted_HOR:
        i = h[0]
        pattern = i.split('_')
        # 独立重新命名 R index + L
        repeat_number = 0
        for j in final_HORs[i]:
            repeat_number += j[5]
        out_pattern_file.write('R' + str(new_index) + 'L' + str(len(pattern)) + '\t' + i + '\t' + str(repeat_number) + '\n')
        pattern_names[i] = 'R' + str(new_index) + 'L' + str(len(pattern))
        out_pattern_file.write('P(read_name,start,end,strand,pattern,repeat_number):')
        for j in final_HORs[i]:
            out_pattern_file.write('\t' + str(j[0]) + ',' + str(j[1]) + ',' + str(j[2]) + ',' + str(j[3]))
            pattern = ''
            for l in j[4]:
                pattern += str(l) + '_'
            pattern = pattern[:-1]
            out_pattern_file.write(',' + pattern)
            out_pattern_file.write(',' + str(j[5]))
        out_pattern_file.write('\n')
        new_index += 1
    out_pattern_file.close()

    # 输出按照reads reads_HOR[r] = filter_HORs
    outLayer(read_blocks,reads_HOR,pattern_names, outdir)

    # 输出HOR fa
    pattern_file = outdir + '/' + 'out_final_hor.xls'
    patterns, patterns_info = readPattern(pattern_file)
    buildHORFile(patterns, patterns_info, base_sequences, reads_mono, read_blocks, outdir)


def main():
    parser = argparse.ArgumentParser(description="Annotation HOR for centromere reads")
    parser.add_argument("-i", "--hicat_reads_dir", help="HiCAT reads output path, required", required=True)
    parser.add_argument("-r", "--hicat_ref_dir", help="HiCAT reference path, required", required=True)

    parser.add_argument("-ms", "--min_similarity",
                        help="The lower bound for similarity threshold which used to define monomer in reads, default is 0.9",
                        type=float, default=0.9, required=False)
    parser.add_argument("-mh", "--max_hor_len",
                        help="An upper bound for the length of the tandem repeat unit by default 40 monomers for improving efficiency",
                        type=int, default=40, required=False)
    parser.add_argument("-th", "--thread", help="The number of threads, default is 1", type=int, default=1,
                        required=False)

    args = parser.parse_args()

    hicat_reads_dir = args.hicat_reads_dir
    hicat_ref_dir = args.hicat_ref_dir
    thread = args.thread
    min_similarity = args.min_similarity
    max_hor_len = args.max_hor_len

    merge_cen_fasta_file = hicat_reads_dir + '/merge.fa'

    # run stringDecomposer
    script_path = os.path.split(os.path.realpath(__file__))[0]
    monomer_template_file = script_path + '/AlphaSat.fa'
    output_SD_dir = hicat_reads_dir + '/SD'

    if not os.path.exists(output_SD_dir):
        os.mkdir(output_SD_dir)

    cmd = 'python ' + \
          script_path + '/stringdecomposer/bin/stringdecomposer' + ' ' + \
          merge_cen_fasta_file + ' ' + monomer_template_file + ' -o ' + output_SD_dir + ' -t ' + str(thread)

    print('Run stringdecomposer\n')
    print(cmd)
    os.system(cmd)
    sd_file = output_SD_dir + '/final_decomposition.tsv'
    label_file = hicat_reads_dir + '/merge.label.xls'
    chr_reads = readLabel(label_file)
    read_blocks = readSDBlock(sd_file)
    out_result_dir = hicat_reads_dir + '/chr_data'
    for i in chr_reads.keys():
        outfile = out_result_dir + '/' + i + '.sd.xls'
        outfile = open(outfile, 'w')
        for j in chr_reads[i]:
            for k in read_blocks[j]:
                outfile.write(str(j) + '\t' + str(k[0]) + '\t' + str(k[1]) + '\t' + str(k[2]) + '\n')
        outfile.close()

    outdir = hicat_reads_dir + '/hicat_reads'

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    print('transform read block sequences to monomer sequences and mining HOR')
    for i in chr_reads.keys():
        print(i)
        hicat_monomer_file = hicat_ref_dir + '/' + i[3:] + '/out/out_monomer_represent.fa'
        monomers = readHiCATMonomers(hicat_monomer_file)
        sd_path = out_result_dir + '/' + i + '.sd.xls'
        if not os.path.exists(sd_path):
            continue
        read_blocks = readSDFile(sd_path)
        reads_path = out_result_dir + '/' + i + '.fa'
        base_sequences = readSequences(reads_path)
        if len(base_sequences) == 0:
            continue
        out_chr_dir = outdir + '/' + i[3:]
        if not os.path.exists(out_chr_dir):
            os.mkdir(out_chr_dir)
        # 得到全部block
        All_blocks = []
        for j in read_blocks.keys():
            for k in read_blocks[j]:
                All_blocks.append(j + '@' + str(k[0]) + '@' + str(k[1]) + '@' + str(k[2]))
        print('thread: ' + str(thread))
        edit_distance_matrix, block_name_index, monomer_name_index = calculateED(All_blocks,
                                                                                 monomers, base_sequences,
                                                                                 thread)

        eddis = pd.DataFrame(edit_distance_matrix, index=block_name_index, columns=monomer_name_index)
        #eddis.to_csv(out_chr_dir + '/eddis.xls', sep='\t')

        uniq_monomer_ID = -1
        block2monomers = []
        for j in range(len(block_name_index)):
            block_name = block_name_index[j]
            minindex = np.argmin(edit_distance_matrix[j])
            min_monomer = monomer_name_index[minindex]
            max_identity = 1 - edit_distance_matrix[j][minindex]
            if max_identity < min_similarity:
                block2monomers.append([block_name, uniq_monomer_ID])
                uniq_monomer_ID = uniq_monomer_ID - 1
            else:
                block2monomers.append([block_name, min_monomer])
        out_read_blocks = {}
        for j in block2monomers:
            block_name = j[0]
            min_monomer = j[1]
            block_info = block_name.split('@')
            readname = block_info[0]
            if readname not in out_read_blocks.keys():
                out_read_blocks[readname] = [min_monomer]
            else:
                out_read_blocks[readname].append(min_monomer)

        outmonseqfile = open(out_chr_dir + '/monseq.xls', 'w')
        for j in out_read_blocks.keys():
            outmonseqfile.write(j + '\t')
            for k in out_read_blocks[j]:
                outmonseqfile.write(str(k) + ' ')
            outmonseqfile.write('\n')
        outmonseqfile.close()
        buildHORForReads(read_blocks, base_sequences,out_chr_dir + '/monseq.xls', out_chr_dir,
                         max_hor_len = max_hor_len)





if __name__ == '__main__':
    main()
