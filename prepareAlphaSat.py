# 根据连续区间获取拼接的CENSat seq
import os
import argparse

def readFa(file):
    seq = ''
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if line.startswith('>'):
                continue
            seq += line
    return seq

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-f", "--file", required=True)
    parser.add_argument("-w", "--workdir", required=True)

    args = parser.parse_args()

    file = args.file
    ref = args.ref
    workdir = args.workdir
    outdir = workdir + '/censeq'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    chr_seq = {}
    print('start')
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            chr = items[0]
            start = int(items[1])
            end = int(items[2])
            length = int(items[3])
            if chr not in chr_seq.keys():
                chr_seq[chr] = [[start,end,length]]
            else:
                chr_seq[chr].append([start,end,length])

    print('chr_seq')
    for i in chr_seq.keys():
        print(i)
        out_chr = outdir + '/' + i
        if not os.path.exists(out_chr):
            os.mkdir(out_chr)
        count = 1
        file_list = []
        for j in chr_seq[i]:
            cmd = 'samtools faidx ' + ref + ' ' + i+':'+str(j[0])+'-'+str(j[1]) +' > ' + out_chr + '/' + str(count) + '.fa'
            file_list.append(out_chr + '/' + str(count) + '.fa')
            count += 1
            os.system(cmd)
        cen = ''
        regions = []
        start = 0
        for j in file_list:
            i_cen = readFa(j)
            cen += i_cen
            end = len(i_cen) + start
            regions.append([start,end]) # [start,end)
            start = end

        out_seq_file = out_chr + '/' + i+'.cen.fa'
        out_seq_file = open(out_seq_file,'w')
        out_seq_file.write('>' + i +'\n')
        out_seq_file.write(cen+'\n')
        out_seq_file.close()
        out_region_file = out_chr + '/' + i + '.region.xls'
        out_region_file = open(out_region_file,'w')
        for j in regions:
            out_region_file.write(str(j[0]) + '\t' + str(j[1])+'\n')
        out_region_file.close()

if __name__ == '__main__':
    main()