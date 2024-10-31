import pandas as pd
import argparse
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-fi", "--fai", required=True)
    parser.add_argument("-w", "--workdir", required=True)

    args = parser.parse_args()

    fai = args.fai
    workdir = args.workdir

    outfile = workdir + '/lastz_regions.bed'

    outfile = open(outfile,'w')
    chrs = []
    with open(fai, 'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            if line.startswith('chr'):
                chrs.append(items[0])
    # 连续区间间隔5K以内合并，最终区间超过10K保留
    for i in chrs:
        if i == 'chrM':
            continue
        lastz_file = workdir + '/' + i + '.xls'
        print(lastz_file)
        beds = []
        with open(lastz_file,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                if line.startswith('#'):
                    continue
                items = line.split('\t')
                start = int(items[4])
                end = int(items[5])
                beds.append([start,end])
        beds = sorted(beds,key=lambda x:x[0])
        continue_region = []
        if len(beds) == 0:
            continue
        init = beds[0]
        for j in range(len(beds) - 1):
            if beds[j + 1][0] - init[1] < 5000:
                init[1] = beds[j + 1][1]
            else:
                if init[1] - init[0] > 10000:
                    continue_region.append(init)
                init = beds[j + 1]
        if init[1] - init[0] > 10000:
            continue_region.append(init)
        for j in continue_region:
            outfile.write(i + '\t' + str(j[0]) + '\t' + str(j[1]) + '\t' + str(j[1] - j[0]) + '\n')
    outfile.close()




if __name__ == '__main__':
    main()