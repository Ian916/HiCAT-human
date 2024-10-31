import argparse

import os

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-fi", "--fai", required=True)
    parser.add_argument("-w", "--workdir", required=True)
    parser.add_argument("-a", "--alpha", required=True)

    args = parser.parse_args()

    ref = args.ref
    fai = args.fai
    workdir = args.workdir
    alpha = args.alpha

    chrs = []
    with open(fai,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            if line.startswith('chr'):
                chrs.append(items[0])

    for i in chrs:
        if i == 'chrM':
            continue
        cmd = 'samtools faidx ' + ref + ' ' + i + ' > ' + workdir + '/' + i + '.fa'
        os.system(cmd)
        cmd = 'lastz ' + workdir + '/' + i + '.fa' + ' ' + alpha + ' ' + \
              '--format=general:score,name1,strand1,size1,start1,' \
              'end1,name2,strand2,identity,' \
              'length1,align1' + ' > ' + workdir + '/' + i + '.xls'
        os.system(cmd)


if __name__ == '__main__':
    main()