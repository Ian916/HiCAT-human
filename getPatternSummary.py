import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="get pattern summary")
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
            ID = line.strip().split('\t')[0]
            samples.append(ID)
    outfile = workdir + '/pattern.summary.xls'
    outfile = open(outfile, 'w')
    for i in samples:
        print(i)
        chr_dir = workdir + '/' + i + '/reads_out/hicat_reads'
        chrs = os.listdir(chr_dir)
        for j in chrs:
            print(j)
            file = chr_dir + '/' + j + '/out_hor.normal.merge.fa'
            with open(file, 'r') as f:
                while True:
                    line = f.readline()[:-1]
                    if not line:
                        break
                    if line.startswith('>'):
                        items = line.split(' ')
                        info = items[0][1:].split('::')[0]
                        pattern = items[1].split('::')
                        nHOR = pattern[0][5:]
                        rHOR = pattern[1][5:]
                        number = items[2][14:]
                        outfile.write(i+'\t'+j+'\t'+info+'\t'+
                                      nHOR+'\t'+rHOR+'\t'+number+'\n')
    outfile.close()


if __name__ == '__main__':
    main()
