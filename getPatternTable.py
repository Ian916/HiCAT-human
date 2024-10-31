import argparse

def main():
    parser = argparse.ArgumentParser(description="get pattern summary")
    parser.add_argument("-i", "--all_sample_dir", help="All sample result path, required", required=True)

    args = parser.parse_args()

    workdir = args.all_sample_dir
    file = workdir + '/pattern.summary.xls'

    pattern_table = {}
    with open(file, 'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            pattern = items[1] + '_' + items[2]
            mon_pattern = items[3]
            pattern_table[pattern] = mon_pattern
    outfile = workdir + '/pattern.table.xls'
    outfile = open(outfile, 'w')
    for i in pattern_table.keys():
        outfile.write(i+'\t'+pattern_table[i]+'\n')
    outfile.close()


if __name__ == '__main__':
    main()
