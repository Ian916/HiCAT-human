import pandas as pd


class FeatureBuilder:
    def __init__(self, kmers_file):
        self.reads = {}
        self.features = []
        with open(kmers_file, 'r') as kf:
            while True:
                line = kf.readline()[:-1]
                if not line:
                    break
                if line not in self.features:
                    self.features.append(line)

    def __reverse(self, sequence):
        base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        new_sequence = ''
        for i in sequence[::-1]:
            new_sequence += base_map[i]
        return new_sequence

    def loadReads(self, fa_file):
        self.reads = {}
        header = ''
        with open(fa_file, 'r') as rf:
            while True:
                line = rf.readline()[:-1]
                if not line:
                    break
                if line.startswith('>'):
                    header = line[1:]
                    self.reads[header] = ''
                else:
                    self.reads[header] += line
        print('reads number')
        print(len(self.reads.keys()))

    def generateFeatures(self):
        reads_features = {}
        count = 0
        for i in self.reads.keys():
            if count % 100 == 0:
                print('processed ' + str(count) + ' reads')
            reads_features[i] = {}
            for j in self.features:
                reads_features[i][j] = 0
            read = self.reads[i]
            for j in range(len(read) - 5):
                kmer = read[j:j + 5]
                r_kmer = self.__reverse(kmer)
                if kmer in reads_features[i].keys():
                    reads_features[i][kmer] += 1
                if r_kmer in reads_features[i].keys():
                    reads_features[i][r_kmer] += 1
            # 根据长度标准化
            read_len = len(read)
            for j in reads_features[i].keys():
                reads_features[i][j] = reads_features[i][j] / read_len
            count += 1
        features = pd.DataFrame(reads_features).T
        return features
