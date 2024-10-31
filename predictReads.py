from buildPredictor import Predictor
from buildFeatures import FeatureBuilder
import argparse
import pandas as pd

def reverse(sequence):
    base_map = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    new_sequence = ''
    for i in sequence[::-1]:
        new_sequence += base_map[i]
    return new_sequence

# 全部load进来，考虑不内存，预测
def main():
    # 两步，先计算feature，然后预测,两步model
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-r", "--reads_file")
    parser.add_argument("-m1", "--model1_file")
    parser.add_argument("-m2", "--model2_file")
    parser.add_argument("-f", "--featurename_file")
    parser.add_argument("-p", "--prefix")
    parser.add_argument("-o", "--outdir")

    args = parser.parse_args()

    reads_file = args.reads_file
    model1_file = args.model1_file
    model2_file = args.model2_file
    featurename_file = args.featurename_file
    outdir = args.outdir
    prefix = args.prefix

    fb = FeatureBuilder(featurename_file)
    fb.loadReads(reads_file)
    features = fb.generateFeatures()
    features.to_csv(outdir+'/'+prefix+'.feature.xls', sep='\t')
    feature_list = pd.read_csv(outdir+'/'+prefix+'.feature.xls', sep='\t', index_col=0)
    read_name = feature_list.index.tolist()
    predictor = Predictor(model1_file)
    label = predictor.predict(feature_list)
    cen_readname = []
    for i in range(len(label)):
        if label[i] == 'CEN':
            cen_readname.append(read_name[i])
    if len(cen_readname) != 0:
        cen_feature_list = feature_list[feature_list.index.isin(cen_readname)]
        cen_readname = cen_feature_list.index.tolist()
        predictor = Predictor(model2_file)
        label = predictor.predict(cen_feature_list)
        outfile = outdir + '/' + prefix + '.label.xls'
        out_cen_readfile = outdir + '/' + prefix + '.cen.fa'
        out_cen_readfile = open(out_cen_readfile,'w')
        outfile = open(outfile,'w')
        outfile.write('readname\tlabel\n')
        for i in range(len(label)):
            outfile.write(cen_readname[i]+'\t' + label[i]+'\n')
            out_cen_readfile.write('>'+cen_readname[i]+'\n')
            out_cen_readfile.write(fb.reads[cen_readname[i]] +'\n')
        out_cen_readfile.close()
        outfile.close()
    else:
        outfile = outdir + '/' + prefix + '.label.xls'
        out_cen_readfile = outdir + '/' + prefix + '.cen.fa'
        out_cen_readfile = open(out_cen_readfile, 'w')
        outfile = open(outfile, 'w')
        outfile.write('readname\tlabel\n')
        out_cen_readfile.close()
        outfile.close()


if __name__ == '__main__':
    main()