import pickle

class Predictor:
    def __init__(self, model_file):
        with open(model_file, 'rb') as infile:
            self.loaded = pickle.load(infile)
        self.model = self.loaded['svm']
        self.pca = self.loaded['pca']

    def predict(self, features):
        features_pca = self.pca.transform(features)
        y = self.model.predict(features_pca)
        return y

    def outPredict(self):
        pass
