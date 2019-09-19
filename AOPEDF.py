# -*- coding: utf-8 -*-

import numpy as np
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.metrics import average_precision_score
import scipy.io as sio
from gcforest.gcforest import GCForest
import matplotlib.pyplot as plt
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler,StandardScaler


def load_csv_data(filename):
    data = []
    datafile = open(filename)
    for line in datafile:
        fields = line.strip().split('\t')
        data.append([float(field) for field in fields])
    data = np.array(data)
    return data

def get_toy_config():
    config = {}
    ca_config = {}
    ca_config["random_state"] = 0
    ca_config["max_layers"] = 100
    ca_config["early_stopping_rounds"] = 3
    ca_config["n_classes"] = 2
    ca_config["estimators"] = []
    ca_config["estimators"].append(
            {"n_folds": 5, "type": "XGBClassifier", "n_estimators": 500, "max_depth": 5,
             "objective": "multi:softprob", "silent": True, "nthread": -1, "learning_rate": 0.1,"num_class": 2} )
    ca_config["estimators"].append(
            {"n_folds": 5, "type": "XGBClassifier", "n_estimators": 500, "max_depth": 5,
             "objective": "multi:softprob", "silent": True, "nthread": -1, "learning_rate": 0.1,"num_class": 2} )
    ca_config["estimators"].append({"n_folds": 5, "type": "RandomForestClassifier", "n_estimators": 500, "max_depth": None, "n_jobs": -1})
    ca_config["estimators"].append({"n_folds": 5, "type": "RandomForestClassifier", "n_estimators": 500, "max_depth": None, "n_jobs": -1})
    ca_config["estimators"].append({"n_folds": 5, "type": "ExtraTreesClassifier", "n_estimators": 500, "max_depth": None, "n_jobs": -1})
    ca_config["estimators"].append({"n_folds": 5, "type": "ExtraTreesClassifier", "n_estimators": 500, "max_depth": None, "n_jobs": -1})
    config["cascade"] = ca_config
    return config

#Feature extracted from 15 networks by AROPE
drugFeature=np.loadtxt('drugFeature.txt')
proteinFeature=np.loadtxt('proteinFeature.txt')
interaction=np.loadtxt('dataset/Networks/drugProtein.txt')


positive_feature=[]
negative_feature=[]
alldata=[]
for i in range(np.shape(interaction)[0]):
    for j in range(np.shape(interaction)[1]):
        temp=np.append(drugFeature[i],proteinFeature[j])
        if int(interaction[i][j])==1:
            positive_feature.append(temp)
        elif int(interaction[i][j])==0:
            negative_feature.append(temp)
negative_sample_index = np.random.choice(np.arange(len(negative_feature)),size=len(positive_feature),replace=False)
negative_sample_feature=[]
for i in negative_sample_index:
    negative_sample_feature.append(negative_feature[i])
feature=np.vstack((positive_feature,negative_sample_feature))
label1=np.ones((len(positive_feature),1))
label0=np.zeros((len(negative_sample_feature),1))
label=np.vstack((label1,label0))


rs = np.random.randint(0, 1000, 1)[0]
kf = StratifiedKFold(label[:,0], n_folds=5, shuffle=True, random_state=rs)

test_auc_fold = []
test_aupr_fold = []
for train_index, test_index in kf:
    Xtrain, Xtest = feature[train_index], feature[test_index]
    Ytrain, Ytest = label[train_index], label[test_index]

    config = get_toy_config()
    rf=GCForest(config)
    Ytrain=Ytrain.flatten()
    rf.fit_transform(Xtrain, Ytrain)

    # deep forest
    predict_y = rf.predict(Xtest)
    acc = accuracy_score(Ytest, predict_y)
    print("Test Accuracy of GcForest = {:.2f} %".format(acc * 100))
    prob_predict_y = rf.predict_proba(Xtest)  # Give a result with probability values，the probability sum is 1 
    predictions_validation = prob_predict_y[:, 1]
    fpr, tpr, _ = roc_curve(Ytest, predictions_validation)
    roc_auc = auc(fpr, tpr)
    aupr = average_precision_score(Ytest, predictions_validation)
    print(roc_auc)
    print(aupr)
    test_auc_fold.append(roc_auc)
    test_aupr_fold.append(aupr)
    plt.figure()
    plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
mean_auc=np.mean(test_auc_fold)
mean_pr=np.mean(test_aupr_fold)
print('mean auc aupr', mean_auc, mean_pr)
