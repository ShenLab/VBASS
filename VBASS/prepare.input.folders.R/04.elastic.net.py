# Load the LogisticRegression classifier
# Note, use CV for cross-validation as requested in the question
from sklearn.linear_model import LogisticRegressionCV

# Load some other sklearn functions
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

# Import other libraries
import pandas as pd, numpy as np
import sys

X = pd.read_csv(f'{sys.argv[1]}/score.csv', index_col=0)

# do SVD on X before training
Xu, Xs, Xvt = np.linalg.svd(X)
def calc_k(S, percentge):
  '''
  identify the minimum k value to reconstruct
  '''
  k = 0
  total = sum(np.square(S))
  svss = 0 #singular values square sum
  for i in range(np.shape(S)[0]):
      svss += np.square(S[i])
      if (svss/total) >= percentge:
          k = i+1
          break
  return k

k = calc_k(Xs, 0.95) # get the number of k to reconstruct 0.9 square sum

def buildSD(S, k):
  '''
  reconstruct k singular value diag matrix
  '''
  SD = np.eye(k) * S[:k]
  return SD

def buildSD_1(S, k):
  '''
  reconstruct k singular value diag matrix
  '''
  SD = np.eye(k) * 1/S[:k]
  return SD


Xu_new = pd.DataFrame(Xu[:len(Xu), :k], index=X.index)
Xvt_new = Xvt[:k, :len(Xvt)]
Xs_new = buildSD(Xs, k)
Xs_new_1 = buildSD_1(Xs, k)
X_new = pd.DataFrame(np.dot(np.dot(Xu_new, Xs_new), Xvt_new), index=X.index)

labels = pd.read_csv(f'{sys.argv[1]}/label.csv')
y = labels['conditions']
y[y=='SPARK.risk'] = 0
y[y=='control.LGD'] = 1
y = y.astype(int)

X_train, X_test, y_train, y_test = train_test_split(X_new, y, stratify=y, random_state=0)

elastic_net_classifier = LogisticRegressionCV(cv=5, penalty='elasticnet', l1_ratios=np.linspace(0, 1, 21), solver='saga', tol=0.0001, max_iter=5000, n_jobs=64)

elastic_net_classifier.fit(X_train, y_train)

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=0)

with open(f'{sys.argv[1]}/elastic.net.report.txt', 'w') as outputfile:
  outputfile.write('Training Summary:\n')
  outputfile.write(classification_report(y_train, elastic_net_classifier.predict(X_train)))
  outputfile.write('\nTesting Summary:\n')
  outputfile.write(classification_report(y_test, elastic_net_classifier.predict(X_test)))
  outputfile.write('\n')
  outputfile.write("Training accuracy: {}\n\n".format(elastic_net_classifier.score(X_train, y_train)))
  outputfile.write("Testing accuracy: {}".format(elastic_net_classifier.score(X_test, y_test)))


import pickle
with open(f'{sys.argv[1]}/elasticnet.pkl', 'wb') as outputfile:
  pickle.dump(elastic_net_classifier, outputfile)
  
prediction=labels
prediction['is.train']=np.isin(prediction['genes'], X_train.index.values)
prediction['prediction']=elastic_net_classifier.predict(X_new)
prediction['prediction_prob']=elastic_net_classifier.predict_proba(X_new)[:,0]
prediction.to_csv(f'{sys.argv[1]}/prediction.csv')

all_X = pd.read_csv(f'{sys.argv[1]}/all.scores.csv', index_col=0)
all_prediction = pd.DataFrame(elastic_net_classifier.predict_proba(all_X), index=all_X.index)
all_prediction.to_csv(f'{sys.argv[1]}/all.prediction.csv')
# PC_loadings = pd.read_csv(f'{sys.argv[1]}/PC.loadings.csv', index_col=0)
# Xs_new_Xvt_new = np.dot(Xs_new, Xvt_new)
# Xs_new_Xvt_new_1 = np.dot(Xvt_new.T, Xs_new_1)
# elastic_net_loadings = pd.DataFrame(np.dot(Xs_new_Xvt_new_1, elastic_net_classifier.coef_.T), index=X.columns.values)
elastic_net_loadings = pd.DataFrame(elastic_net_classifier.coef_.T, index=X.columns.values)
elastic_net_loadings.to_csv(f'{sys.argv[1]}/elastic.net.loadings.csv')
elastic_net_intercept = pd.DataFrame(elastic_net_classifier.intercept_)
elastic_net_intercept.to_csv(f'{sys.argv[1]}/elastic.net.intercept.csv')


