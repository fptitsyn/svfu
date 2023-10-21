import numpy as np
import data
from sklearn.naive_bayes import GaussianNB

labels = np.array(data.getLabels())
side_effects = np.array(data.getSideEffects())
test_labels = np.array(data.getTestLabels())
test_side_effects = np.array(data.getTestSideEffects())
labels = np.array(np.split(labels, 15), dtype="float64")
test_labels = np.array(np.split(test_labels, 10), dtype="float64")
n = 22

for i in range(n):
    gnb = GaussianNB()
    current_effect = list()
    for j in range(len(side_effects)):
        current_effect.append(side_effects[j][i])
    gnb.fit(labels, current_effect)
    predict = gnb.predict(test_labels)

    print(f"{i + 1}. {data.get_effect_name(i)}: {predict}")
