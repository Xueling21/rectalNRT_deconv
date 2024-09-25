import array

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import itertools  # Python内置模块，用于高效循环创建迭代器
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from sklearn.metrics import confusion_matrix
import warnings
warnings.filterwarnings('ignore')

#--------------------------方法变量修改-------------------------------
path_input='./'
path_output='./result/'
name_inputfile=''
name_outputfile='bulk_result'

if not os.path.isdir(path_output):
    os.makedirs(path_output)

#--------------------------自定义函数---------------------------------
# 绘制混淆矩阵
def plot_confusion_matrix(cm, classes,title='Confusion matrix', cmap=plt.cm.Blues):
    plt.figure(figsize=(6, 6))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title, fontsize=20)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, fontsize=20, rotation=0)
    plt.yticks(tick_marks, classes, fontsize=20)
    # 不显示网格
    plt.grid(False)
    threshold = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        # 若对应格子上面的数量超过阈值，则上面的字体为白色，为了方便查看
        plt.text(j, i, cm[i, j],
                 horizontalalignment="center",
                 fontsize=25,
                 color="white" if cm[i, j] > threshold else "black")
    plt.ylabel('True', fontsize=25)
    plt.xlabel('Predicted', fontsize=25)
    plt.tight_layout()
    plt.show()

# 计算灵敏度、特异度
def CalculateSensitivitySpecificity(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    # 切片操作，获取每一个类别各自的 tn, fp, tp, fn
    tn_sum = cm[0, 0]  # True Negative
    fp_sum = cm[0, 1]  # False Positive
    tp_sum = cm[1, 1]  # True Positive
    fn_sum = cm[1, 0]  # False Negative
    Condition_positive = tp_sum + fn_sum + 1e-6
    Condition_negative = tn_sum + fp_sum + 1e-6
    sensitivity = tp_sum / Condition_positive
    specificity = tn_sum / Condition_negative
    return sensitivity, specificity

# 打印模型性能报告
def check_perform(model, X_train, X_test, y_train, y_test):
    y_train_predict = model.predict(X_train)
    print(f'\nOn Trainning Datasets:\n')
    print(f'accuracy_score is:{accuracy_score(y_train, y_train_predict)}')  # 准确率
    print(f'precison_score is:{precision_score(y_train, y_train_predict)}')  # 精确率
    print(f'recall_score is:{recall_score(y_train, y_train_predict)}')  # 召回率
    print(f'sensitivity_score is:{CalculateSensitivitySpecificity(y_train, y_train_predict)[0]}')  # 灵敏度
    print(f'specificity_score is:{CalculateSensitivitySpecificity(y_train, y_train_predict)[1]}')  # 特异度
    y_score = model.predict_proba(X_train)[:, 1]
    print(f'auc is:{roc_auc_score(y_train, y_score)}')
    # AUC(Area under Curve)：Roc曲线下的面积，介于0.1和1之间。Auc作为数值可以直观的评价分类器的好坏，值越大越好。
    print(f'\n===========================\n')
    y_test_predict = model.predict(X_test)
    print(f'\nOn Testing Datasets:\n')
    print(f'accuracy_score is:{accuracy_score(y_test, y_test_predict)}')
    print(f'precison_score is:{precision_score(y_test, y_test_predict)}')
    print(f'recall_score is:{recall_score(y_test, y_test_predict)}')
    print(f'sensitivity_score is:{CalculateSensitivitySpecificity(y_test, y_test_predict)[0]}')
    print(f'specificity_score is:{CalculateSensitivitySpecificity(y_test, y_test_predict)[1]}')
    y_score=model.predict_proba(X_test)[:, 1]
    print(f'auc is:{roc_auc_score(y_test, y_score)}')


# --------------------------数据准备阶段----------------------------
np.random.seed(42)
data_raw = pd.read_csv("cs_predict.txt", delimiter='\t')

data=data_raw.copy()

df = ["NR" for _ in range(25)] + ["R" for _ in range(42 - 25)]
df = pd.DataFrame(df, columns=['Target'])
print(df)
data=pd.concat([data, df], axis=1)

df=list()
for i in range(0,len(data)):
    if(data["NR"][i]>data["R"][i]):
        df.append("NR")
    else:
        df.append("R")

df = pd.DataFrame(df, columns=['Predict'])
data=pd.concat([data, df], axis=1)


cm = confusion_matrix(data['Target'], data['Predict'])
plot_confusion_matrix(cm,["NR", "R"], title="csDEG Confusion Matrix")

print(f'accuracy_score is:{accuracy_score(data["Target"], data["Predict"])}')
print(f'sensitivity_score is:{CalculateSensitivitySpecificity(data["Target"], data["Predict"])[0]}')
print(f'specificity_score is:{CalculateSensitivitySpecificity(data["Target"], data["Predict"])[1]}')






