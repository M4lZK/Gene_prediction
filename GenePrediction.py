import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt

filename = "12864_2006_660_MOESM1_ESM.csv"

all = pd.read_csv(filename)

all["5_MAR_presence"][all["5_MAR_presence"] == "no"] = 0.0
all["5_MAR_presence"][all["5_MAR_presence"] == "yes"] = 1.0
all["3_MAR_presence"][all["3_MAR_presence"] == "no"] = 0.0
all["3_MAR_presence"][all["3_MAR_presence"] == "yes"] = 1.0
all["5_polyA_18_presence"][all["5_polyA_18_presence"] == "no"] = 0.0
all["5_polyA_18_presence"][all["5_polyA_18_presence"] == "yes"] = 1.0
all["5_CCGNN_2_5_presence"][all["5_CCGNN_2_5_presence"] == "no"] = 0.0
all["5_CCGNN_2_5_presence"][all["5_CCGNN_2_5_presence"] == "yes"] = 1.0
all["is_hk"][all["is_hk"] == "no"] = 0.0
all["is_hk"][all["is_hk"] == "yes"] = 1.0
del all["EMBL_transcript_id"]

#all["is_hk"] num_nan = 46459 all = 47229
#all["perc_go_ts_match] num_nan = 21272 all = 47229

data = np.array(all[["cDNA_length","cds_length","exons_nr","5_MAR_presence","3_MAR_presence","5_polyA_18_presence","5_CCGNN_2_5_presence","perc_go_hk_match","perc_go_ts_match","is_hk"]].values)

''' prepare test_set = 10% of known result and train_set = 90% (remaining) '''
num_yes = 0
num_no = 0
for i in data:
    if(i[9] == 1):
        num_yes += 1
    if(i[9] == 0):
        num_no += 1
test_set = []
train_set = []
num_yes_10 = 0
num_no_10 = 0
for i in data:
    if(i[9] == 1 and num_yes_10 <= 0.1*num_yes):
        test_set.append(i)
        num_yes_10 += 1
    elif (i[9] == 0 and num_no_10 <= 0.1*num_no):
        test_set.append(i)
        num_no_10 += 1
    else :
        train_set.append(i)
test_set = np.array(test_set)
train_set = np.array(train_set)

''' prepare sup_train_set by filter only known result from train_set and prepare unsup_train_set by filter nan result from train_set '''
unsup_train_set = []
sup_train_set = []
for i in train_set:
    if(np.isnan(i[9])):
        unsup_train_set.append(i)
    else:
        sup_train_set.append(i)
unsup_train_set = np.array(unsup_train_set)
sup_train_set = np.array(sup_train_set)

''' discretization '''

train_set_cDNA_no_nan = []
for i in train_set:
    if(not np.isnan(i[0])):
        train_set_cDNA_no_nan.append(i[0])
train_set_cDNA_no_nan.sort()
#train_set_cDNA_no_nan = np.array(train_set_cDNA_no_nan)
bin_edge_cDNA = train_set_cDNA_no_nan[0::70]
bin_edge_cDNA = np.unique(bin_edge_cDNA)
#print(bin_edge_cDNA.shape)

train_set_cds_no_nan = []
for i in train_set:
    if(not np.isnan(i[1])):
        train_set_cds_no_nan.append(i[1])
train_set_cds_no_nan.sort()
#train_set_cds_no_nan = np.array(train_set_cds_no_nan)
bin_edge_cds = train_set_cds_no_nan[0::32]
bin_edge_cds = np.unique(bin_edge_cds)

train_set_exons_no_nan = []
for i in train_set:
    if(not np.isnan(i[2])):
        train_set_exons_no_nan.append(i[2])
train_set_exons_no_nan.sort()
#train_set_exons_no_nan = np.array(train_set_exons_no_nan)
bin_edge_exons = train_set_exons_no_nan[0::3000]
bin_edge_exons = np.unique(bin_edge_exons)

train_set_perc_go_hk_match_no_nan = []
for i in train_set:
    if(not np.isnan(i[7])):
        train_set_perc_go_hk_match_no_nan.append(i[7])
train_set_perc_go_hk_match_no_nan.sort()
#train_set_perc_go_hk_match_no_nan = np.array(train_set_perc_go_hk_match_no_nan)
bin_edge_perc_go_hk_match = train_set_perc_go_hk_match_no_nan[0::1500]
bin_edge_perc_go_hk_match = np.unique(bin_edge_perc_go_hk_match)

train_set_perc_go_ts_match_no_nan = []
for i in train_set:
    if(not np.isnan(i[8])):
        train_set_perc_go_ts_match_no_nan.append(i[8])
train_set_perc_go_ts_match_no_nan.sort()
#train_set_perc_go_ts_match_no_nan = np.array(train_set_perc_go_ts_match_no_nan)
bin_edge_perc_go_ts_match = train_set_perc_go_ts_match_no_nan[0::5000]
bin_edge_perc_go_ts_match = np.unique(bin_edge_perc_go_ts_match)

''' MLE '''
cDNA_isHK = []
cDNA_notHK = []
cds_isHK = []
cds_notHK = []
exons_isHK = []
exons_notHK = []
I5_MAR_isHK = []
I5_MAR_notHK = []
I3_MAR_isHK = []
I3_MAR_notHK = []
I5_poly_isHK = []
I5_poly_notHK = []
I5_CCGN_isHK = []
I5_CCGN_notHK = []
perc_hk_isHK = []
perc_hk_notHK = []
perc_ts_isHK = []
perc_ts_notHK = []
num_isHK = 0
num_notHK = 0
for i in sup_train_set:
    if(i[9] == 1):
        cDNA_isHK.append(i[0])
        cds_isHK.append(i[1])
        exons_isHK.append(i[2])
        I5_MAR_isHK.append(i[3])
        I3_MAR_isHK.append(i[4])
        I5_poly_isHK.append(i[5])
        I5_CCGN_isHK.append(i[6])
        perc_hk_isHK.append(i[7])
        perc_ts_isHK.append(i[8])
        num_isHK += 1
    if(i[9] == 0):
        cDNA_notHK.append(i[0])
        cds_notHK.append(i[1])
        exons_notHK.append(i[2])
        I5_MAR_notHK.append(i[3])
        I3_MAR_notHK.append(i[4])
        I5_poly_notHK.append(i[5])
        I5_CCGN_notHK.append(i[6])
        perc_hk_notHK.append(i[7])
        perc_ts_notHK.append(i[8])
        num_notHK += 1

cDNA_isHK = np.array(cDNA_isHK)
cDNA_isHK.sort()
cDNA_notHK = np.array(cDNA_notHK)
cDNA_notHK.sort()
prob_cDNA_isHK = np.bincount(np.digitize(cDNA_isHK,bin_edge_cDNA))/num_isHK
prob_cDNA_notHK= np.bincount(np.digitize(cDNA_notHK,bin_edge_cDNA))/num_notHK

cds_isHK = np.array(cds_isHK)
cds_isHK.sort()
cds_notHK = np.array(cds_notHK)
cds_notHK.sort()
prob_cds_isHK = np.bincount(np.digitize(cds_isHK,bin_edge_cds))/num_isHK
prob_cds_notHK = np.bincount(np.digitize(cds_notHK,bin_edge_cds))/num_notHK

exons_isHK = np.array(exons_isHK)
exons_isHK.sort()
exons_notHK = np.array(exons_notHK)
exons_notHK.sort()
prob_exons_isHK = np.bincount(np.digitize(exons_isHK,bin_edge_exons))/num_isHK
prob_exons_notHK = np.bincount(np.digitize(exons_notHK,bin_edge_exons))/num_notHK

I5_MAR__isHK = np.array(I5_MAR_isHK)
I5_MAR__isHK.sort()
I5_MAR__notHK = np.array(I5_MAR_notHK)
I5_MAR__notHK.sort()
prob_I5_MAR_isHK = np.bincount(np.digitize(I5_MAR_isHK,[0,1]))/num_isHK
prob_I5_MAR_notHK = np.bincount(np.digitize(I5_MAR_notHK,[0,1]))/num_notHK

I3_MAR__isHK = np.array(I3_MAR_isHK)
I3_MAR__isHK.sort()
I3_MAR__notHK = np.array(I3_MAR_notHK)
I3_MAR__notHK.sort()
prob_I3_MAR_isHK = np.bincount(np.digitize(I3_MAR_isHK,[0,1]))/num_isHK
prob_I3_MAR_notHK = np.bincount(np.digitize(I3_MAR_notHK,[0,1]))/num_notHK

I5_poly__isHK = np.array(I5_poly_isHK)
I5_poly__isHK.sort()
I5_poly__notHK = np.array(I5_poly_notHK)
I5_poly__notHK.sort()
prob_I5_poly_isHK = np.bincount(np.digitize(I5_poly_isHK,[0,1]))/num_isHK
prob_I5_poly_notHK = np.bincount(np.digitize(I5_poly_notHK,[0,1]))/num_notHK

I5_CCGN__isHK = np.array(I5_CCGN_isHK)
I5_CCGN__isHK.sort()
I5_CCGN__notHK = np.array(I5_CCGN_notHK)
I5_CCGN__notHK.sort()
prob_I5_CCGN_isHK = np.bincount(np.digitize(I5_CCGN_isHK,[0,1]))/num_isHK
prob_I5_CCGN_notHK = np.bincount(np.digitize(I5_CCGN_notHK,[0,1]))/num_notHK

perc_hk__isHK = np.array(perc_hk_isHK)
perc_hk__isHK.sort()
perc_hk__notHK = np.array(perc_hk_notHK)
perc_hk__notHK.sort()
prob_perc_hk_isHK = np.bincount(np.digitize(perc_hk_isHK,bin_edge_perc_go_hk_match))/num_isHK
prob_perc_hk_notHK = np.bincount(np.digitize(perc_hk_notHK,bin_edge_perc_go_hk_match))/num_notHK

perc_ts__isHK = np.array(perc_ts_isHK)
perc_ts__isHK.sort()
perc_ts__notHK = np.array(perc_ts_notHK)
perc_ts__notHK.sort()
prob_perc_ts_isHK = np.bincount(np.digitize(perc_ts_isHK,bin_edge_perc_go_ts_match))/num_isHK
prob_perc_ts_notHK = np.bincount(np.digitize(perc_ts_notHK,bin_edge_perc_go_ts_match))/num_notHK

''' fix data '''
prob_I3_MAR_isHK = prob_I3_MAR_isHK[1:]
prob_I3_MAR_notHK = prob_I3_MAR_notHK[1:]

prob_I5_MAR_isHK = prob_I5_MAR_isHK[1:]
prob_I5_MAR_notHK = prob_I5_MAR_notHK[1:]

prob_I5_CCGN_isHK = prob_I5_CCGN_isHK[1:]
prob_I5_CCGN_notHK = prob_I5_CCGN_notHK[1:]

prob_I5_poly_isHK = prob_I5_poly_isHK[1:]
prob_I5_poly_notHK = prob_I5_poly_notHK[1:]

prob_cDNA_isHK.resize(bin_edge_cDNA.size)

prob_cds_isHK.resize(bin_edge_cds.size)

prob_exons_isHK = prob_exons_isHK[1:]
prob_exons_notHK = prob_exons_notHK[1:]

prob_perc_hk_isHK = prob_perc_hk_isHK[1:]
prob_perc_hk_notHK = prob_perc_hk_notHK[1:]

prob_perc_ts_isHK = prob_perc_ts_isHK[1:]
prob_perc_ts_notHK = prob_perc_ts_notHK[1:]

#plt.bar(bin_edge_cDNA[:-(bin_edge_cDNA.size-x1_sum.size)],x1_sum, color='b',width=75,alpha=0.2,label="is_hk")
#plt.bar(bin_edge_cDNA,x2_sum, color='r',width=75,alpha=0.2,label="not_hk")
#plt.legend()

#plt.bar(bin_edge_cds[:-(bin_edge_cds.size-x1_sum.size)],x1_sum, color='b',width=75,alpha=0.2,label="is_hk")
#plt.bar(left=bin_edge_cds,height=x2_sum, color='r',width=75,alpha=0.2,label="not_hk")
#plt.legend()

#plt.bar(bin_edge_exons,x1_sum[1:], color='b',width=1,alpha=0.4,label="is_hk")
#plt.bar(bin_edge_exons,x2_sum[1:], color='r',width=1,alpha=0.4,label="not_hk")
#plt.legend()

#plt.bar([0,1],prob_I5_poly_isHK, color='b',width=0.3,alpha=0.4,label="is_hk")
#plt.bar([0,1],prob_I5_poly_notHK, color='r',width=0.3,alpha=0.4,label="not_hk")
#plt.legend()
#plt.show()

#plt.bar([0,1],x1_sum[1:], color='b',width=0.3,alpha=0.4,label="is_hk")
#plt.bar([0,1],x2_sum, color='r',width=0.3,alpha=0.4,label="not_hk")
#plt.legend()

#plt.bar([0,1],x1_sum[1:], color='b',width=0.3,alpha=0.4,label="is_hk")
#plt.bar([0,1],x2_sum, color='r',width=0.3,alpha=0.4,label="not_hk")
#plt.legend()

#plt.bar([0,1],x1_sum[1:], color='b',width=0.3,alpha=0.4,label="is_hk")
#plt.bar([0,1],x2_sum, color='r',width=0.3,alpha=0.4,label="not_hk")
#plt.legend()

#plt.bar(bin_edge_perc_go_hk_match,x1_sum[1:], color='b',width=0.1,alpha=0.4,label="is_hk")
#plt.bar(bin_edge_perc_go_hk_match,x2_sum[1:], color='r',width=0.1,alpha=0.4,label="not_hk")
#plt.legend()

#plt.bar(bin_edge_perc_go_ts_match,x1_sum[1:], color='b',width=0.1,alpha=0.4,label="is_hk")
#plt.bar(bin_edge_perc_go_ts_match,x2_sum[1:], color='r',width=0.1,alpha=0.4,label="not_hk")
#plt.legend()

''' Naive Bayes classification'''
def check_zero(n):
    if(n == 0):
        return 0.00001
    return n

def naive_bayes_class(cDNA,cds,exons,I5_MAR,I3_MAR,I5_poly,I5_CCGN,perc_hk,perc_ts,threshold):
    _cDNA_isHK = prob_cDNA_isHK[np.digitize(cDNA, bin_edge_cDNA)-1]
    _cDNA_notHK = prob_cDNA_notHK[np.digitize(cDNA, bin_edge_cDNA)-1]
    _cds_isHK = prob_cds_isHK[np.digitize(cds, bin_edge_cds)-1]
    _cds_notHK = prob_cds_notHK[np.digitize(cds, bin_edge_cds)-1]
    _exon_isHK = prob_exons_isHK[np.digitize(exons, bin_edge_exons)-1]
    _exon_notHK = prob_exons_notHK[np.digitize(exons, bin_edge_exons)-1]
    _I5_MAR_isHK = prob_I5_MAR_isHK[np.digitize(I5_MAR, [0,1])-1]
    _I5_MAR_notHK = prob_I5_MAR_notHK[np.digitize(I5_MAR, [0, 1])-1]
    _I3_MAR_isHK = prob_I3_MAR_isHK[np.digitize(I3_MAR, [0, 1])-1]
    _I3_MAR_notHK = prob_I3_MAR_notHK[np.digitize(I3_MAR, [0, 1])-1]
    _I5_poly_isHK = prob_I5_poly_isHK[np.digitize(I5_poly, [0, 1])-1]
    _I5_poly_notHK = prob_I5_poly_notHK[np.digitize(I5_poly, [0, 1])-1]
    _I5_CCGN_isHK = prob_I5_CCGN_isHK[np.digitize(I5_CCGN, [0, 1])-1]
    _I5_CCGN_notHK = prob_I5_CCGN_notHK[np.digitize(I5_CCGN, [0, 1])-1]
    _perc_hk_isHK = prob_perc_hk_isHK[np.digitize(perc_hk,bin_edge_perc_go_hk_match)-1]
    _perc_hk_notHK = prob_perc_hk_notHK[np.digitize(perc_hk, bin_edge_perc_go_hk_match)-1]
    _perc_ts_isHK = prob_perc_ts_isHK[np.digitize(perc_ts, bin_edge_perc_go_ts_match)-1]
    _perc_ts_notHK = prob_perc_ts_notHK[np.digitize(perc_ts, bin_edge_perc_go_ts_match)-1]
    p_isHK = num_isHK/(num_isHK+num_notHK)
    p_notHK = num_notHK/(num_isHK+num_notHK)
    _cDNA_isHK = check_zero(_cDNA_isHK)
    _cDNA_notHK = check_zero(_cDNA_notHK)
    _cds_isHK = check_zero(_cds_isHK)
    _cds_notHK = check_zero(_cds_notHK)
    _exon_isHK = check_zero(_exon_isHK)
    _exon_notHK = check_zero(_exon_notHK)
    _I5_MAR_isHK = check_zero(_I5_MAR_isHK)
    _I5_MAR_notHK = check_zero(_I5_MAR_notHK)
    _I3_MAR_isHK = check_zero(_I3_MAR_isHK)
    _I3_MAR_notHK = check_zero(_I3_MAR_notHK)
    _I5_poly_isHK = check_zero(_I5_poly_isHK)
    _I5_poly_notHK = check_zero(_I5_poly_notHK)
    _I5_CCGN_isHK = check_zero(_I5_CCGN_isHK)
    _I5_CCGN_notHK = check_zero(_I5_CCGN_notHK)
    _perc_hk_isHK = check_zero(_perc_hk_isHK)
    _perc_hk_notHK = check_zero(_perc_hk_notHK)
    _perc_ts_isHK = check_zero(_perc_ts_isHK)
    _perc_ts_notHK = check_zero(_perc_ts_notHK)

    result = math.log(p_isHK) + math.log(p_notHK) + math.log(_cDNA_isHK) - math.log(_cDNA_notHK) + math.log(_cds_isHK) - math.log(_cds_notHK) + math.log(_exon_isHK) - math.log(_exon_notHK) + math.log(_I5_MAR_isHK) - math.log(_I5_MAR_notHK) + math.log(_I3_MAR_isHK) - math.log(_I3_MAR_notHK) + math.log(_I5_poly_isHK) - math.log(_I5_poly_notHK) + math.log(_I5_CCGN_isHK) - math.log(_I5_CCGN_notHK) + math.log(_perc_hk_isHK) - math.log(_perc_hk_notHK) + math.log(_perc_ts_isHK) - math.log(_perc_ts_notHK)
    if(result >= threshold):
        return 1
    return 0

def predict(input,threshold):
    result = []
    for x in input:
        for i in x:
            if(np.isnan(i)):
                i = 0
        result.append(naive_bayes_class(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],threshold))
    return result

def evaluate(expected,predicted):
    expected2 = []
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    predicted_positive = 0
    for x in expected:
        expected2.append(x[9])
    for i in range(len(expected2)):
        if(expected2[i] == 1 and predicted[i] == 1):
            tp += 1
            predicted_positive += 1
        elif(expected2[i] == 1 and predicted[i] == 0):
            fn += 1
        elif(expected2[i] == 0 and predicted[i] == 1):
            fp += 1
            predicted_positive += 1
        elif(expected2[i] == 0 and predicted[i] == 0):
            tn += 1
    actual_yes = tp + fn
    actual_no = fp + tn
    if(actual_yes == 0):
        actual_yes = 1
    if(predicted_positive == 0):
        predicted_positive = 1
    recall = tp/actual_yes
    precision = tp/predicted_positive
    accuracy = tn/len(expected2)
    return [accuracy,recall,precision,(2*recall*precision)/(recall+precision),tp/actual_yes,fp/actual_no]

''' main of naive bayes classification'''
t = np.arange(-5,5,0.05)
roc = []
max_precision = -10000
max_precision_thres = -10000
max_f_score = -10000
max_f_score_thres = -10000
max_recall = -10000
max_recall_thres = -10000
max_accuracy = -10000
max_accuracy_thres = -10000

for i in t:
    predicted = predict(test_set,i)
    result = evaluate(test_set,predicted)
    roc.append([result[4],result[5]])
    if(result[0] > max_accuracy):
        max_accuracy = result[0]
        max_accuracy_thres = i
    if(result[1] > max_recall):
        max_recall = result[1]
        max_recall_thres = i
    if(result[2] > max_precision):
        max_precision = result[2]
        max_precision_thres = i
    if(result[3] > max_f_score):
        max_f_score = result[3]
        max_f_score_thres = i

'''
print(max_accuracy)
print(max_accuracy_thres)
print(max_recall)
print(max_recall_thres)
print(max_precision)
print(max_precision_thres)
print(max_f_score)
print(max_f_score_thres)
'''

''' plot graph roc '''
'''
x = []
y = []
for i in roc:
    x.append(i[1])
    y.append(i[0])
x = np.array(x)
y = np.array(y)
plt.plot(x, y)
plt.ylabel("True positive rate")
plt.xlabel("False positive rate")
plt.show()
'''

''' predicted unsup '''
result = predict(unsup_train_set,-4.9)
cDNA = []
cds = []
exon = []
_5_MAR = []
_3_MAR = []
_5_poly = []
_5_CCGN = []
perc_go = []
perc_ts = []


''' predicted unsup to csv '''
to_csv = pd.DataFrame(columns = ["cDNA_length","cds_length","exon_nr","5_MAR_presence","3_MAR_presence","5_polyA_18_presence","5_CCGNN_2_5_presence","perc_go_hk_match","perc_go_ts_match","unsup_predicted"])
for i in unsup_train_set:
    cDNA.append(i[0])
    cds.append(i[1])
    exon.append(i[2])
    _5_MAR.append(i[3])
    _3_MAR.append(i[4])
    _5_poly.append(i[5])
    _5_CCGN.append(i[6])
    perc_go.append(i[7])
    perc_ts.append(i[8])

to_csv["cDNA_length"] = cDNA
to_csv["cds_length"] = cds
to_csv["exon_nr"] = exon
to_csv["5_MAR_presence"] = _5_MAR
to_csv["3_MAR_presence"] = _3_MAR
to_csv["5_polyA_18_presence"] = _5_poly
to_csv["5_CCGNN_2_5_presence"] = _5_CCGN
to_csv["perc_go_hk_match"] = perc_go
to_csv["perc_go_ts_match"] = perc_ts
to_csv["unsup_predicted"] = result
to_csv.to_csv('unsup_predicted.csv',index=False)


'''
random_expected = []
one_sum = 0
zero_sum = 0

for i in test_set:
    if(one_sum <= zero_sum):
        random_expected.append(1)
        one_sum += 1
    elif(zero_sum < one_sum):
        random_expected.append(0)
        zero_sum += 1
result = evaluate(test_set,random_expected)
print(result)
'''
