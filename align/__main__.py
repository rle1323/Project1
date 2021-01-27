import numpy as np
import matplotlib.pyplot as plt
from align import algs
import statistics as stats
import sklearn.metrics as metrics
import seaborn as sns


# read in pospairs and keep them as a list of lists 
pos_pairs = []
with open("scoring_matrices/pospairs.txt") as f:
    for line in f.readlines():
        line_items = line.strip().split()
        pos_pairs.append(line_items)
# do the same with negpairs
neg_pairs = []
with open("scoring_matrices/negpairs.txt") as f:
    for line in f.readlines():
        line_items = line.strip().split()
        neg_pairs.append(line_items)

'''
################Question 1#######################
#make a local aligner
local_alignment_tool = algs.SmithWaterman(sub_matrix_file = "scoring_matrices/BLOSUM50.mat", gap_start_penalty = 11, gap_extension_penalty = 3)

# iterate through the pairs and align them to generate alignment scores
alignment_scores = []
for pair in pos_pairs:
    _,score = local_alignment_tool.align(pair[0], pair[1])
    alignment_scores.append(score)
for pair in neg_pairs:
    _,score = local_alignment_tool.align(pair[0], pair[1])
    alignment_scores.append(score)

# produce histogram of the alignment scores
plt.hist(alignment_scores, bins = 25)
plt.xlabel("Optimal local alignment score")
plt.ylabel("Count")
plt.show()


################Question 2#######################
#make a local aligner
local_alignment_tool = algs.SmithWaterman(sub_matrix_file = "scoring_matrices/BLOSUM50.mat", gap_start_penalty = 11, gap_extension_penalty = 3)

# iterate through the pairs and align them to generate alignment scores
alignment_scores = []
ground_truth_scores = []
for pair in pos_pairs:
    _,score = local_alignment_tool.align(pair[0], pair[1])
    alignment_scores.append(score)
    ground_truth_scores.append(1)
for pair in neg_pairs:
    _,score = local_alignment_tool.align(pair[0], pair[1])
    alignment_scores.append(score)
    ground_truth_scores.append(0)

# Average alignment score is threshold of "positive" or "negative" scoring
threshold = stats.mean(alignment_scores)

# score each alignment as a positive or negative
predicted_scores = []
for score in alignment_scores:
    if score >= threshold:
        predicted_scores.append(1)
    else:
        predicted_scores.append(0)

# now indicate whether a predictions is fp, tp, fn, or tn
fp = 0
tp = 0
fn = 0
tn = 0
for pred, true in zip(predicted_scores, ground_truth_scores):
    if pred == 1 and true == 1:
        tp += 1
    elif pred == 1 and true == 0:
        fp += 1
    elif pred == 0 and true == 1:
        fn += 1
    elif pred == 0 and true == 0:
        tn += 1
print(threshold)
print(tp, fp)
print(fn, tn)

################Question 3#######################
#make a local aligner
local_alignment_tool = algs.SmithWaterman(sub_matrix_file = "scoring_matrices/BLOSUM50.mat", gap_start_penalty = 11, gap_extension_penalty = 3)

# iterate through the pairs and align them to generate alignment scores
alignment_scores = []
ground_truth_scores = []
for pair in pos_pairs:
    _,score = local_alignment_tool.align(pair[0], pair[1])
    alignment_scores.append(score)
    ground_truth_scores.append(1)
for pair in neg_pairs:
    _,score = local_alignment_tool.align(pair[0], pair[1])
    alignment_scores.append(score)
    ground_truth_scores.append(0)

# Average alignment score is threshold of "positive" or "negative" scoring
threshold = stats.mean(alignment_scores)

# score each alignment as a positive or negative
predicted_scores = []
for score in alignment_scores:
    if score >= threshold:
        predicted_scores.append(1)
    else:
        predicted_scores.append(0)

# compute needed fpr and tpr for roc plot
fpr, tpr, _ = metrics.roc_curve(ground_truth_scores,  predicted_scores)
auc = metrics.roc_auc_score(ground_truth_scores, predicted_scores)

# make ROC plot
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(fpr,tpr,label="Local alignment classifier, AUROC="+str(auc))
plt.legend(loc=4)
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.xlim(0,1)
plt.ylim(0,1)
ax.set_aspect('equal', adjustable='box')
plt.show()



################Question 5#######################

auc_array = np.zeros([5, 20])
for ge in range(1,6):
    print("HELLO")
    for go in range(1,21):
        local_alignment_tool = algs.SmithWaterman(sub_matrix_file = "scoring_matrices/BLOSUM62.mat", gap_start_penalty = go, gap_extension_penalty = ge)

        # iterate through the pairs and align them to generate alignment scores
        alignment_scores = []
        ground_truth_scores = []
        for pair in pos_pairs:
            _,score = local_alignment_tool.align(pair[0], pair[1])
            alignment_scores.append(score)
            ground_truth_scores.append(1)
        for pair in neg_pairs:
            _,score = local_alignment_tool.align(pair[0], pair[1])
            alignment_scores.append(score)
            ground_truth_scores.append(0)

        # Average alignment score is threshold of "positive" or "negative" scoring
        threshold = stats.mean(alignment_scores)

        # score each alignment as a positive or negative
        predicted_scores = []
        for score in alignment_scores:
            if score >= threshold:
                predicted_scores.append(1)
            else:
                predicted_scores.append(0)

        # compute needed fpr and tpr for roc plot
        auc = metrics.roc_auc_score(ground_truth_scores, predicted_scores)
        auc_array[ge - 1,go - 1] = auc

sns.heatmap(auc_array, vmin=0, vmax=1, annot = True)
plt.xlabel("(Gap open penalty - 1)")
plt.ylabel("(Gap extension penalty - 1)")
plt.show()




################Question 7#######################
# ROC plots for LOCAL ALIGNMENT
matrices = ["BLOSUM50", "BLOSUM62", "PAM100", "PAM250"]
for matrix in matrices:
    mat_path = "scoring_matrices/" + matrix + ".mat"
    local_alignment_tool = algs.SmithWaterman(sub_matrix_file = mat_path, gap_start_penalty = 4, gap_extension_penalty = 5)

    # iterate through the pairs and align them to generate alignment scores
    alignment_scores = []
    ground_truth_scores = []
    for pair in pos_pairs:
        _,score = local_alignment_tool.align(pair[0], pair[1])
        alignment_scores.append(score)
        ground_truth_scores.append(1)
    for pair in neg_pairs:
        _,score = local_alignment_tool.align(pair[0], pair[1])
        alignment_scores.append(score)
        ground_truth_scores.append(0)

    # Average alignment score is threshold of "positive" or "negative" scoring
    threshold = stats.mean(alignment_scores)

    # score each alignment as a positive or negative
    predicted_scores = []
    for score in alignment_scores:
        if score >= threshold:
            predicted_scores.append(1)
        else:
            predicted_scores.append(0)

    # compute needed fpr and tpr for roc plot
    fpr, tpr, _ = metrics.roc_curve(ground_truth_scores,  predicted_scores)
    auc = metrics.roc_auc_score(ground_truth_scores, predicted_scores)

    # make ROC plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(fpr,tpr,label="Local alignment classifier, AUROC="+str(auc))
    plt.title("ROC Curve, local alignment with " + matrix)
    plt.legend(loc=4)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.xlim(0,1)
    plt.ylim(0,1)
    ax.set_aspect('equal', adjustable='box')
    plt.show()
'''

# ROC plots for GLOBAL ALIGNMENT
matrices = ["BLOSUM50", "BLOSUM62", "PAM100", "PAM250"]
for matrix in matrices:
    mat_path = "scoring_matrices/" + matrix + ".mat"
    global_alignment_tool = algs.NeedlemanWunsch(sub_matrix_file = mat_path, gap_start_penalty = 4, gap_extension_penalty = 5)

    # iterate through the pairs and align them to generate alignment scores
    alignment_scores = []
    ground_truth_scores = []
    for pair in pos_pairs:
        _,score = global_alignment_tool.align(pair[0], pair[1])
        alignment_scores.append(score)
        ground_truth_scores.append(1)
    for pair in neg_pairs:
        _,score = global_alignment_tool.align(pair[0], pair[1])
        alignment_scores.append(score)
        ground_truth_scores.append(0)

    # Average alignment score is threshold of "positive" or "negative" scoring
    threshold = stats.mean(alignment_scores)

    # score each alignment as a positive or negative
    predicted_scores = []
    for score in alignment_scores:
        if score >= threshold:
            predicted_scores.append(1)
        else:
            predicted_scores.append(0)

    # compute needed fpr and tpr for roc plot
    fpr, tpr, _ = metrics.roc_curve(ground_truth_scores,  predicted_scores)
    auc = metrics.roc_auc_score(ground_truth_scores, predicted_scores)

    # make ROC plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(fpr,tpr,label="Global alignment classifier, AUROC="+str(round(auc,2)))
    plt.title("ROC Curve, global alignment with " + matrix)
    plt.legend(loc=4)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.xlim(0,1)
    plt.ylim(0,1)
    ax.set_aspect('equal', adjustable='box')
    plt.show()