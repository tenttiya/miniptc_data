import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pylab as plt

plt.rc('font', family='Arial')

def generate_ptcs(seq, start1, end1, start2, end2, start3, end3):
    '''
    Returns ptcs from the given sequence and intervals
    Ex:
        seq = 'cgaucgc'
        start1 = 43
        end1 = 45
        start2 = 250
        end2 = 251
        start3 = 400
        end3 = 401
        ptcs = ['c43', 'g44', 'a45', 'u250', 'c251', 'g400', 'c401']
    '''
    l1 = list(range(start1, end1+1))
    l2 = list(range(start2, end2+1))
    l3 = list(range(start3, end3+1))
    l = l1 + l2 + l3
    ptcs = []
    for idx in range(len(l)):
        ptc = seq[idx].upper() + str(l[idx])
        ptcs.append(ptc)
    return ptcs

def read_reactivity(fn):
    with open(fn) as f:
        text = f.read()
    lines = text.split('\n')
    result = []
    for line in lines:
        if not line == '':
            result.append(float(line))
    return result

def get_subset(reacts, start, end):
    '''
    get items from [start, end] inclusive
    '''
    react = reacts[start - 1:end]
    return react

def get_combined_subsets(react_files, start, end):
    reacts_subsets = []
    for react_file in react_files:
        reacts = read_reactivity(react_file)
        reacts_subset = get_subset(reacts, start, end)
        reacts_subsets.append(reacts_subset)
    combine = np.array(reacts_subsets)
    return combine

def normalize(value_old, start_old, end_old, start_new, end_new):
    '''
    https://imgur.com/a/Dy3anFJ
    '''
    value_new = (value_old - start_old)/(end_old - start_old)*(end_new - start_new) + start_new
    return value_new

def get_cmap(anchor_values, anchor_colors, vmin, vmax):
    clist = []
    for anchor_value, anchor_color in zip(anchor_values, anchor_colors):
        normalized_anchor_value = normalize(anchor_value, vmin, vmax, 0, 1)
        clist.append((normalized_anchor_value, anchor_color))
    cm = matplotlib.colors.LinearSegmentedColormap.from_list('test', clist, N=256)
    return cm

if __name__ == '__main__':
    # generate combined subsets
    react_files = [
        #'data/predicted_norm.txt'    comment predicted.txt because it has only 284 nt due to the lack of 5' and 3' flanking regions that present in mini-PTC constructs.
        'data/reacts.txt',
        'data/reacts_1-1.txt',
        'data/reacts_1-2.txt',
        'data/reacts_1-3.txt',
        'data/reacts_1-4.txt',
        'data/reacts_1-5.txt',
        'data/reacts_1-6.txt',
        'data/reacts_1-7.txt'
        ]
    combine = get_combined_subsets(react_files, 24, 307)

    # generate ptcs for xticks
    seq = 'cccgcggcaagacggaaagaccccguaaaccuuuacuauagcuugacacaaaagugucugguggguaguuugacuggggcggucacggauaaaagguacuccggggauaacaggcugauaccgcccaagaguucauaucgacggcgguguuuggcaccucgaugucggcucaucacauccuggggcugaaguaggucccaaggguauggcuguucgccauuuaaagugguacgcgagcuggguuuagaacgucgugagacaguucggucccuaucugccguggg'
    start1 = 2043
    end1 = 2095
    start2 = 2228
    end2 = 2258
    start3 = 2426
    end3 = 2625
    ptcs = generate_ptcs(seq, start1, end1, start2, end2, start3, end3)

    # generate color map
    vmin = 0.5
    vmax = 1.5
    anchor_values = [0.5, 1.0, 1.5]
    anchor_colors = ['#ffffff', '#ffc107', '#f44336']
    cm = get_cmap(anchor_values, anchor_colors, vmin, vmax)

    # plot heatmap
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(12, 12))
    yticklabels = ['WT','1.1','1.2','1.3','1.4','1.5','1.6','1.7']

    combine_small = combine[:, 0:71]
    ptcs_small = ptcs[0:71]
    ax = sns.heatmap(combine_small, vmin=vmin, vmax=vmax, linewidth=0.5, linecolor='white', cmap=cm, annot=False, square=True, ax=axs[0], xticklabels=ptcs_small, yticklabels=yticklabels)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=9)

    combine_small = combine[:, 71:71*2]
    ptcs_small = ptcs[71:71*2]
    ax = sns.heatmap(combine_small, vmin=vmin, vmax=vmax, linewidth=0.5, linecolor='white', cmap=cm, annot=False, square=True, ax=axs[1], xticklabels=ptcs_small, yticklabels=yticklabels)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=9)

    combine_small = combine[:, 71*2:71*3]
    ptcs_small = ptcs[71*2:71*3]
    ax = sns.heatmap(combine_small, vmin=vmin, vmax=vmax, linewidth=0.5, linecolor='white', cmap=cm, annot=False, square=True, ax=axs[2], xticklabels=ptcs_small, yticklabels=yticklabels)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=9)

    combine_small = combine[:, 71*3:71*4]
    ptcs_small = ptcs[71*3:71*4]
    ax = sns.heatmap(combine_small, vmin=vmin, vmax=vmax, linewidth=0.5, linecolor='white', cmap=cm, annot=False, square=True, ax=axs[3], xticklabels=ptcs_small, yticklabels=yticklabels)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=9)

    fig.tight_layout()
    plt.savefig("040523_heatmap_overall.png",dpi=600)
    plt.show()