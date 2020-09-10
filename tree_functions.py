import numpy as np
import pandas as pd
import re
import random
import matplotlib.pyplot as plt
import first_functions as funs
from tqdm import tqdm

def estimates_for_all_trees(file_path):
    trees = open(file_path)
    trees_text = trees.read()
    trees = re.findall("tree [\_\'\ \(\)\=\:\.\,\[\]\{\}\%\&\-0-9a-zA-Z]+\;", trees_text)

    lln_estimates = []
    exp_estimates = []
    for text in trees:
        text = text[:-1]
        ile = sum([d=="(" for d in text ])
        times = []
        for k in range(ile-1):
            text_temp = text.split("(")[-1]
            left = float(text_temp.split(")")[0].split(":")[-1])
            right = "?"

            text2 = ")".join(text_temp.split(")")[1:])
            
            for i,j in enumerate(text2):
                if j=="," or j==")":
                    right = float(re.findall("[0-9Ee\-\.]+", text2)[0])
                    text2 = text2[i:]
                    break

            text = "(".join(text.split("(")[:-1])+"a:"+str(left+right)+str(text2)

            times.append(left)
            if k==ile-2:
                times.append(left+right)

        #times = list(np.unique([round(t,4) for t in times]))
        lln_estimates.append(funs.estimate_theta_lln(times))
        exp_estimates.append(funs.estimate_theta_exp(times))

    return(lln_estimates, exp_estimates)
    
    
def estimates_for_single_tree(file_path):
    trees = open(file_path)
    trees_text = trees.read()
    trees = re.findall("tree [\_\'\ \(\)\=\:\.\,\[\]\{\}\%\&\-0-9a-zA-Z]+\;", trees_text)

    text = trees[0][:-1]

    ile = sum([d=="(" for d in text ])
    times = []
    for k in range(ile-1):
        text_temp = text.split("(")[-1]
        left = float(text_temp.split(")")[0].split(":")[-1])
        right = "?"

        text2 = ")".join(text_temp.split(")")[1:])
        
        for i,j in enumerate(text2):
            if j=="," or j==")":
                right = float(re.findall("[0-9Ee\-\.]+", text2)[0])
                text2 = text2[i:]
                break

        text = "(".join(text.split("(")[:-1])+"a:"+str(left+right)+str(text2)

        times.append(left)
        if k==ile-2:
            times.append(left+right)

    M = len(times)+1
    
    return( funs.estimate_theta_lln(times), funs.estimate_theta_exp(times), M, times ) 
 
    
# ----------------------------------------------- same functions with individuals sampling --------------------------------

def get_number(code):
    return int(re.findall("[0-9]+", code)[0])

def estimates_for_all_trees_with_sampling(file_path, sample_size=None):
    """
    sample_choice is boolean vector. If none all of the individuals are sampled.
    """
    
    trees = open(file_path)
    trees_text = trees.read()
    trees = re.findall("tree [\_\'\ \(\)\=\:\.\,\[\]\{\}\%\&\-0-9a-zA-Z]+\;", trees_text)
    ile = sum([d=="(" for d in trees[0] ])
    if sample_size==None:
        sample_size = ile+1
    if sample_size==-1:
        sample_size = ile
    elif sample_size>ile+1:
        raise(Exception("Sample size is bigger than the number of individuals, which is "+str(ile+1)+"!"))
    
    sample_choice=[True]*sample_size+[False]*(ile+1-sample_size)
    random.shuffle(sample_choice)
    
    lln_estimates = []
    exp_estimates = []
    for text in trees:
        text = text[:-1]
        times = []
        for k in range(ile-1):
            text_temp = text.split("(")[-1]
            left = float(text_temp.split(")")[0].split(":")[-1])
            name1 = text_temp.split(")")[0].split(":")[0]
            name2 = text_temp.split(")")[0].split(":")[-2].split(",")[-1]
            right = "?"

            text2 = ")".join(text_temp.split(")")[1:])
            
            for i,j in enumerate(text2):
                if j=="," or j==")":
                    right = float(re.findall("[0-9Ee\-\.]+", text2)[0])
                    text2 = text2[i:]
                    break
            
            if sample_choice[get_number(name1)-1]:
                new_name = name1
            else:
                new_name = name2
                
            text = "(".join(text.split("(")[:-1])+new_name+":"+str(left+right)+str(text2)
                
            if sample_choice[get_number(name1)-1] and sample_choice[get_number(name2)-1]:
                times.append(left)
            
            if k==ile-2 and len(times) < sample_size-1:
                times.append(left+right)

        lln_estimates.append(funs.estimate_theta_lln(times))
        exp_estimates.append(funs.estimate_theta_exp(times))
    return(lln_estimates, exp_estimates)

def estimates_for_single_tree_with_sampling(file_path, sample_size=None, tree_text=None):
    if tree_text==None:
        trees = open(file_path)
        trees_text = trees.read()
        trees = re.findall("tree [\_\'\ \(\)\=\:\.\,\[\]\{\}\%\&\-0-9a-zA-Z]+\;", trees_text)
        text = trees[0][:-1]
    else:
        text = tree_text
    ile = sum([d=="(" for d in text ])
    if sample_size==None:
        sample_size = ile+1
    if sample_size==-1:
        sample_size = ile
    elif sample_size>ile+1:
        raise(Exception("Sample size is bigger than the number of individuals, which is "+str(ile+1)+"!"))
    
    sample_choice=[True]*sample_size+[False]*(ile+1-sample_size)
    random.shuffle(sample_choice)
    
    times = []
    for k in range(ile-1):
        text_temp = text.split("(")[-1]
        left = float(text_temp.split(")")[0].split(":")[-1])
        name1 = text_temp.split(")")[0].split(":")[0]
        name2 = text_temp.split(")")[0].split(":")[-2].split(",")[-1]
        right = "?"

        text2 = ")".join(text_temp.split(")")[1:])

        for i,j in enumerate(text2):
            if j=="," or j==")":
                right = float(re.findall("[0-9Ee\-\.]+", text2)[0])
                text2 = text2[i:]
                break

        if sample_choice[get_number(name1)-1]:
            new_name = name1
        else:
            new_name = name2

        text = "(".join(text.split("(")[:-1])+new_name+":"+str(left+right)+str(text2)

        if sample_choice[get_number(name1)-1] and sample_choice[get_number(name2)-1]:
            times.append(left)

        if k==ile-2 and len(times) < sample_size-1:
            times.append(left+right)

    M = len(times)+1
    
    return( funs.estimate_theta_lln(times), funs.estimate_theta_exp(times), M, times ) 

# ---------------------------- generating input to hybrid-lambda and seq-gen ------------------------------

def generate_topology(times):
    M = len(times) + 1
    branches = [ ("S"+str(i),0) for i in range(M) ]
    for i in range(M-1):
        m = M-i
        to_join = random.sample(range(m), 2)
        b1 = branches[to_join[0]]
        b2 = branches[to_join[1]]
        joined_branch = ("("+b1[0]+":"+str(times[i]-b1[1])+","+b2[0]+":"+str(times[i]-b2[1])+")",times[i])
        branches = [ j for i,j in enumerate(branches) if i not in to_join ] + [joined_branch]
    return branches[0]

def gen_hybrid_input(numer_testu, N, theta=1, disturb=False, disturb_variance=0.002):
    times = funs.generate_kingman_times(N, theta)
    if disturb:
        times = funs.disturb_times_beta(times, disturb_variance)
    print(times)
#     tree_text = "("*(N-1) + "S0"
#     pre_time = times[0]
#     for i,t in enumerate(times):
#         if i==0:
#             tree_text += ":"+str(pre_time)+",S"+str(i+1)+":"+str(times[i])+")"
#         else:
#             diff = times[i]-pre_time
#             tree_text += ":"+str(diff)+",S"+str(i+1)+":"+str(pre_time+diff)+")"
#             pre_time += diff
#     tree_text += "r;"
    tree_text = generate_topology(times)[0] + "r;"
    
    print("cd ~/Documents/Studia/PracaMag/BEAST_files/ULTIMATE/trees_hybrid/")
    text = "~/Documents/Github/hybrid-Lambda/hybrid-Lambda  -spcu '"+tree_text+"' -o test"+str(numer_testu)+" -pop 2500 -mu 0.0001 -num 20 -sim_mut_unit"
    print(text)
    print("cd ~/Documents/Studia/PracaMag/BEAST_files/ULTIMATE/seqs_seqgen/")
    print("seq-gen < ../trees_hybrid/test"+str(numer_testu)+"_mut_unit > test"+str(numer_testu)+".nex -m HKY -l 500 -n 40 -on")

def gen_times_from_tree(text):
    ile = sum([d=="(" for d in text ]) 
    times = []
    for k in range(ile-1):
        text_temp = text.split("(")[-1]
        left = float(text_temp.split(")")[0].split(":")[-1])
        name1 = text_temp.split(")")[0].split(":")[0]
        name2 = text_temp.split(")")[0].split(":")[-2].split(",")[-1]

        text2 = ")".join(text_temp.split(")")[1:])

        for i,j in enumerate(text2):
            if j=="," or j==")":
                right = float(re.findall("[0-9Ee\-\.]+", text2)[0])
                text2 = text2[i:]
                break

        text = "(".join(text.split("(")[:-1])+"a:"+str(left+right)+str(text2)

        times.append(left)

        if k==ile-2:
            times.append(left+right)
    return times
    
# ---------------------------------- printing info about tests --------------------------------
    
def print_test_info(numer_testu):
    if numer_testu>9:
        nr_str = str(numer_testu)
    else:
        nr_str = "0"+str(numer_testu)

    lln, exp, M, times = estimates_for_single_tree("../BEAST_files/ULTIMATE/exported/test"+nr_str+"_single")
    print("Liczba linii: "+str(M))
    print("EXP: ", exp)
    print("LLN: ", lln)

    lln_estimates, exp_estimates = estimates_for_all_trees("../BEAST_files/ULTIMATE/exported/test"+nr_str+"_all")
    f = 50
    plt.scatter(lln_estimates[f:], exp_estimates[f:])
    plt.ylabel("EXP estimations")
    plt.xlabel("LLN estimations")
    plt.show()

    print("EXP:   mean: ", round(np.mean(exp_estimates[f:]),4), ", Variance: ", round(np.std(exp_estimates[f:])**2,4) )
    print("LLN:   mean: ", round(np.mean(lln_estimates[f:]),4), ", Variance: ", round(np.std(lln_estimates[f:])**2,4) )
    print("EXP MSE: ", round((np.mean(exp_estimates[f:])-1)**2+np.std(exp_estimates[f:])**2,5))
    print("LLN MSE: ", round((np.mean(lln_estimates[f:])-1)**2+np.std(lln_estimates[f:])**2,5))
    
# -------------------- printing info about sampling -----------------

def get_tree_size(file_path):
    trees = open(file_path)
    trees_text = trees.read()
    trees = re.findall("tree [\_\'\ \(\)\=\:\.\,\[\]\{\}\%\&\-0-9a-zA-Z]+\;", trees_text)
    text = trees[0][:-1]
    return sum([d=="(" for d in text ])+1

def print_sampling_info(numer_testu):
    file_path2 = "../BEAST_files/ULTIMATE/exported/test"+str(numer_testu)+"_single"
    
    # data for plot for a single run
    x = np.arange(2,get_tree_size(file_path2)+1)
    lln = []
    exp = []
    for x1 in x:
        lln1, exp1, M, t = estimates_for_single_tree_with_sampling(file_path2,x1)
        lln.append(lln1)
        exp.append(exp1)
    
    # data for the averaged runs
    results = pd.DataFrame(columns=["exp"+str(i) for i in x]+["lln"+str(i) for i in x])
    for i in tqdm(range(100)):
        lln_temp = []
        exp_temp = []
        for x1 in x:
            lln1, exp1, M, t = estimates_for_single_tree_with_sampling(file_path2, x1)
            lln_temp.append(lln1)
            exp_temp.append(exp1)
        results.loc[i,:] = lln_temp+exp_temp
    print(M)
    print(t)
    print(funs.estimate_theta_exp(t))
    print(funs.estimate_theta_lln(t))
    
    means = results.mean(axis=0)
    stds = results.std(axis=0)
    lln2 = means[:(len(means)//2)]
    exp2 = means[(len(means)//2):]
    lln_ci = stds[:(len(stds)//2)]
    exp_ci = stds[(len(stds)//2):]

    # plotting
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(15,5))
    ax1.plot(x,lln,c="r",linewidth=2)
    ax1.plot(x,exp,c="b",linewidth=2)
    ax1.hlines(1,2,len(lln)+1,linestyle="--",linewidth=2)
    ax1.set_xlabel("Sample size")
    ax1.set_ylabel("$\\hat{\\theta}$")
    ax1.set_title("Single run")
    ax1.legend(loc="lower right", labels=["LLN method","EXP method"])

    ax2.plot(x,lln2,c="r")
    ax2.plot(x,exp2,c="b")
    ax2.hlines(1,2,len(lln)+1,linestyle="--")
    ax2.fill_between(x, (lln2-lln_ci), (lln2+lln_ci), color='r', alpha=.1)
    ax2.fill_between(x, (exp2-exp_ci), (exp2+exp_ci), color='b', alpha=.1)
    ax2.set_xlabel("Sample size")
    ax2.set_ylabel("$\\hat{\\theta}$")
    ax2.set_title("Average over 100 runs")
    ax2.legend(loc="lower right", labels=["LLN method","EXP method"])

    plt.tight_layout()
    plt.style.use('fivethirtyeight')
    plt.show()
    
def print_multiple_sampling_info(M, tests_count=6, theta=1, disturb=False, disturb_variance=0.002, single_times=False):
    file_path2 = "../BEAST_files/ULTIMATE/exported/test60_single"
    
    rows = int(np.ceil(tests_count/2))
    fig, ax = plt.subplots(rows, 2, figsize=(15,5*np.ceil(tests_count/2)))
    
    times = funs.generate_kingman_times(M, theta)
    if disturb:
        times = funs.disturb_times_beta(times, disturb_variance)
    
    for test in tqdm(range(tests_count)):
        tree_text = generate_topology(times)[0]
        if not single_times:
            times = funs.generate_kingman_times(M, theta)
            if disturb:
                times = funs.disturb_times_beta(times, disturb_variance)
        
        # data for the averaged runs
        x = np.arange(2,M+1)
        results = pd.DataFrame(columns=["exp"+str(i) for i in x]+["lln"+str(i) for i in x])
        for i in range(100):
            lln_temp = []
            exp_temp = []
            for x1 in x:
                lln1, exp1, M, t = estimates_for_single_tree_with_sampling(file_path2, x1, tree_text)
                lln_temp.append(lln1)
                exp_temp.append(exp1)
            results.loc[i,:] = lln_temp+exp_temp

        means = results.mean(axis=0)
        stds = results.std(axis=0)
        lln2 = means[:(len(means)//2)]
        exp2 = means[(len(means)//2):]
        lln_ci = stds[:(len(stds)//2)]
        exp_ci = stds[(len(stds)//2):]

        # plotting
        ax[test//2][test%2].plot(x,lln2,c="r")
        ax[test//2][test%2].plot(x,exp2,c="b")
        ax[test//2][test%2].hlines(theta, 2, len(lln2)+1, linestyle="--")
        ax[test//2][test%2].fill_between(x, (lln2-lln_ci), (lln2+lln_ci), color='r', alpha=.1)
        ax[test//2][test%2].fill_between(x, (exp2-exp_ci), (exp2+exp_ci), color='b', alpha=.1)
        ax[test//2][test%2].set_xlabel("Sample size")
        ax[test//2][test%2].set_ylabel("$\\hat{\\theta}$")
        ax[test//2][test%2].legend(loc="lower right", labels=["LLN method","EXP method"])

    #plt.tight_layout()
    #plt.style.use('fivethirtyeight')
    plt.show()