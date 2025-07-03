#!/usr/bin/env python

import numpy as np
from tqdm import tqdm
import json
import matplotlib.pyplot as plt
import argparse
import os

np.random.seed(42)

def get_args():
    parser = argparse.ArgumentParser(description='EM proportions estimations')
    parser.add_argument('-o', '--output_dir', help='Output directory',  default='output')
    parser.add_argument('-n', '--num_iter', help='Output directory',  default=1000, type=int)
    args = parser.parse_args()

    return args
    

# Base EM algorithm
def base_em(L, num_iter=100, init_props=None):
    num_srcs = L.shape[1]
    # Initialize Pr(S|R) := Pr(R|S)Pr(S)/Pr(R)
    w = np.zeros_like(L)
    # Initialize Pr(S)for each EM iteration
    props = np.zeros((num_iter+1, num_srcs))

    # Log likelihoods
    lls = np.zeros(num_iter+1)
    # Set starting proportions to uniform or some given proportions (like Naive MLEs)
    if isinstance(init_props, np.ndarray):
        props[0] = init_props
    else:
        props[0] = np.ones(num_srcs)*(1/num_srcs)

    # Start EM
    for i in tqdm(range(num_iter)):
        # E-step: updating Pr(S|R)
        w = props[i]*L
        clik = np.sum(w, axis=1)[:, np.newaxis]
        w /= clik
        # M-step: New estimates of Pr(S) 
        props[i+1] = np.mean(w, axis=0)
        
        lls[i] = np.sum(np.log(clik))

    w = props[-1]*L
    clik = np.sum(w, axis=1)[:, np.newaxis]
    lls[-1] = np.sum(np.log(clik))


    # Return all props
    return props,lls

def get_plots(lls, num_iter, outfile):
    x = range(num_iter+1)
    fig, ax = plt.subplots(1, 1)
    ax.set_title("Loglikelihoods")
    ax.plot(x, lls)
    ax.label_outer()
    ax.set(xlabel='Iteration', ylabel='Loglikelihood')
    plt.savefig(outfile)

def write_ml_matches(results, out_file="ml-matches.txt"):
    with open(out_file, "w") as f:
        for (read,ref, ll) in results:
            f.write(f"{read}\t{ref}\t{ll}\n")


def main():

    args = get_args()

    os.makedirs(args.output_dir, exist_ok=True)
    num_iter = args.num_iter

    # Getting read index position maps
    with open('read_idxs.json') as f:
        read_idx = json.load(f)
        read_idx_rev = {v: k for k, v in read_idx.items()}
    # Getting ref index position maps
    with open('ref_idxs.json') as f:
        ref_idx = json.load(f)
        ref_idx_rev = {v: k for k, v in ref_idx.items()}

    num_dists = len(ref_idx)
    ll_array = np.load("ll_array.npy")
    ll_array[ll_array==0] = ll_array.min()

    x = np.load("ll_array.npy")
    x[x==0] = -np.inf

    iter_props,lls = base_em(np.exp(ll_array), num_iter)
    em_props = iter_props[-1]

    w = em_props*np.exp(ll_array)
    clik = np.sum(w, axis=1)[:, np.newaxis]
    w /= clik

    results = [(read_idx_rev[n], ref_idx_rev[i], f"{ll_array[n,i]:.5f}") for n,i in enumerate(w.argmax(axis=1))]

    get_plots(lls, num_iter, os.path.join(args.output_dir, "em-iterations.png"))

    write_ml_matches(results, os.path.join(args.output_dir, "ml-matches.txt"))

if __name__=="__main__":
    main()