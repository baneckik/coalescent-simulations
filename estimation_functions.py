import numpy as np
import pandas as pd
import random
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import minimize

def generate_kingman_times(M, theta):
    '''
    Function generates random kingman process for M starting lineages 
    and constant mutation parameter theta.
    Function returnes M-1 coalescence times (cumulative) for generated process.
    '''
    t = []
    for m in range(2, M+1):
        t.append(np.random.exponential( 1/(m*(m-1)/theta) ))
    tau = np.cumsum(list(reversed(t)))
    return list(tau)

def disturb_times_gamma(times, variance):
    return [ random.gammavariate(alpha=t**2/variance, beta=t/variance) for t in times ]

def generate_kingman_intervals(M, theta, lam=1):
    '''
    Function generates random kingman process for M starting lineages 
    and constant mutation parameter theta.
    Function returnes M-1 intervals for coalescence times (cumulative) for generated process.
    '''
    tau = generate_kingman_times(M, theta)
    intervals = [ [random.uniform(0,t), t+np.random.exponential(lam*t) ] for t in tau ]
    return list(intervals)

def estimate_theta_exp(tau_i, M=0, tau=None):
    '''
    Function estimates mutation parameter theta 
    via maximum likelihood method (analitycal formula).
    all the times from tau_i should be < tau
    M>N
    '''
    tau_i = list(tau_i)
    if M==0:
        M = len(tau_i)+1
    
    N = M-len(tau_i)
    tau_i.sort(reverse=False)
    tau_i = [0] + tau_i
    t_i = [ tau_i[i]-tau_i[i-1] for i in range(1, len(tau_i)) ]
    if N==1:
        theta = sum([i*(i-1)*t_i[M-i] for i in range(2,M+1)])/(M-1)
        return theta
    else:
        S_t = tau - sum(t_i)
        theta = ( sum([i*(i-1)*t_i[M-i] for i in range(N+1,M+1)]) + N*(N-1)*S_t )/(M-N)
        return theta

def estimate_theta_exp_from_intervals(intervals):
    '''
    Function estimates mutation parameter theta 
    via maximum likelihood method (analitycal formula).
    all the times from tau_i should be < tau
    M>N
    '''
    
    tau_i = []
    for interval in intervals:
        tau_i += [(interval[1]+interval[0])/2]

    return estimate_theta_exp(tau_i)

    
def estimate_theta_lln(tau_i, M=0, tau=None, kappa0=1, eps=0.0001):
    '''
    Function estimates mutation parameter theta by maximum likelihood via LLN method.
    '''
    
    tau_i = list(tau_i)
    if M==0:
        M = len(tau_i)+1
    N = M-len(tau_i)
    if N==1:
        # Newton-Raphson method
        kappa = kappa0
        g = 2*sum([i/(i+1/kappa0) for i in tau_i])-M+1
        g_prim = 2*sum([i/(1+kappa0*i)**2 for i in tau_i])
        kappa = kappa0-g/g_prim
        if kappa<0:
            kappa = kappa0/2
        while np.abs(kappa-kappa0)>eps:
            kappa0 = kappa
            g = 2*sum([i/(i+1/kappa0) for i in tau_i])-M+1
            g_prim = 2*sum([i/(1+kappa0*i)**2 for i in tau_i])
            kappa = kappa0-g/g_prim
            if kappa<0:
                kappa = kappa0/2
        return M/kappa
    else:
        S_tau = tau - max(tau_i)
        # Newton-Raphson method
        kappa = kappa0
        g = 2*sum([i/(i+1/kappa0) for i in tau_i])+kappa0*tau/(1+kappa0*tau)*(N-1)-(M-N)
        g_prim = 2*sum([i/(1+kappa0*i)**2 for i in tau_i])+(N-1)/(1+kappa0*tau)**2
        kappa = kappa0-g/g_prim
        if kappa<0:
            kappa = kappa0/2
        while np.abs(kappa-kappa0)>eps:
            kappa0 = kappa
            g = 2*sum([i/(i+1/kappa0) for i in tau_i])+kappa0*tau/(1+kappa0*tau)*(N-1)-(M-N)
            g_prim = 2*sum([i/(1+kappa0*i)**2 for i in tau_i])+(N-1)/(1+kappa0*tau)**2
            kappa = kappa0-g/g_prim
            if kappa<0:
                kappa = kappa0/2
        return M/kappa

    
def estimate_theta_lln_from_intervals(intervals, kappa0=1, eps=0.0001):
    '''
    Function estimates mutation parameter theta by maximum likelihood via LLN method.
    '''
    
    M = len(intervals)+1
    
    # Newton-Raphson method
    kappa = kappa0
    g = sum([ kappa0*i[0]/(1+kappa0*i[0])+kappa0*i[1]/(1+kappa0*i[1]) for i in intervals])-M+1
    g_prim = sum([i[0]/(1+kappa0*i[0])**2 + i[1]/(1+kappa0*i[1])**2 for i in intervals])
    kappa = kappa0-g/g_prim
    if kappa<0:
        kappa = kappa0/2
    while np.abs(kappa-kappa0)>eps:
        kappa0 = kappa
        g = sum([ kappa0*i[0]/(1+kappa0*i[0])+kappa0*i[1]/(1+kappa0*i[1]) for i in intervals])-M+1
        g_prim = sum([i[0]/(1+kappa0*i[0])**2 + i[1]/(1+kappa0*i[1])**2 for i in intervals])
        kappa = kappa0-g/g_prim
        if kappa<0:
            kappa = kappa0/2
    return M/kappa

    
def erase_times_from_intervals(times, intervals):
    indexes_to_remove = []
    for i, t in enumerate(times):
        if any([left<t<right for (left,right) in intervals]):
            indexes_to_remove += [i]
    return [t for i, t in enumerate(times) if i not in indexes_to_remove]





