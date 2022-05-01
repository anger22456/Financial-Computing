#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 15:28:12 2022

@author: anger22456
"""

import numpy as np
import math
from scipy.stats import norm


def black_scholes(s0, r, q, var, t, k):
    try:
        d1 = (np.log(s0/k) + (r - q + (var**2)/2)*t)/(var*pow(t, 0.5))
        d2 = (np.log(s0/k) + (r - q - (var**2)/2)*t)/(var*pow(t, 0.5))
        
        call_payoff = s0*np.exp((r-q) * t)*norm.cdf(d1)-k*norm.cdf(d2)
        put_payoff = k*(1-norm.cdf(d2)) - s0*np.exp((r-q) * t)*(1-norm.cdf(d1))
        c0 = call_payoff*np.exp(-r * t)
        p0 = put_payoff*np.exp(-r * t)
        
        print('European Call:{: .4f}\n'.format(c0)+\
              'European Put:{: .4f}\n'.format(p0))
        
        #return output
    
    
    except ZeroDivisionError:
        d1 = d2 = 0

#S:Simulation, R: Repetitions
def call_monte_carlo(s0, r, q, var, t, k, S, R):
    mean = np.log(s0)+(r-q-(var**2)/2)*t
    std = var*pow(t, 0.5)
    
    values_mean = []
    for i in range(R):
        v_sum = 0
        lnSt_simulation = np.random.normal(mean, std, S).tolist()
        for j in range(len(lnSt_simulation)):
            st = np.exp(lnSt_simulation[j])
            if st < k:
                v = 0
            elif st>=k:
                v = np.exp(-r*t)*(st-k)
            
            v_sum += v
        
        values_mean.append(v_sum/10000)
    
    Avg = np.mean(values_mean)
    STD = np.std(values_mean)
    upper = Avg+2*STD
    lower = Avg-2*STD
    
    print('--European Call--\nAvg:{: .4f}\n'.format(Avg)+\
          'STD:{: .4f}\n'.format(STD)+\
          'Upper:{: .4f}\n'.format(upper)+\
          'Lower:{: .4f}\n'.format(lower))
    #return Avg, STD, upper, lower

def put_monte_carlo(s0, r, q, var, t, k, S, R):
    mean = np.log(s0)+(r-q-(var**2)/2)*t
    std = var*pow(t, 0.5)
    
    values_mean = []
    for i in range(R):
        v_sum = 0
        lnSt_simulation = np.random.normal(mean, std, S).tolist()
        for j in range(len(lnSt_simulation)):
            st = np.exp(lnSt_simulation[j])
            if st > k:
                v = 0
            elif st <= k:
                v = np.exp(-r*t)*(k-st)
            
            v_sum += v
        
        values_mean.append(v_sum/10000)
    
    Avg = np.mean(values_mean)
    STD = np.std(values_mean)
    upper = Avg+2*STD
    lower = Avg-2*STD
    
    print('--European Put--\nAvg:{: .4f}\n'.format(Avg)+\
          'STD:{: .4f}\n'.format(STD)+\
          'Upper:{: .4f}\n'.format(upper)+\
          'Lower:{: .4f}\n'.format(lower))


def CRR_2D(s0, r, q, var, t, k, n):
    #先畫出二維矩陣
    Eur_call = [[0 for x in range(n+1)] for y in range(n+1)]
    Amer_call = [[0 for x in range(n+1)] for y in range(n+1)]
    Eur_put = [[0 for x in range(n+1)] for y in range(n+1)]
    Amer_put = [[0 for x in range(n+1)] for y in range(n+1)]
    St = [[0 for x in range(n+1)] for y in range(n+1)]
    
    #按照公式解u, d, p
    u = np.exp(var * pow(t/n, 0.5))
    d = np.exp(-var * pow(t/n, 0.5))
    p = (np.exp((r-q)*(t/n))-d)/(u - d)
    
    #每一個時段的股價
    for i in range(n+1):
        for j in range(i+1):
            St[j][i] = s0*(u**(i-j))*(d**j)
    
    
    #先填上最後一期的選擇權價值
    for j in range(n+1):
        Eur_call[j][n] = max(s0*(u**(n-j))*(d**j)-k, 0)
        Amer_call[j][n] = max(s0*(u**(n-j))*(d**j)-k, 0)
        Eur_put[j][n] = max(k - s0*(u**(n-j))*(d**j), 0)
        Amer_put[j][n] = max(k - s0*(u**(n-j))*(d**j), 0)
    
    for i in range(n-1, -1, -1):
        for j in range(i+1):
            Eur_call[j][i] = np.exp(-r*(t/n))*(p*Eur_call[j][i+1] + (1-p)*Eur_call[j+1][i+1])
            Amer_call[j][i] = max(np.exp(-r*(t/n))*(p*Amer_call[j][i+1] + (1-p)*Amer_call[j+1][i+1]), St[j][i]-k)
            Eur_put[j][i] = np.exp(-r*(t/n))*(p*Eur_put[j][i+1] + (1-p)*Eur_put[j+1][i+1])
            Amer_put[j][i] = max(np.exp(-r*(t/n))*(p*Amer_put[j][i+1] + (1-p)*Amer_put[j+1][i+1]), k-St[j][i])
    
    
    Eur_call_value = Eur_call[0][0]
    Amer_call_value = Amer_call[0][0]
    Eur_put_value = Eur_put[0][0]
    Amer_put_value = Amer_put[0][0]
    
    
    print('--European Call Price :{: .4f}\n'.format(Eur_call_value)+\
          '--European Put Price :{: .4f}\n'.format(Eur_put_value)+\
          '--American Call Price :{: .4f}\n'.format(Amer_call_value)+\
          '--American Put Price :{: .4f}\n'.format(Amer_put_value))
    
    #return Eur_call_value, Amer_call_value, Eur_put_value, Amer_put_value

def CRR_1D(s0, r, q, var, t, k, n):
    Eur_call = [0 for x in range(n+1)]
    Amer_call = [0 for x in range(n+1)]
    Eur_put = [0 for x in range(n+1)]
    Amer_put = [0 for x in range(n+1)]
    St = [[0 for x in range(n+1)] for y in range(n+1)]
    
    
    #按照公式解u, d, p
    u = np.exp(var * pow(t/n, 0.5))
    d = np.exp(-var * pow(t/n, 0.5))
    p = (np.exp((r-q)*(t/n))-d)/(u - d)
    
    #每一個時段的股價
    for i in range(n+1):
        for j in range(i+1):
            St[j][i] = s0*(u**(i-j))*(d**j)
    
    #先填上最後一期的選擇權價值
    for j in range(n+1):
        Eur_call[j] = max(s0*(u**(n-j))*(d**j)-k, 0)
        Amer_call[j] = max(s0*(u**(n-j))*(d**j)-k, 0)
        Eur_put[j] = max(k - s0*(u**(n-j))*(d**j), 0)
        Amer_put[j] = max(k - s0*(u**(n-j))*(d**j), 0)
    
    for i in range(n-1, -1, -1):
        for j in range(i+1):
            Eur_call[j] = np.exp(-r*(t/n))*(p*Eur_call[j] + (1-p)*Eur_call[j+1])
            Amer_call[j] = max(np.exp(-r*(t/n))*(p*Amer_call[j] + (1-p)*Amer_call[j+1]), St[j][i]-k)
            Eur_put[j] = np.exp(-r*(t/n))*(p*Eur_put[j] + (1-p)*Eur_put[j+1])
            Amer_put[j] = max(np.exp(-r*(t/n))*(p*Amer_put[j] + (1-p)*Amer_put[j+1]), k-St[j][i])
    
    Eur_call_value = Eur_call[0]
    Amer_call_value = Amer_call[0]
    Eur_put_value = Eur_put[0]
    Amer_put_value = Amer_put[0]
    
    
    print('--European Call Price :{: .4f}\n'.format(Eur_call_value)+\
          '--European Put Price :{: .4f}\n'.format(Eur_put_value)+\
          '--American Call Price :{: .4f}\n'.format(Amer_call_value)+\
          '--American Put Price :{: .4f}\n'.format(Amer_put_value))
    #return Eur_call_value, Amer_call_value, Eur_put_value, Amer_put_value


#bonus2

def binnomial(n, j, p):
    
    temp = 0
    
    for i in range(n, n-j, -1):
        temp += np.log(i)

    for i in range(1, j+1):
        temp -= np.log(i)
    
    temp = temp + (n-j)*np.log(p) + j*np.log(1-p)
    
    return np.exp(temp)
    


def Combinatorial(s0, r, q, var, t, k, n):
    #按照公式解u, d, p
    u = np.exp(var * pow(t/n, 0.5))
    d = np.exp(-var * pow(t/n, 0.5))
    p = (np.exp((r-q)*(t/n))-d)/(u - d)
    
    call_price = 0
    put_price = 0
    
    for j in range(n+1):
        call_price += binnomial(n, j, p)*max(s0* (u**(n-j)) * (d**j) - k, 0)
        put_price +=  binnomial(n, j, p)*max(k - s0* (u**(n-j)) * (d**j), 0)
    
    c0 = call_price * np.exp(-r*t)
    p0 = put_price * np.exp(-r*t)
    
    print('--European Call Price :{: .4f}\n'.format(c0)+\
          '--European Put Price :{: .4f}\n'.format(p0))
    #return c0, p0


def comprehensive(s0, k, r, q, var, t, S, R, n):
    
    print('-----------------------------------')
    print('Black Scholes Formulas\n')
    black_scholes(s0, r, q, var, t, k)
    print('-----------------------------------')
    print('Monte Carlo Simulation\n')
    call_monte_carlo(s0, r, q, var, t, k, S, R)
    put_monte_carlo(s0, r, q, var, t, k, S, R)
    print('-----------------------------------')
    print('CRR 2D\n')
    CRR_2D(s0, r, q, var, t, k, n)
    print('-----------------------------------')
    print('CRR 1D\n')
    CRR_1D(s0, r, q, var, t, k, n)
    print('-----------------------------------')
    print('Combinatorial Method\n')
    Combinatorial(s0, r, q, var, t, k, n)
    print('-----------------------------------')
    
    
comprehensive(50, 50, 0.1, 0.05, 0.4, 0.5, 10000, 20, 500)

##測試用##
CRR_2D(50, 0.1, 0.05, 0.4, 0.5, 50, 1000)
CRR_1D(50, 0.1, 0.05, 0.4, 0.5, 50, 1000)
black_scholes(50, 0.1, 0.05, 0.4, 0.5, 50)
call_monte_carlo(50, 0.1, 0.05, 0.4, 0.5, 50, 10000, 20)
put_monte_carlo(50, 0.1, 0.05, 0.4, 0.5, 50, 10000, 20)
Combinatorial(50, 0.1, 0.05, 0.4, 0.5, 50, 1000)

