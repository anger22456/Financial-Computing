#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 15:35:57 2022

@author: anger22456
"""
import numpy as np
import math
from scipy.stats import norm


import tkinter as tk
from tkinter.ttk import *


def option_value(s0, r, q, var, t, k1, k2, k3, k4):
     try:

        d1 = (np.log(s0/k1) + (r-q+(var**2)/2)*t)/(var*pow(t, 0.5))
        d2 = (np.log(s0/k2) + (r-q+(var**2)/2)*t)/(var*pow(t, 0.5))
        d3 = (np.log(s0/k1) + (r-q-(var**2)/2)*t)/(var*pow(t, 0.5))
        d4 = (np.log(s0/k2) + (r-q-(var**2)/2)*t)/(var*pow(t, 0.5))
        d5 = (np.log(s0/k3) + (r-q-(var**2)/2)*t)/(var*pow(t, 0.5))
        d6 = (np.log(s0/k4) + (r-q-(var**2)/2)*t)/(var*pow(t, 0.5))
        d7 = (np.log(s0/k3) + (r-q+(var**2)/2)*t)/(var*pow(t, 0.5))
        d8 = (np.log(s0/k4) + (r-q+(var**2)/2)*t)/(var*pow(t, 0.5))
            
        payoff_1 = s0*np.exp((r-q)*t)*(norm.cdf(d1) - norm.cdf(d2))
        payoff_2 = -1*k1*(norm.cdf(d3) - norm.cdf(d4))
        payoff_3 = (k2 - k1)*(norm.cdf(d4) - norm.cdf(d5))
        payoff_4 = (k4*(k2-k1)/(k4-k3))*(norm.cdf(d5) - norm.cdf(d6))
        payoff_5 = -1*((k2-k1)/(k4-k3))*(norm.cdf(d7) - norm.cdf(d8))*np.exp((r-q)*t)*s0
        
        c0 = (payoff_1+payoff_2+payoff_3+payoff_4+payoff_5)*np.exp(-r*t)
        return c0
    
     except ZeroDivisionError:
         d1=d2=d3=d4=d5=d6=d7=d8=0

def monte_carlo(s0, r, q, var, t, k1, k2, k3, k4):
    
    mean = np.log(s0)+(r-q-(var**2)/2)*t
    std = var*pow(t, 0.5)
    
    values_mean = []
    for i in range(20):
        v_sum = 0
        lnSt_simulation = np.random.normal(mean, std, 10000).tolist()
        for j in range(len(lnSt_simulation)):
            st = np.exp(lnSt_simulation[j])
            if st < k1:
                v = 0
            elif (st > k1) & (st <= k2):
                v = np.exp(-r*t)*(st-k1)
            elif (st > k2) & (st <= k3):
                v = np.exp(-r*t)*(k2-k1)
            elif (st > k3) & (st <= k4):
                v = np.exp(-r*t)*((k2-k1)/(k4-k3))*(k4-st)
            else:
                v = 0
            
            v_sum += v
        
        values_mean.append(v_sum/10000)
    
    Avg = np.mean(values_mean)
    STD = np.std(values_mean)
    upper = Avg+2*STD
    lower = Avg-2*STD
    
    return Avg, STD, upper, lower
        


class HW1():

    
    def __init__(self):
        
        self.window = tk.Tk()
        self.window.title('HW1')
        self.window.geometry('400x600')
        self.window.configure(background = 'white')
        
        Label(self.window, text = 'Basic Requirement').pack(side = tk.TOP)
        
        #S0
        self.s0_frame = tk.Frame(self.window)
        self.s0_frame.pack(side = tk.TOP)
        
        Label(self.s0_frame, text = 's0:').pack(side = tk.LEFT)
        self.s0 = tk.DoubleVar()
        Entry(self.s0_frame, textvariable = self.s0, state='normal').pack(side = tk.LEFT)
        
        
        
        #r
        self.r_frame = tk.Frame(self.window)
        self.r_frame.pack(side = tk.TOP)
        
        Label(self.r_frame, text = 'r:').pack(side = tk.LEFT)
        self.r = tk.DoubleVar()
        Entry(self.r_frame, textvariable = self.r, state='normal').pack(side = tk.LEFT)
        
        
        
        #q
        self.q_frame = tk.Frame(self.window)
        self.q_frame.pack(side = tk.TOP)
        
        Label(self.q_frame, text = 'q:').pack(side = tk.LEFT)
        self.q = tk.DoubleVar()
        Entry(self.q_frame, textvariable = self.q, state='normal').pack(side = tk.LEFT)
        
        
        
        #var
        self.var_frame = tk.Frame(self.window)
        self.var_frame.pack(side = tk.TOP)
        
        Label(self.var_frame, text = 'var:').pack(side = tk.LEFT)
        self.var = tk.DoubleVar()
        Entry(self.var_frame, textvariable = self.var, state='normal').pack(side = tk.LEFT)
        
        
        
        #t
        self.t_frame = tk.Frame(self.window)
        self.t_frame.pack(side = tk.TOP)
        
        Label(self.t_frame, text = 't:').pack(side = tk.LEFT)
        self.t = tk.DoubleVar()
        Entry(self.t_frame, textvariable = self.t, state='normal').pack(side = tk.LEFT)
        
        
        
        #k1
        self.k1_frame = tk.Frame(self.window)
        self.k1_frame.pack(side = tk.TOP)
        
        Label(self.k1_frame, text = 'k1:').pack(side = tk.LEFT)
        self.k1 = tk.DoubleVar()
        Entry(self.k1_frame, textvariable = self.k1, state='normal').pack(side = tk.LEFT)
        
        
        
        #k2
        self.k2_frame = tk.Frame(self.window)
        self.k2_frame.pack(side = tk.TOP)
        
        Label(self.k2_frame, text = 'k2:').pack(side = tk.LEFT)
        self.k2 = tk.DoubleVar()
        Entry(self.k2_frame, textvariable = self.k2, state='normal').pack(side = tk.LEFT)
        
        
        
        #k3
        self.k3_frame = tk.Frame(self.window)
        self.k3_frame.pack(side = tk.TOP)
        
        Label(self.k3_frame, text = 'k3:').pack(side = tk.LEFT)
        self.k3 = tk.DoubleVar()
        Entry(self.k3_frame, textvariable = self.k3, state='normal').pack(side = tk.LEFT)
        
        
        
        #k4
        self.k4_frame = tk.Frame(self.window)
        self.k4_frame.pack(side = tk.TOP)
        
        Label(self.k4_frame, text = 'k4:').pack(side = tk.LEFT)
        self.k4 = tk.DoubleVar()
        Entry(self.k4_frame, textvariable = self.k4, state='normal').pack(side = tk.LEFT)
        
        
        
        #button
        self.btn_frame = tk.Frame(self.window)
        self.btn_frame.pack(side = tk.TOP)
        
        Button(self.btn_frame, text = 'Close Form', command = self.basic).pack(side = tk.LEFT)
        Button(self.btn_frame, text = 'MonteCarlo Simulation', command = self.bonus).pack(side = tk.LEFT)
        
        
        #Result
        self.result_frame = tk.Frame(self.window)
        self.result_frame.pack(side = tk.TOP)
        self.result_1 = tk.DoubleVar()
        self.result_1.set(0)
        Label(self.result_frame, text = 'Option Value:').pack(side = tk.LEFT)
        Label(self.result_frame, textvariable = self.result_1).pack(side = tk.TOP)
        
        
        self.monte_frame = tk.Frame(self.window)
        self.monte_frame.pack(side = tk.TOP)
        Label(self.monte_frame, text = '---Monte Carlo Simulation---').pack(side = tk.TOP)
        
        self.bonus_frame = tk.Frame(self.window)
        self.bonus_frame.pack(side = tk.TOP)
        self.bonus_mean = tk.DoubleVar()
        self.bonus_mean.set(0)
        self.bonus_std = tk.DoubleVar()
        self.bonus_std.set(0)
        self.bonus_upper = tk.DoubleVar()
        self.bonus_upper.set(0)
        self.bonus_lower = tk.DoubleVar()
        self.bonus_lower.set(0)
        
        
        Label(self.bonus_frame, text = 'Mean:').pack(side = tk.TOP)
        Label(self.bonus_frame, textvariable = self.bonus_mean).pack()
        Label(self.bonus_frame, text = 'STD:').pack(side = tk.TOP)
        Label(self.bonus_frame, textvariable = self.bonus_std).pack()
        Label(self.bonus_frame, text = 'Upper Bound:').pack(side = tk.TOP)
        Label(self.bonus_frame, textvariable = self.bonus_upper).pack()
        Label(self.bonus_frame, text = 'Lower Bound:').pack(side = tk.TOP)
        Label(self.bonus_frame, textvariable = self.bonus_lower).pack()
        
        
        self.window.mainloop()
        
    def basic(self):
        
        
        #self.result_1.set(self.result_1.get()+1)
        
        s0 = self.s0.get()
        r = self.r.get()
        q = self.q.get()
        var = self.var.get()
        t = self.t.get()
        k1 = self.k1.get()
        k2 = self.k2.get()
        k3 = self.k3.get()
        k4 = self.k4.get()
        
            
        c0 = option_value(s0, r, q, var, t, k1, k2, k3, k4)
        self.result_1.set(c0)
    
    def bonus(self):
        s0 = self.s0.get()
        r = self.r.get()
        q = self.q.get()
        var = self.var.get()
        t = self.t.get()
        k1 = self.k1.get()
        k2 = self.k2.get()
        k3 = self.k3.get()
        k4 = self.k4.get()
        
        mean, std, upper, lower = monte_carlo(s0, r, q, var, t, k1, k2, k3, k4)
        self.bonus_mean.set(mean)
        self.bonus_std.set(std)
        self.bonus_upper.set(upper)
        self.bonus_lower.set(lower)
        
            
    #option_value(100, 0.05, 0.02, 0.5, 0.4, 90, 98, 102, 104)

if __name__ == '__main__':
    HW1()
        

#option_value(100, 0.05, 0.02, 0.5, 0.4, 90, 98, 102, 104)
