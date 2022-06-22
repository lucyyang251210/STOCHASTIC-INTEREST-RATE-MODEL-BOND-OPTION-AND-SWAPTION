import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
# from tqdm import tqdm  #### loop
import plotly
import plotly.graph_objects as go
import scipy.stats as stats
from sklearn.neighbors import KernelDensity
from scipy.stats import norm, gaussian_kde
import seaborn as sns
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import math

#np.random.seed(0) little curve
#np.random.seed(1) skew

#np.random.seed(5)
np.random.seed(6) #nice skew
# define parameters
r_0 = 0.02
alpha = 3
sigma = 0.01

theta_0 = 0.03
beta = 1
phi = 0.05
eta = 0.005

nsims=1000
t = np.linspace (0,6,601)


def simtheta(t, beta, phi, eta, nsims):
    ndt = len(t)
    dt = t[1] - t[0]
    sqrt_dt = np.sqrt(dt)

    theta_t = np.zeros((nsims, ndt))
    theta_t[:, 0] = theta_0

    for i in range(len(t) - 1):
        theta_t[:, i + 1] = theta_t[:, i] + beta * (phi - theta_t[:, i]) * dt + eta * sqrt_dt * np.random.randn(nsims)

    return theta_t


theta = simtheta(t,beta,phi,eta,nsims)


def simr(t, alpha, theta, sigma, nsims):
    ndt = len(t)
    dt = t[1] - t[0]
    sqrt_dt = np.sqrt(dt)

    rt = np.zeros((nsims, ndt))
    rt[:, 0] = r_0

    for i in range(len(t) - 1):
        rt[:, i + 1] = rt[:, i] + alpha * (theta[:, i] - rt[:, i]) * dt + sigma * sqrt_dt * np.random.randn(nsims)

    return rt


rt = simr(t, alpha, theta, sigma, nsims)


def B(alpha, t):
    B_t = (1 - np.exp(-alpha * t)) / alpha
    return B_t


def C(alpha, beta, t, B_t):
    C_t = (alpha / (alpha - beta)) * ((1 - np.exp(-beta * t)) / beta - B_t)
    return C_t


def At(alpha, beta, sigma, phi, eta, t, B_t, C_t):
    O1 = t / beta ** 2 - (2 * (B_t + C_t)) / beta ** 2 + (1 - np.exp(-2 * alpha * t)) / (
                2 * alpha * (alpha - beta) ** 2)

    O21 = -(2 * alpha / (beta * (alpha - beta) ** 2)) * (1 - np.exp(-(alpha + beta) * t)) / (alpha + beta)
    O22 = ((alpha / (beta * (alpha - beta))) ** 2) * ((1 - np.exp(-2 * beta * t)) / (2 * beta))
    O2 = O21 + O22

    A_t = (B_t - t) * (phi - sigma ** 2 / (2 * alpha ** 2)) + C_t * phi - (sigma * B_t) ** 2 / (4 * alpha) + 0.5 * (
                eta ** 2) * (O1 + O2)
    return A_t


def BondPriceSim(rt, theta, t, T, dt):
    delt = T-t
    Bdt = B(alpha, delt)
    Cdt = C(alpha, beta, delt, Bdt)
    Adt = At(alpha, beta, sigma, phi, eta, delt, Bdt, Cdt)
    bondP = np.zeros(1000)
    for i in range(1000):
        bondP[i] = np.exp(Adt-Bdt*rt[i, int(t/dt)]-Cdt*theta[i, int(t/dt)])
    return bondP


def getP0(r_rn, dt):
        
    #bond_valT0 = []
    # calculate for r_tao0
    tau = np.linspace(3, 6,13)
    bond_valT0 = np.zeros((len(tau), r_rn.shape[1]))
    
    i = 0
    while i < len(tau): 
        #for tau_val in tau:
        #print(tau_val)
        #every bon_val is a list of simulated values for bond price at tau
        # every row 
        bond_val = BondPriceSim(rt, theta, 0, tau[i], dt)
        bond_valT0[i,] = bond_val
        i += 1
    
    return bond_valT0


def getSwapRateT0(bond_valT0):
    part = bond_valT0[0,:] - bond_valT0[-1,:]
    a0 = sum(bond_valT0[1:,:]) * 0.25
    
    S0 = part/a0
    return(a0,S0)


def getP3(r_rn, dt):
    tau = np.linspace(3, 6,13)
    bond_valT3 = np.zeros((len(tau), r_rn.shape[1]))
    
    i = 0
    while i < len(tau): 
    #for tau_val in tau:
        #print(tau_val)
        #every bon_val is a list of simulated values for bond price at tau
        # every row 
        bond_val = BondPriceSim(rt, theta, 3, tau[i], dt)
        bond_valT3[i,] = bond_val
        i += 1
    return (bond_valT3)


def getSwapRateT3(bond_valT3):    
    # bond_valT3 each row is an estimate of value P3 at tao
    
    # n: number of columns i.e. mc samples 
    
    #n = bond_valT3.shape[1]
    part = bond_valT3[0,:] - bond_valT3[-1,:]
    #print("part" ,part)
    a3 = sum(bond_valT3[1:, :]) * 0.25
    #print("a3" ,a3)
    s3 = part/a3
    
    return(a3,s3)

# def getSwapRateT3_v2(bond_valT3):    
#     # bond_valT3 each row is an estimate of value P3 at tao
    
#     # n: number of columns i.e. mc samples 
    
#     n = bond_valT3.shape[1]
#     part = no,zeros(n)
#     while j < n:
#         part[j] = bond_valT3[0,j] - bond_valT3[-1,j]
#     #print("part" ,part)
#         a3[j] = sum(bond_valT3[1:, j) * 0.25
#     #print("a3" ,a3)
#         s3[j] = part[j]/a3[j]
#     return(a3,s3)
    
# def OptionPrice(r, T1, T2, dt, ai):
#     K = np.zeros(r.shape[1])
#     bondP = np.zeros(r.shape[1])
#     payoff = np.zeros(r.shape[1])
#     optionP = np.zeros(r.shape[1])
#     for i in range(r.shape[1]):
#         K[i] = np.exp(-r[1:int(T2/dt), i].sum() * dt) / np.exp(-r[1:int(T1/dt), i].sum() * dt)
#         bondP[i] = np.exp(-r[int(T1/dt)+1:int(T2/dt), i].sum() * dt)
#         print(bondP[i])
#         payoff[i] = max(0, bondP[i] - ai * K[i])
#         optionP[i] = np.exp(-r[1:int(T1/dt), i].sum() * dt) * payoff[i]
#     return optionP.mean()


def SwapOptionPrice(r_rn,a0,a3,s3,K):
    #T1 start
    #a3, s3 = getSwapRateT3(bond_valT3)
    #A = np.linspace(0.95, 1.05, 11)
    #K_list = A * K
    
    #bondP = np.zeros(r.shape[1])
    payoff = np.zeros(r_rn.shape[1])
    v0 = np.zeros(r_rn.shape[1])
    dt = 0.01
    
    i = 0
    while i < r_rn.shape[1]:
        #print(bondP[i])
        payoff[i] = max(0, s3[i] - K)
        #print(payoff[i])

        v0[i] = np.exp(-r_rn[0: 299, i].sum() * dt) * payoff[i] * a3[i]
        #v0[i] =  payoff[i]
        #print(v0[i])
        i+=1
    #v0 = v0 * np.exp(-r_rn[1:int(3/dt), i].sum() * dt)
    return (np.mean(v0))

def SwapOptionPriceQb(r_rn,a0,a3,s3,K):
    #T1 start
    #a3, s3 = getSwapRateT3(bond_valT3)
    #A = np.linspace(0.95, 1.05, 11)
    #K_list = A * K
    
    #bondP = np.zeros(r.shape[1])
    payoff = np.zeros(r_rn.shape[1])
    v0 = np.zeros(r_rn.shape[1])
    dt = 0.01
    
    i = 0
    while i < r_rn.shape[1]:
        #print(bondP[i])
        payoff[i] = max(0, s3[i] - K)
        #print(payoff[i])

        v0[i] = payoff[i]
        #v0[i] =  payoff[i]
        #print(v0[i])
        i+=1
    #v0 = v0 * np.exp(-r_rn[1:int(3/dt), i].sum() * dt)
    return (np.mean(v0)*a0)
    
        

def LSM_SwaptionValue(a0,S0,K,sigma):

    # i < n mc samples
       
    # 6-3 or 6-0?
    omega = sigma * math.sqrt(3 - 0)
    # note here r =0  q = 0 ,sigma should be fid as implied vol
    
    d1 = (np.log(S0/K) + (0.5 * omega**2)) / omega
    #print(S0/K)
    d2 = (np.log(S0/K) - (0.5 * omega**2)) / omega
    
    # note we know S = K
    #LSM_V = S * norm.cdf(d1) - norm.cdf(d2)
    LSM_V = (S0 * norm.cdf(d1) - K * norm.cdf(d2)) * a0
    return LSM_V



def brutalforce(simSwapP, a0,S0,K):
    #t0_swap, sim_swapV = getSwaptionValue(n, T, m, dt, K)
    # simulated_swaption_price 
   
    #simSwapP = SwapOptionPrice(r, T1, T2, dt, ai)
    
    # possible volatility 
    volatility_candidates = np.arange(0.0001,1.9999,0.0001)
    # initialize to hold price difference 
    price_differences = np.zeros_like(volatility_candidates) 

    for i in range(len(volatility_candidates)):
    
        candidate = volatility_candidates[i]
    
        price_differences[i] = simSwapP - LSM_SwaptionValue(a0,S0,K,candidate)
        
    idx = np.argmin(abs(price_differences))
    implied_volatility = volatility_candidates[idx]
    
    return(implied_volatility)
    

def getImpVol(v0,s0,a0):
    
    value =  (v0 / (a0 * s0) +1) * 0.5

    #value = numerator / denom
    omega = 2 * norm.ppf(value)
    
    impVol = math.sqrt(omega/3)
    return impVol

def getImpVol_test(v0,s0,a0, scale):
    
    volatility_candidates = np.arange(0.0001,1.9999,0.0001)
    
    diff = np.zeros_like(volatility_candidates) 
    
    for i in range(len(volatility_candidates)):
    
        omega = volatility_candidates[i] * np.sqrt(3)
        
        sim_LSM_value = v0 / (a0 * s0)
        d1 = (np.log(1/scale) + 0.5 * omega ** 2) / omega
        d2 = (np.log(1/scale) - 0.5 * omega ** 2) / omega
        
        LSM_val = norm.cdf(d1) -  scale * norm.cdf(d2)
    
        diff[i] = sim_LSM_value - LSM_val
        
    idx = np.argmin(abs(diff))
    implied_volatility = volatility_candidates[idx]
    return(implied_volatility)





#------ to run


if __name__ == '__main__':
    
    T1 = 3
    T2 = 6
    dt = 0.01
    mc = 1000

    # below is for interest rate simulation
    rt = simr(t, alpha, theta, sigma, nsims)
    r_rn = np.transpose(rt)
    #bond = BondPrice(r_rn, T1, T2, dt)
    
    # below is to simulate a0, S0
    bond_valT0 = getP0(r_rn, dt)
    a0, S0 = getSwapRateT0(bond_valT0)
    
    a0 = np.mean(a0)
    S0 = np.mean(S0)
    #print(a0)
    #print(S0)
    
    # set inital K = S0


    # this is to simulate bond price P_3(0), P_3(1) .....
    bond_valT3 = getP3(r_rn, dt)
    
    #print("T3", bond_valT3)
    #print(bond_valT3)
    #bond_valT3_needed = bond_valT3[1:,]
    
    a3, s3 = getSwapRateT3(bond_valT3)
    
    #K = np.mean(s3)
    K = np.mean(S0)
    A = np.linspace(0.95, 1.05, 11)
    K_list = A * K
    

    
    ii = 0
    
    sim_v0 = []
    imp_vol_brutal = []
    imp_vol_a = []
    imp_vol_b = []
    

    while ii < len(A):
        
        sim_v = SwapOptionPriceQb(r_rn,a0, a3,s3,K_list[ii])

        sim_v0.append(sim_v)
        
        imp = getImpVol_test(sim_v,S0,a0, A[ii])
        imp_vol_brutal.append(imp)
        
        imp_a = getImpVol(sim_v,S0,a0)
        imp_vol_a.append(imp_a)
        
        imp_b = brutalforce(sim_v, a0,S0,K_list[ii])
        imp_vol_b.append(imp_b)
        
        
        #sim_v0.append(SwapOptionPrice(r_rn,a3,s3,K_list[ii]))
        ii +=1
        
    plt.plot(K_list,imp_vol_brutal , label = "Implied Volatility Curve")
    #plt.plot(K_list, imp_vol_a, label = "imp_vol_analytical")
    #plt.plot(K_list, imp_vol_b, label = "imp_vol_brrrutal")
    #plt.plot(imp_vol_brutal, sim_v0)
    #plt.axvline(x=np.mean(s3))
    #plt.plot(np.mean(s3), 0, marker = '*', clip_on =False)
    #plt.vlines(x=np.mean(s3), ymin = np.min(imp_vol_brutal), ymax=np.min(imp_vol_brutal)+0.0004, colors='green', ls='-', lw=2, label="K = S_T")
    #plt.vlines(x=S0, ymin = np.min(imp_vol_brutal), ymax=np.min(imp_vol_brutal)+0.0004, colors='r', ls='-', lw=2, label="K = S_0")
    plt.xlabel("Strike")
    plt.ylabel("Implied Volatility")
    plt.title("Implied Volatility Curve against Strike Price")
    plt.legend()
    plt.show()
    print(sim_v0)
    
    plt.plot(imp_vol_brutal, sim_v0, label = "Implied Volatility Curve")
    #plt.plot(K_list, imp_vol_a, label = "imp_vol_analytical")
    #plt.plot(K_list, imp_vol_b, label = "imp_vol_brrrutal")
    #plt.plot(imp_vol_brutal, sim_v0)
    plt.xlabel("Implied Volatility")
    plt.ylabel("Simulated Swaption Price")
    plt.title("Simulated Swaption Price against Implied Volatility")
    plt.show()
    print(sim_v0)
    
    plt.plot(K_list,sim_v0)
    #plt.plot(K_list, imp_vol_a, label = "imp_vol_analytical")
    #plt.plot(K_list, imp_vol_b, label = "imp_vol_brrrutal")
    #plt.plot(imp_vol_brutal, sim_v0)
    plt.ylabel("Simulated swaption prices")
    plt.xlabel("Strikes")
    plt.title("Simulated swaption price changes with Strikes ")
    plt.show()
    print(sim_v0)
        
        

    
