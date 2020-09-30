from scipy.stats import norm, t
import numpy as np
from math import sqrt
from random import randint

def Zt_test(two_tail: bool,upper: bool,mean: float, std: float,mu_0: float, n: int,alpha=0.05,sample=None):
    def Z_test(two_tail,upper,mean,std,n,alpha):
        s = std/sqrt(n)
        Z = (mean - mu_0)/s
        if two_tail:
            Z = (-2)*abs(Z)
            return (norm.cdf(Z), norm.cdf(Z) <= alpha)
        else:
            if upper:
                return (norm.cdf(Z), norm.cdf(Z) <= alpha)
            else:
                Z = (-1)*Z
                return (norm.cdf(Z), norm.cdf(Z) <= alpha)
    def t_test(two_tail,upper,mean,std,n,alpha):
        s = std/sqrt(n)
        t_stat = (mean - mu_0)/s
        d = n - 1
        if two_tail:
            t_stat = (-2)*abs(t_stat)
            return (t.cdf(t_stat,df=d), t.cdf(t_stat,df=d) <= alpha)
        else:
            if upper:
                return (t.cdf(t_stat,df=d), t.cdf(t_stat,df=d) <= alpha)
            else:
                t_stat = (-1)*t_stat
                return (t.cdf(t_stat,df=d), t.cdf(t_stat,df=d) <= alpha)
    def print_summary(Z,two_tail,upper,mean,std,mu_0,alpha,n,reject,p):
        if Z:
            print('Z Test Summary for comparison of mean against μ_0')
            if two_tail:
                print('H_0 (Null Hypothesis): μ = {0}'.format(mu_0))
                print('H_a (Alternative Hypothesis): μ ≠ {0}'.format(mu_0))
            else:
                if upper:
                    print('H_0 (Null Hypothesis): μ ≥ {0}'.format(mu_0))
                    print('H_a (Alternative Hypothesis): μ < {0}'.format(mu_0))
                else:
                    print('H_0 (Null Hypothesis): μ ≤ {0}'.format(mu_0))
                    print('H_a (Alternative Hypothesis): μ > {0}'.format(mu_0))
            print('Sample mean:',mean)
            print('Sample standard deviation:',std)
            print('Sample size:',n)
            print('Z-score:',(mean-mu_0)/(std/sqrt(n)))
            print('significance level:',round((1-alpha)*100,2),'%')
            print('p-value:',p)
            if reject:
                print('Conclusion: Reject null hypothesis')
            else:
                print('Conclusion: Do not reject null hypothesis')
        else:
            print('t Test Summary for comparison of mean against μ_0')
            if two_tail:
                print('H_0 (Null Hypothesis): μ = {0}'.format(mu_0))
                print('H_a (Alternative Hypothesis): μ ≠ {0}'.format(mu_0))
            else:
                if upper:
                    print('H_0 (Null Hypothesis): μ ≥ {0}'.format(mu_0))
                    print('H_a (Alternative Hypothesis): μ < {0}'.format(mu_0))
                else:
                    print('H_0 (Null Hypothesis): μ ≤ {0}'.format(mu_0))
                    print('H_a (Alternative Hypothesis): μ > {0}'.format(mu_0))
            print('Sample mean:',mean)
            print('Sample standard deviation:',std)
            print('Sample size:',n)
            print('Z-score:',(mean-mu_0)/(std/sqrt(n)))
            print('significance level:',round((1-alpha)*100,2),'%')
            print('p-value:',p)
            if reject:
                print('Conclusion: Reject null hypothesis')
            else:
                print('Conclusion: Do not reject null hypothesis')
    if sample is not None:
        if n > 30:
            result = Z_test(two_tail,upper,np.mean(sample),np.std(sample),len(sample),alpha)
            p = result[0]
            reject = result[1]
            print_summary(True,two_tail,upper,np.mean(sample),np.std(sample),mu_0,alpha,len(sample),reject,p)
        else:
            result = t_test(two_tail,upper,np.mean(sample),np.std(sample),len(sample),alpha)
            p = result[0]
            reject = result[1]
            print_summary(False,two_tail,upper,np.mean(sample),np.std(sample),mu_0,alpha,len(sample),reject,p)
    else:
        if n > 30:
            reject = Z_test(two_tail,upper,mean,std,n,alpha)
            p = result[0]
            reject = result[1]
            print_summary(True,two_tail,upper,mean,std,mu_0,alpha,n,reject,p)
        else:
            reject = t_test(two_tail,upper,mean,std,n,alpha)
            p = result[0]
            reject = result[1]
            print_summary(False,two_tail,upper,mean,std,mu_0,alpha,n,reject,p)

#replace parameters as you wish
mean = randint(-1000,1001)
std = sqrt(randint(2,8000)) #must be positive
sample_size = randint(2,70) # must be > 1
# generate sample from normal distribution with parameters
sample = np.random.normal(loc=mean,scale=std,size=sample_size)
if randint(0,1) == 0:
    two_tail = True
else:
    two_tail = False
if randint(0,1) == 0:
    upper = True
else:
    upper = False
Zt_test(two_tail,upper,mean,std,mean+randint(-5,5),sample_size,sample=sample)
