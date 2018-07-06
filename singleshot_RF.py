import numpy as np
import pylab as plt
from numpy import linalg as LA
import scipy.sparse as sparse
def rflt(phi): 
    th=-0.1 #The threshold for the restarting part of the algorithm
    cntmax=5 # cntmax is the max of "cnt" in the restarting part of the algorithm. In other words "cntmax" number of measurements is done and if the estimated phase "mu" is not found to be accurate the algorithm restarts by setting "mu" to initial mu and "sigma" to initial sigma 
    nums=1000 #num of samples drawn from the Gaussian
    nm=100 #a limit that restricts the number of measurements. However the number of measurements could go beyond 100 because if the estimate is not accurate the counter of the "while" loop initializes. 
    tau=12.5 #minimum interaction time
    T2=1300.0 #decoherence time
    muin=0.0 # mean of the initial Gaussian
    mu=muin
    mu0=0.0
    sigma00=0.0
    #sigma10=0.0
    sigmain=np.pi**2/3.0 #variance of the initial Gaussian. The uniform distribution between [-pi,pi] is approximated by a Gaussian with mean zero and this variance
    sigma=sigmain
    t=tau
    theta=0.0
    sigmav=0.0
    muv=0.0
    muf=0.0
    sigmaf=0.0
    muv=[] #an array to keep track of the means
    sigmav=[] #an array to keep trakc of the variances
    sigmav.append(sigmain)
    muv.append(mu)
    cnt=0
    n1=1
    tott=0.0
    while n1 < nm:
        sigma00=0.0
        sigma10=0.0
        mu0=0.0
        tott=tott+t #This calculates total interaction time
        tr = t/T2
        vis=np.exp(-tr**2)
        xr=np.random.uniform(0,1)
        p0=0.5*(1.0+vis*np.cos(t*phi/tau-theta))
        if xr<=p0: #u1 is the measurement result
            u1=0.0 
        else:
            u1=1.0
        n3=0
        phis=np.random.normal(mu,np.sqrt(sigma),nums) #samples drawn from the Gaussian with mean "mu" and variance "sigma"
        pps=0.5*(1.0+(-1.0)**u1*vis*np.cos(t*phis/tau-theta)) # The probability that the drawn samples "phis" gives the measurement result "u1"
        maxpps=np.amax(pps)
        kappae=np.random.uniform(maxpps,1) #this scales the above probability "pps" and the samples are accepted based on pps/kappae. It is a random number between max of p(phis|u1) and 1 (see PRL 117, 010503 (2016)). Setting "kappae" to 1 does not significantly change the final result.
        for n2 in range(1,nums):
            u2=np.random.uniform(0,1)
            spps=pps[n2]/kappae
            if u2<=spps: #accepting the samples based on pps/kappae
                mu0=mu0+phis[n2]    #preparing the mean of the accepted samples
                sigma00=sigma00+phis[n2]**2 #and this prepares the variance of the accepted samples
                #sigma10=sigma10+phisp(n2)**2
                n3=n3+1
       #if n3>1: 
        mu=mu0/n3 #The mean of the accepted samples
        muv.append(mu) #storing the means in the array "muv"
        sigma1=(sigma00-n3*mu**2)/(n3) #The variance of the accepted samples
        #sigma2=(sigma10-n3*mu**2)/(n3-1.0)
        sigma=sigma1 
        sigmav.append(sigma) #Storing the variance to the array "sigmav"
        #else:
        #    print('n3')
        #    print(n3)
        #    mu=muin
        #    sigma=sigmain
        #    n1=1
        ################################################################ Restarting part of the algorithm (see supplement of PRL 117, 010503 (2016))
        #if cnt>=cntmax: 
        dv=np.log10(np.sqrt(sigmav[n1]))-np.log10(np.sqrt(sigmav[n1-1]))
        if dv>=th and cnt>=cntmax:
            pt=0.5*(1.0+np.exp(-(tau/T2)**2)*np.cos(phi-mu)) #This checks if the current estimate "mu" is accurate. For t=tau=12.5 and T2=1300 vis~1 therefore if mu is close to the system phase "phi" then cos(phi-mu)~1 and pt~1. As a result we should get ut=0.
            xrt=np.random.uniform(0,1)
            if xrt<=pt:
                ut=0
            else:
                ut=1
            if ut==0: # Verifies that the estimate i.e., "mu" is accurate 
                cnt=cnt+1
            else: #In this case the estimate is not accurate and the variance "sigma" and the mean "mu" needs to be reset. Here we set the variance to the initial variance "sigmain".
                cnt=0
                sigma=sigmain#100*np.amin(sigmav)
                mu=muin#muv[np.argmin(sigmav)]
                n1=1
        else:
            cnt=cnt+1
        #########################################################
        t1=np.ceil(1.25/np.sqrt(sigma))*tau #This calculates the interaction time for the next measurement based on the current variance "sigma".
        t=np.minimum(t1,T2) #If the interaction time is bigger than the decoherence time we choose t=T2
        theta=mu+np.random.normal(mu,np.sqrt(sigma)) #The controlled phase for the next measurement
        n1=n1+1
    sigmaf=sigma#np.amin(sigmav) 
    #sfl=np.argmin(sigmav)#minloc[sigmav]
    muf=mu#muv[sfl]
    #op=(muf, sigmaf, sfl, tott)
    op=(muf, sigmaf, tott)
    return op
sigv=[]
mse0=[]
dimtot=1000
sphiv=[]
sharp=[]
timev=[]
#np.random.seed(305)
fName = r'C:\HTD\data\singleshotrf.txt'
tf = open(fName, 'w+')
for cntr in range(0,dimtot):
    sphi=np.random.uniform(-np.pi,np.pi)#-2.1662379052008944
    #print('system phase')
    #print(sphi)
    sphiv.append(sphi)
    musig=rflt(sphi)
    #print('mu and sigma')
    #print(musig[0])
    #print(musig[1])
    #int(input('input an integer'))
    sigv.append(musig[1])
    mse0.append((musig[0]-sphi)**2)
    sphistr=str(sphi)
    tf.write(sphistr + '   ')
    estphi=musig[0]
    estphstr=str(estphi)
    tf.write(estphstr + '   ')
    #nmsmstr=str(musig[2])
    abser=np.abs(estphi-sphi)
    abserstr=str(abser)
    time=musig[2]
    timestr=str(time)
    timev.append(time)
    tf.write(abserstr + '    ')
    #tf.write(nmsmstr + '    ')
    tf.write(timestr + '\r\n')
    sharp.append(np.exp(1j*(sphi-estphi)))
tf.close()
#plt.figure(figsize=(8.1,5))
#plt.plot (, markerfacecolor='crimson', markeredgecolor='crimson', marker='o', markersize = 6)
#plt.plot (t2, p3, 'crimson',  linewidth=2)
#plt.show()
totalt=np.sum(timev)/dimtot
print('tot avg time')
print(totalt)
avgsharp=np.absolute(np.sum(sharp)/dimtot)
holv=avgsharp**(-2)-1.0
print('holevo variance')
print(holv)
sigavg=np.sum(sigv)/dimtot
mse=np.sum(mse0)/dimtot
print('variance avg')
print(sigavg)
print('average mse')
print(mse)

    #print('sigma')
    #print(musig[1])
