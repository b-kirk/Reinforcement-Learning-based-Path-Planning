import numpy as np
import math
from matplotlib import path
import plotly.express as px
import plotly.graph_objects as go
import webbrowser
import torch
firefox_path = "C:\\Program Files\\Mozilla Firefox\\firefox.exe" #define the Path to firefox
webbrowser.register('firefox', None,webbrowser.BackgroundBrowser(firefox_path))

def cholcov(s):
    s = np.atleast_2d(np.asarray(s))
    n,m = np.size(s,0), np.size(s, 1)
    if n==m:
        s = torch.Tensor(s)
        try:
            T = torch.linalg.cholesky(s, upper=True).numpy()
        except:
            U, D = np.linalg.eig((s+s.H)/2)
            maxind = U.index(max(U,))
            # Must account for non-positive-definite
    else:
        T = np.zeros(0, like=s)
        p = np.nan(like=s)
    
    return T.flatten()


def prior(s):
    q = np.size(s["x"])
    if s["x"].ndim == 1: 
        m, n = [q, 1]
    else:
        m, n = q
    tempx = np.zeros([m,])
    tempx[(s["Q"]>0) & (s["u"] > 0) & (s["ci"] > 0) & (s["cii"] > 0) & (s["z"] > 0)] = 4
    px = np.argwhere(tempx!=4)

    return px

PI = torch.pi

def normaliseWeight(theta, Wp, N):
    pri = prior(theta)

    if len(pri) >= 1:
        Wp[pri.conj().T]=0
    Wp[Wp==None]=0
    
    Wpnorm = Wp/sum(Wp)
    if sum(Wpnorm)==0:
        Wpnorm = np.ones(N,)/N
    
    return Wpnorm

def Pasquil_Gaussian_Plume(s, p):
    x = p["x_matrix"]
    y = p["y_matrix"]
    z = p["z_matrix"]
    D = s["ci"]
    t = s["cii"]

    lmbda = np.sqrt(np.asarray(((D*t)/(1+(s["u"]**2)*t/(4*D)))))

    module_dist = np.sqrt(((s['x'] - x))**2+((s["y"]-y))**2+((s["z"]-z))**2)
    angazimuth = np.arctan2((x-s["x"]),(y-s["y"]))
    angazimuth = np.asarray(angazimuth)
    angazimuth[angazimuth < 0] = angazimuth[angazimuth < 0]+2*PI
    ang_diff_azimuth = angazimuth+s["phi"]-PI/2
    ang_diff_azimuth = np.asarray(ang_diff_azimuth)
    ang_diff_azimuth[ang_diff_azimuth < 0] = ang_diff_azimuth[ang_diff_azimuth < 0]+2*PI
    xr = module_dist * np.cos(ang_diff_azimuth)

    # ----- Chemical Sources

    # C = s["Q"]/(4*np.pi*D*module_dist)*np.exp(1*(xr)*s["u"]/(2*D)+(-1*module_dist/lmbda)) # This may have /0 issues?
    C = s["Q"]/(4*PI*D*module_dist)*np.exp(-module_dist/lmbda)*np.exp(-(x-s["x"])*s["u"]*np.cos(s["phi"])/2*D)*np.exp(-(y-s["y"])*s["u"]*np.sin(s["phi"])/2*D)

    # ----- Radioactive Sources

    # C = 1/(4*PI*module_dist**2)

    return C

def Likelihood_Like_Yee(C, D, Wpnorm, thresh):
    PH0 = 0.3
    PH1 = 1-PH0
    ProbBackground = 1
    
    NDsigma = 1e-4 + C
    
    sigma = 1e-4

    if D<=thresh:
        Wp = Wpnorm * ((PH0 * ProbBackground) + (PH1*(1/2)*(1+torch.erf(torch.as_tensor((thresh-C)/(NDsigma*np.sqrt(2)))).numpy())))
        Wp[Wp==None]=0
    else:
        Wp = Wpnorm * 1/sigma*np.sqrt(2*np.pi)*np.exp((-(C-D)**2)/(sigma**2))
        Wp[Wp==None]=0

    return Wp

def resampleSystematic(w, N=None):
    """ Systematic resampling method for particle filtering
        Author: Tiancheng Li, 
        Ref: T. Li, M. Bolic, P. Djuric
        Resampling methods for particle filtering
        submit to IEEE Signal Processing Magazine, August 2013

         Input:
            w    the input weight sequence 
            N    the desired length of the output sequence(i.e. the desired number of resampled particles)
        Output:
            indx the resampled index according to the weight sequence
    """
    if N == None:
        N = len(w)
    M = len(w)
    w = w/sum(w)
    Q = np.cumsum(w)
    indx = np.zeros((N,))
    T = (np.linspace(0, 1-1/N,N)+np.random.rand(1, N)/N).flatten()
    
    i = j = 0
    while i < N and j < M:
        while Q[j] < T[i]:
            j += 1
        indx[i] = j
        i += 1

    return indx.astype(int)

def resampleParticles(theta, Wpnorm, N):
    indx = resampleSystematic(Wpnorm)
    
    theta["x"] = theta["x"][indx]
    theta["y"] = theta["y"][indx]

    theta["Q"] = theta["Q"][indx]

    Covx = np.cov(theta["x"])
    Covy = np.cov(theta["y"])
    
    CovQ = np.cov(theta["Q"])
    
    dk = {
        "x": cholcov(Covx),
        "y": cholcov(Covy),
        "Q": cholcov(CovQ),
    }

    mm = 3 
    A=(4/(mm+2))**(1/(mm+4))
    cx = 4*np.pi/3
    hopt = A*(N**(-1/(mm+4)))

    return theta, dk, hopt

def mcmcResampleStep_Memory(theta,Wpnorm,N, D_k_store,thresh, dk, hopt,pos,P_k_store,PF_Memory):
    keep = []
    for zz in range(1):
        keep = []
        n = {
            "x": theta["x"] + (hopt*dk["x"]*np.random.normal(size=(N,))),
            "y": theta["y"] + (hopt*dk["y"]*np.random.normal(size=(N,))),
            "z": theta["z"],
            "Q": theta["Q"] + (hopt*dk["Q"]*np.random.normal(size=(N,))),
            "u": theta["u"],
            "phi": theta["phi"],
            "ci": theta["ci"],
            "cii": theta["cii"]
        }
        for loop in range(2):
            pri = prior(n).flatten()
            numPri = len(pri)
            for index, item in enumerate(pri):
                n["x"][item]  = theta["x"][item] + (hopt*dk["x"]*np.random.normal())
                n["y"][item] = theta["y"][item] + (hopt*dk["y"]*np.random.normal())
                n["z"][item] = theta["z"][item] 
                n["Q"][item] = theta["Q"][item] + (hopt*dk["Q"]*np.random.normal())
                n["u"][item] = theta["u"][item] 
                n["phi"][item] = theta["phi"][item] 
                n["ci"][item] = theta["ci"][item]
                n["cii"][item] = theta["cii"][item]  
        pri = prior(n).flatten()
        n["x"][pri] = theta["x"][pri]
        n["y"][pri] = theta["y"][pri]
        n["z"][pri] = theta["z"][pri] 
        n["Q"][pri] = theta["Q"][pri]
        n["u"][pri] = theta["u"][pri]
        n["phi"][pri] = theta["phi"][pri]
        n["ci"][pri] = theta["ci"][pri]
        n["cii"][pri] = theta["cii"][pri]

        nWpnorm = Wpnorm
        r = np.atleast_2d(D_k_store).reshape(-1,1)
        r = np.asarray([np.size(r,0), np.size(r, 1)])
        npP_k_store = np.asarray(P_k_store)
        if PF_Memory==1 or len(npP_k_store[0,:]) == 1:
            pos["X"] = P_k_store[-1, 0]
            pos["Y"] = P_k_store[-1, 1]
            D = D_k_store[-1]

            nC = Pasquil_Gaussian_Plume(n,pos)

            nWp = Likelihood_Like_Yee(nC, D, Wpnorm, thresh)

            nWpnorm = normaliseWeight(n,nWp,N)
        else:
            for mem in range(PF_Memory):
                ind = int(np.ceil(np.random.rand()*(r[0]-1)))
                pos["X"] = npP_k_store[ind, 0]
                pos["Y"] = npP_k_store[ind, 1]
                D = D_k_store[ind]

                nC = Pasquil_Gaussian_Plume(n, pos)

                nWp = Likelihood_Like_Yee(nC, D, nWpnorm, thresh)

                nWpnorm = normaliseWeight(n,nWp,N)

    alpha = nWpnorm/(Wpnorm) # (Wpnorm(indx))
    mcrand = np.random.rand(N,1)
    keep = np.argwhere(alpha > mcrand)
    notkeep = np.argwhere((alpha < np.random.rand()))

    theta["x"][keep] = n["x"][keep]
    theta["y"][keep] = n["y"][keep]
    theta["z"][keep] = n["z"][keep]
    theta["Q"][keep] = n["Q"][keep]
    theta["u"][keep] = n["u"][keep]
    theta["phi"][keep] = n["phi"][keep]
    theta["ci"][keep] = n["ci"][keep]
    theta["cii"][keep] = n["ci"][keep]

    Wpnorm = np.ones(N,)/N

    return theta

def UpdatePFPlume( D_k_store, theta, Wpnorm, pos, P_k_Store, thresh, N, PF_Memory, domain):
    # Estimated from particle filter

    # C = simpleGaussianPlume(theta,thresh,pos)
    C = Pasquil_Gaussian_Plume(theta, pos)
    # C[C<thresh]=0

    Wp = Likelihood_Like_Yee(C, D_k_store[-1], Wpnorm, thresh)
    pointstack = np.swapaxes(np.array([theta["x"],theta["y"]]), 0, 1)
    keep = path.Path([(domain[0],domain[2]), (domain[0], domain[3]), (domain[1], domain[3]), (domain[1], domain[2])]).contains_points(pointstack)
    
    Wp[~keep]=0
    Wpnorm = normaliseWeight(theta, Wp, N)
    Neff = 1/sum(Wpnorm**2)

    if Neff < 0.95*N : # For low release rate of sensing area
        indx = resampleSystematic(Wpnorm)
        theta, dk, hopt = resampleParticles(theta, Wpnorm, N)
        theta = mcmcResampleStep_Memory(theta, Wpnorm, N, D_k_store, thresh, dk, hopt, pos, P_k_Store, PF_Memory)
        Wpnorm = np.ones(N,)/N
    
    return theta, Wpnorm

def resampleStratified(w, N=None):
    """ Stratified resampling method for particle filtering. Author: T. Li, Ref:
        T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
        submit to IEEE Signal Processing Magazine, August 2013

        Input:
                w    the input weight sequence 
                N    the desired length of the output sequence(i.e. the desired number of resampled particles)
        Output:
                indx the resampled index according to the weight sequence
    """
    if N == None:
        N = len(w)
    M = len(w)
    w = w/sum(w)
    Q = np.cumsum(w)
    indx = np.zeros((N,))
    T = (np.linspace(0, 1-1/N, N) + np.random.rand(1, N)/N).flatten()

    i = j = 0
    while i < N and j < M:
        while Q[j] < T[i]:
            j += 1
        indx[i] = j
        i += 1

    return indx.astype(int)



UAVVel=2
sampleTime=10
bLim=900
# Simulated source parameters
# true source


s = { "Q": 5, # Release rate per t
"x": 25, # source coodinates
"y": 37.5,
"z": 1,
"u": 4, # wind speed
"phi": -270 * np.pi/180,
"ci": 1, # 0.14;  % Also s.D for Pasquil model
"cii": 8, # 0.53; % Also s.t or tau for Pasquil Model
    # ci = 0.14;
    # cii = 0.53;
    # duration = 0;
}

thresh = 5e-4

# Create rectangular domain area
xmin = 0
xmax = 50
ymin = 0
ymax = 50
zmin = 0
zmax = 50
domain = np.array([xmin, xmax, ymin, ymax]); # Size of search area

# example data
stepsize = 1; # horisontal (x and y) spacing

x_coord = np.arange(xmin,xmax+stepsize,stepsize)
y_coord = np.arange(ymin,ymax+stepsize,stepsize)
z_coord = np.arange(zmin,zmax+stepsize,stepsize)

# Create 3D Grid
X, Y, Z = np.meshgrid(x_coord,y_coord,z_coord,) # Issues with torch.meshgrid means it cannot be used here, no matter indexing

ex = {
    "x_matrix": X,
    "y_matrix": Y,
    "z_matrix": Z
}

StartingPosition = [2.0,2.0,0.0] # Starting position [x,y,z]
moveDist = 2 # How far to move

x, y, z = StartingPosition # Current position
P_k_store = []
P_k_store.append(StartingPosition)

pos = {
    "x_matrix": x,
    "y_matrix": y,
    "z_matrix": z ,
}

# Plot example dispersion from true source
# conc = simpleGaussianPlume(s,m,ex);
conc = Pasquil_Gaussian_Plume(s,ex)
conc[conc <= thresh] = float('nan')

height = 1

concSurf = conc[:,:,height]/s["Q"]
concSurf[concSurf <= thresh] = float('nan')

# Initialize PF
N = 10000 # 10000
PF_Memory = 10
resample = 0


# Uniform prior for location

a = np.ones((N,))*2
b = np.ones((N,))*5

theta = {
    "x": xmin + (xmax-xmin)*np.random.rand(N,),
    "y": ymin + (ymax-ymin)*np.random.rand(N,),
    "z": np.ones([N,])*s["z"],
    "Q": np.random.gamma(a,b), # 200*np.rand(N,1)
    "u": s["u"]*np.ones((N,)), # 2+6*np.random.rand(N,1)    0.75+0.5*np.random.rand(N,1)    0 + randn(N,1)*0.5 ?
    "phi": s["phi"]*np.ones((N,)), # (10 + 30*np.random.rand(N,1)).*np.pi/180;
    "ci": s["ci"]*np.ones((N,)),
    "cii":s["cii"]*np.ones((N,))
}

# Wp refers to particle weights
Wp = np.ones((N,))
Wpnorm = Wp/sum(Wp)
Wp = Wpnorm

timestamp = []
timestamp.append(0)
D=[]
RMSE_hist = []

sampleHistory = [StartingPosition]

for i in range(100):
    # Sim data
    # Dsim = simpleGaussianPlume(s,m,pos)
    Dsim = Pasquil_Gaussian_Plume(s,pos)

    # Dsim[Dsim < thresh] = 0
    ersize = np.size(Dsim)
    if ersize == 1:
        ersize = [1, 1]
    error = 0.6*Dsim*np.random.normal(ersize[0],ersize[1])
    Dsim = np.asarray(Dsim + error)
    Dsim[Dsim < thresh] = 0


    if np.random.rand(1,1) < 0.3:
        Dsim = np.asarray(0)
    
    D.append(Dsim)
    thetaPrev=theta
    theta, Wpnorm = UpdatePFPlume(D, theta, Wpnorm, pos, P_k_store, thresh, N, PF_Memory, domain)

    RMSE_hist.append(np.sqrt(np.mean(np.linalg.norm(np.asarray([theta["x"],theta["y"]]).conj().T-np.asarray([s["x"],s["y"]]).conj().T, axis=1).conj().T**2)))

    if i % 10 == 0:
        fig = go.Figure(data=go.Volume(
        x=ex["x_matrix"].flatten(),
        y=ex["y_matrix"].flatten(),
        z=ex["z_matrix"].flatten(),
        value=conc.flatten(),
        isomin=thresh,
        isomax=0.1,
        opacity=0.5, # needs to be small to see through all surfaces
        surface_count=100, # needs to be a large number for good volume rendering
        ))
        fig.add_scatter3d(x=theta["x"], y=theta["y"], z=theta["z"],opacity=1, marker=dict(color='green',size=2), mode='markers')
        fig.add_scatter3d(x=[s["x"]], y=[s["y"]], z=[s["z"]],opacity=1, marker=dict(color='black',size=5))
        fig.add_scatter3d(x=np.asarray(sampleHistory)[:,0], y=np.asarray(sampleHistory)[:,1], z=np.asarray(sampleHistory)[:,2],opacity=1,marker=dict(color='red',size=5, symbol='x'))
        fig.show(renderer='firefox')


    indx = resampleStratified(Wpnorm)
    t = {
        "x": theta["x"][indx],
        "y": theta["y"][indx]
    }
    
    Xneighbour = np.zeros((8,))
    Yneighbour = np.zeros((8,))
    Zneighbour = np.zeros((8,))

    # All proposed sensor locations are right, up, left, for 1 step, 2 step, & 3 step
    ynew = [moveDist,moveDist,0,-moveDist,-moveDist,-moveDist,0,moveDist]
    xnew = [0,moveDist,moveDist,moveDist,0,-moveDist,-moveDist,-moveDist]
    znew = [0,0,0,0,0,0,0,0,0,0,0,0]

    M = 40 
    MM = 1

    # Entrotaxis reward
    indx = resampleStratified(Wpnorm, M)
    d = {
        "x": theta["x"][indx],
        "y": theta["y"][indx],
        "z": theta["z"][indx],
        "Q": theta["Q"][indx],
        "u": theta["u"][indx],
        "phi": theta["phi"][indx],
        "ci": theta["ci"][indx],
        "cii": theta["cii"][indx],
    }

    # ---------------------- Entropy Reduction prediction at potential locations
    var = []
    theta_RMSE = []
    dist_theta = []
    for k in range(0, 7):
        Xneighbour[k] = pos["x_matrix"]+xnew[k]
        Yneighbour[k] = pos["y_matrix"]+ynew[k]
        Zneighbour[k] = pos["z_matrix"]+znew[k]

        if (Xneighbour[k] < xmin) or (Xneighbour[k] > xmax) or (Yneighbour[k] < ymin) or (Yneighbour[k] > ymax) or (Zneighbour[k] < zmin) or (Zneighbour[k] > zmax):
            var.append(0) # None
            dist_theta.append(0)
            theta_RMSE.append(0)
            continue
        
        npos = {
            "x_matrix": Xneighbour[k],
            "y_matrix": Yneighbour[k],
            "z_matrix": Zneighbour[k]
        }

        pC = Pasquil_Gaussian_Plume(theta, npos)

        entropy = 0
        theta_RMSE.append(0) 
        dist_theta.append(0)

        desC = Pasquil_Gaussian_Plume(d, npos)
        designd = desC + (1*desC*np.random.normal(M,MM))
        designd[np.random.rand(40,)<0.3]=0
        designd[designd<thresh]=0
        for jj in range(M):
            for jjj in range(MM):
                if designd.ndim == 1:
                    dC = designd[jj]
                else:
                    dC = designd[jj,jjj]

                zWp = Likelihood_Like_Yee(pC, dC, Wpnorm, thresh)
                zWpnorm = zWp/sum(zWp)


                theta_mean_xy = N*np.asarray([np.mean(theta["x"]*zWpnorm),np.mean(theta["y"]*zWpnorm)])
                theta_RMSE[k] = theta_RMSE[k] + np.sum(zWpnorm*np.linalg.norm(np.asarray([theta["x"],theta["y"]]).conj().T-theta_mean_xy.conj().T, axis=1).conj().T**2)/(M*MM) # Dimension issue
                dist_theta[k] = dist_theta[k] + np.linalg.norm(theta_mean_xy-[Xneighbour[k],Yneighbour[k]])/(M*MM)

        var.append(dist_theta[k]+theta_RMSE[k])
    
    val, ind = np.min(var), np.argmin(var)
    val2, ind2 = np.min(theta_RMSE), np.argmin(theta_RMSE)
    val3, ind3 = np.min(dist_theta), np.argmin(dist_theta)
    # dualControlJ

    pos["x_matrix"] = Xneighbour[ind]
    pos["y_matrix"] = Yneighbour[ind]
    pos["z_matrix"] = Zneighbour[ind]

    P_k = [pos["x_matrix"], pos["y_matrix"], pos["z_matrix"]]
    sampleHistory.append(P_k)

    move_time = np.floor(np.linalg.norm(np.asarray(P_k)-np.asarray(P_k_store[i]))/UAVVel)+sampleTime
    bLim=bLim-move_time
    if bLim <= 0:
        break
    timestamp.append(timestamp[i]+move_time)
    P_k_store.append(P_k)

    Covar = np.cov(theta["x"],theta["y"])
    Spread = np.sqrt(Covar[0,0]+Covar[1,1])

    if Spread<0:
        break
    
indx = resampleStratified(Wpnorm)

theta["x"] = theta["x"][indx]
theta["y"] = theta["y"][indx]
theta["z"] = theta["z"][indx]
theta["Q"] = theta["Q"][indx]
theta["u"] = theta["u"][indx]
theta["phi"] = theta["z"][indx]
theta["ci"] = theta["ci"][indx]
theta["cii"] = theta["cii"][indx]
    
