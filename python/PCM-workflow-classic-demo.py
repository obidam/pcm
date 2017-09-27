# coding: utf-8

# # Profile Classification Model with python
# ## Classic workflow demonstration
# Created on Sat Feb  4 08:32:46 2017
# 
# We're going to classify Argo data with the following workflow:
#  - Load netcdf Argo profiles on Standard Depth Levels
#  - Normalize profiles
#  - Compress the collection with PCA (dimensionality reduction)
#  - Train a GMM on the reduced dataset
#  - Classify all the dataset
#  - Compute class profile statistics
#  - Create a map with labels by profiles
# 
# Along the way, we'll produce figures to be compared with Maze et al, POC, 2017
# 
# @author: gmaze@ifremer.fr
# 
# More documentations:
# - http://scikit-learn.org
# - http://xarray.pydata.org
# - http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
# - http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
# - http://scikit-learn.org/stable/modules/mixture.html

# ## Import modules and libraries

# In[1]:

import os
import numpy as np
import xarray, dask
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap as bm
import pandas as pd

# http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
from sklearn import preprocessing

# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
from sklearn.decomposition import PCA

# http://scikit-learn.org/stable/modules/mixture.html
from sklearn.mixture import GaussianMixture

import seaborn as sns


# ## Define key parameters of the analysis

# In[2]:

# Number of GMM classes to compute:
K = 8

# Compression level for the dimensionality reduction 
maxvar = 99.9 # in %


# In[3]:

# Define usefull functions:
def vrange(V): # Return the array value range
    return [np.min(V),np.max(V)]
def vrangec(V): # Return the array value centered range 
    xl = np.max(np.abs(vrange(V)))
    return np.array([-xl,xl])


# ## Load data
# Create the X[Np,Nz] array

# In[4]:

# Define where to find the Argo file on STD:
ncroot = '/home5/pharos/venthsal/data/ARGO/copoda_db/setup_H/db_thd_config6_last/gmm'
ncfile = 'NATL_HOMOGENEOUS_variables_7subset_1.nc'
#ncfile = 'NATL_HOMOGENEOUS_variables_7subset_*.nc'

# Load Argo SDL data:
dset = xarray.open_mfdataset(os.path.join(ncroot,ncfile))

# Select the depth layer for training:
#dsub = dset.sel(DEPTH=slice(-100,-1000))
dsub = dset.sel(DEPTH=slice(0,-1400))

# And Create the array X(Nz,Np): The field to classify with a GMM, 
# Np profiles with Nz depth levels.
X = dsub['TEMP'].values
DPTmodel = dsub['DEPTH'].values
lon = dsub['LONGITUDE'].values
lat = dsub['LATITUDE'].values
          
# Size of the training set X:
[Np, Nz] = X.shape
print "Number of raw features (Depth Levels): ", Nz 
print "Number of samples (N profiles): ", Np


# ## Normalization
# We operate along feature dimensions to:
# - Remove the sample mean
# - Divide by the sample std
# 
# one depth level at a time
# 
# Create the Xn[Np,Nz] array
# 
# Doc:
# http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html

# In[5]:

# Fit the scaler object:
scaler = preprocessing.StandardScaler()
scaler = scaler.fit(X)

# The mean and std profiles are in the scaler object properties:
X_ave = scaler.mean_
X_std = scaler.scale_

# Normalize data:
Xn = scaler.transform(X)       

# Here, we only center data:
Xc = preprocessing.StandardScaler(with_std=False).fit(X).transform(X)

# Compute additional statistics, like the observed PDFz:
def diag_pdfz(X,xedges):
    Nz = X.shape[1]
    PDFz = np.zeros((xedges.shape[0]-1,Nz))
    for iz in np.arange(Nz):
        h, hx = np.histogram(X[:,iz],bins=xedges,density=True)
        PDFz[:,iz] = h
    PDFz_axis = hx[0:-1]
    return PDFz, PDFz_axis
PDFz, PDFz_axis = diag_pdfz(X,np.arange(0,30,0.2))


# ### Time for a figure

# In[6]:

fig, ax = plt.subplots(nrows=1, ncols=2, sharey='row', figsize=(10,5), dpi=80, facecolor='w', edgecolor='k')
ax[0].plot(X_ave, DPTmodel, '-', linewidth=2,label='Sample Mean')
ax[1].plot(X_std, DPTmodel, '-', linewidth=2,label='Sample Std')
# tidy up the figure
ax[0].set_ylabel('Model Depth')
for ix in range(0,2):
    ax[ix].legend(loc='lower right')
    ax[ix].grid(True)
    ax[ix].set_xlabel('[degC]')
fig.suptitle('Training Set Standardization profiles', fontsize=12)
plt.show()


# In[7]:

# Select 100 random profiles:
n = 150
ip = np.unique(np.random.randint(0,Np-1,n))

# Random selection of profiles
fig, ax = plt.subplots(nrows=1, ncols=2, sharey='row', figsize=(10,5), dpi=80, facecolor='w', edgecolor='k')
ax[0].plot(X[ip,:].T, np.reshape(np.repeat(DPTmodel,ip.shape[0]),[Nz,ip.shape[0]]), '-', color='gray', linewidth=1)
ax[0].plot(np.mean(X[ip,:].T,axis=1), DPTmodel, '-', color='k', linewidth=2)
ax[0].plot(np.mean(X[ip,:].T,axis=1)-np.std(X[ip,:].T,axis=1), DPTmodel, '--', color='k', linewidth=2)
ax[0].plot(np.mean(X[ip,:].T,axis=1)+np.std(X[ip,:].T,axis=1), DPTmodel, '--', color='k', linewidth=2)
ax[0].grid(True)
ax[0].set_title('Temperature profiles (random selection)')
ax[0].set_xlabel('[degC]')

cmap = plt.get_cmap('Paired',lut=128)
df = xarray.DataArray(PDFz.T, coords=[DPTmodel,PDFz_axis], dims=['Depth','[degC]'])
p = df.plot(cmap=cmap,vmin=0,vmax=0.6,ax=ax[1])
ax[1].set_title("Observed PDFz")

fig.suptitle('Training Set profiles', fontsize=12)
plt.show()


# In[8]:

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
sns.set(context="notebook",style="whitegrid", palette="deep", color_codes=True)
iz = 30
sns.distplot(X[:,iz], norm_hist=True, color="m", axlabel="PDF at %0.2fm depth"%(DPTmodel[iz]))
plt.show()


# In[9]:

iz1 = 1
iz2 = np.argmax(DPTmodel<=-300)
x = X[:,iz1]
y = X[:,iz2]
#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
sns.set(context="notebook",style="whitegrid", color_codes=True)
sns.jointplot(x=x, y=y, kind="kde",              xlim=vrange([x,y]),ylim=vrange([x,y]),              size=7, ratio=4, space=0)            .set_axis_labels("PDF at %0.2fm depth"%(DPTmodel[iz1]), "PDF at %0.2fm depth"%(DPTmodel[iz2]))
plt.suptitle("%0.2fm vs %0.2fm depth"%(DPTmodel[iz1],DPTmodel[iz2]))
plt.show()


# ## Reduction
# Now that we have a standardized collection of profiles, let's compress it with a PCA decomposition
# 
# The goal here, is to reduce the dimensionality of the problem from N depths levels down to a couple of principal components
# 
# \begin{eqnarray}
# 	\mathbf{x}(z) &=& \sum_{j=1}^{Nc} \mathbf{P}(z,j) \mathbf{y}(j)
# \end{eqnarray}
# where $\mathbf{P}\in \mathbb{R}^{Nz\times Nc}$ and $\mathbf{y}\in \mathbb{R}^{Nc\times Np}$ with $Nc\leq Nz$. 
# The first rows of $\mathbf{P}$ contain profiles maximizing the structural variance throughout the collection of profiles.
# 
# Create the Xr[Np,Nc] array
# 
# Doc: http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html

# In[10]:

# Compute P (the EOFs) from x:
reducer = PCA(n_components=maxvar/100,svd_solver='full')
reducer.fit(Xn)

# Reduce the dataset (compute the y):
Xr = reducer.transform(Xn) # Here we compute: np.dot(Xn - reducer.mean_, np.transpose(reducer.components_))

# Variables of the reduced space:
Nc = reducer.n_components_ # Number of components retained
EOFs = reducer.components_ # [Nc , Nz], the P matrix
V = reducer.explained_variance_ratio_ # Explained variance, with 0 to 1 values

# We can also compute EOFs with real units this way:
S = np.sqrt(reducer.explained_variance_*Xn.shape[0]) # These are the singular values
Z = np.dot(Xn - reducer.mean_, np.transpose(reducer.components_)) # This is simply Xr or the principal components
Ztilde = Z/np.sqrt(S) # Normalized PCs
#EOFs_real = np.dot(np.transpose(Ztilde),X)/X.shape[0] # Regression on any collection of profiles
EOFs_realc = np.dot(np.transpose(Ztilde),Xc)/Xc.shape[0] # Regression on any collection of profiles
EOFs_real = np.dot(np.transpose(Ztilde),Xn)/Xn.shape[0] # Regression on any collection of profiles

# Compute the RMS difference between the reconstructed and original dataset:
Xn_reconstructed = reducer.inverse_transform(Xr)
X_reconstructed = scaler.inverse_transform(Xn_reconstructed)
rms = np.sqrt(np.mean(np.square(X_reconstructed-X),axis=0))

#
print "\nWe reduced the dimensionality of the problem from %i depth levels down to %i PCs\n"%(Nz,Nc)


# ### Figures with PCA decomposition details

# In[11]:

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
ax.bar(range(0,Nc),V*100)
ax.set_xlabel('i PCA component')
ax.set_ylabel('% of variance explained')
ax.grid(True)
ax.set_xticklabels(range(0,Nc))
ax.set_title('Variance explained')
plt.show()


# In[12]:

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20,8), dpi=160, facecolor='w', edgecolor='k', sharey='row')

iax = 0
xl = np.max(np.abs([np.min(EOFs),np.max(EOFs)]))
for ie in range(0,4):
    ax[iax].plot(np.transpose(EOFs[ie,:]),DPTmodel,label="PCA-%i"%(ie+1))
ax[iax].set_xlim(1.1*np.array([-xl,xl]))
ax[iax].legend(loc='lower right')
ax[iax].set_xlabel('PCA components (no units)')
ax[iax].set_ylabel('Model Depth')
ax[iax].grid(True)
ax[iax].set_title('PCA 1-4 (no units)')

iax+=1
xl = np.max(np.abs([np.min(EOFs_real),np.max(EOFs_real)]))
for ie in range(0,4):
    ax[iax].plot(np.transpose(EOFs_real[ie,:]),DPTmodel,label="PCA-%i"%(ie+1))
ax[iax].set_xlim(1.1*np.array([-xl,xl]))
ax[iax].legend(loc='lower left')
ax[iax].set_xlabel('PCA components (normalized units)')
ax[iax].set_ylabel('Model Depth')
ax[iax].grid(True)
ax[iax].set_title('PCA 1-4 (real units)')

iax+=1
xl = np.max(np.abs([np.min(EOFs_realc),np.max(EOFs_realc)]))
for ie in range(0,4):
    ax[iax].plot(np.transpose(EOFs_realc[ie,:]),DPTmodel,label="PCA-%i"%(ie+1))
ax[iax].set_xlim(1.1*np.array([-xl,xl]))
ax[iax].legend(loc='lower left')
ax[iax].set_xlabel('PCA components (centered units)')
ax[iax].set_ylabel('Model Depth')
ax[iax].grid(True)
ax[iax].set_title('PCA 1-4 (real units)')
plt.show()


# In[13]:

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
ax.semilogx(rms, DPTmodel, '-', linewidth=2,label="d=%i (%0.2f %%)"%(Nc,np.sum(V*100)))
# tidy up the figure
ax.legend(loc='upper left')
ax.set_xticks([0.0125,0.025,0.05,0.1,0.25,0.5,1])
ax.set_xlim([1e-3,1e1])
ax.set_xlabel('[real unit]')
ax.set_ylabel('Model Depth')
ax.grid(True)
fig.suptitle('RMS difference between reduced and original dataset', fontsize=12)
plt.show()


# ## GMM Classification
# We classify with a GMM the reduce dataset
# 
# Doc: http://scikit-learn.org/stable/modules/mixture.html

# In[14]:

# Set-up and train the classifier:
gmm = GaussianMixture(n_components=K,                      covariance_type='full',                      init_params='kmeans',                      max_iter=1000,                      tol=1e-6)
gmm.fit(Xr) # Training on reduced data

# Extract GMM parameters:
priors = gmm.weights_ # [K,1]
centers= gmm.means_   # [K,Nc]
covars = gmm.covariances_ # [K,Nc,Nc] if 'full'

# Classify the dataset:
LABELS = gmm.predict(Xr) # [Np,1]
POST   = gmm.predict_proba(Xr) # [Np,Nc]


# ## Time for a lot of figures

# In[15]:

def plot_GMMellipse0(gmm,id,ik,col,ax,label=""):
    """
        Plot an 1-STD ellipse for a given component (ik) and 2 dimensions (id) 
        of the GMM model gmm
        This is the class routine, using the matplotlib Ellipse method
        I don't like it because the Ellipse object cannot be labelled...
    """
    covariances = gmm.covariances_[ik][(id[0],id[0],id[1],id[1]),(id[0],id[1],id[0],id[1])].reshape(2,2)
    v, w = np.linalg.eigh(covariances)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan2(u[1], u[0])
    angle = 180 * angle / np.pi 
    v = 2. * np.sqrt(v)
    ell = mpl.patches.Ellipse(gmm.means_[ik,(id[0],id[1])], v[0], v[1],180 + angle,                                         fill=False,label=label)
    ell.set_clip_box(ax.bbox)
    ell.set_alpha(1)
    ell.set_label("Class #%i"%(ik))
#    ell.set_facecolor(color)
    ell.set_edgecolor(col)
    ell.set_linewidth(2)
    ax.add_artist(ell)
    return ell,ax

def plot_GMMellipse(gmm,id,ik,col,ax,label="",std=[1],main_axes=True,**kwargs):
    """
        Plot an 1-STD ellipse for a given component (ik) and 2 given dimensions (id) 
        of the GMM model gmm
        This is my routine, simply working with a matplotlib plot method
        I also added the possiblity to plot the main axes of the ellipse
    """
    covariances = gmm.covariances_[ik][(id[0],id[0],id[1],id[1]),(id[0],id[1],id[0],id[1])].reshape(2,2)
    d, v = np.linalg.eigh(covariances) #  eigenvectors have unit length
    d = np.diag(d)
    theta = np.arange(0,2*np.pi,0.02)
    x = np.sqrt(d[0,0])*np.cos(theta)
    y = np.sqrt(d[1,1])*np.sin(theta)
    xy = np.array((x,y)).T
    ii = 0
    for nstd in np.array(std):
        ii+=1
        ellipse = np.inner(v,xy).T
        ellipse = nstd*ellipse + np.ones((theta.shape[0], 1))*gmm.means_[ik,(id[0],id[1])]
        if ii == 1:
#            p = ax.plot(ellipse[:,0], ellipse[:,1], color=col, axes=ax, label=("%s (%i-std)")%(label,nstd),**kwargs)
            p = ax.plot(ellipse[:,0], ellipse[:,1], color=col, axes=ax, label=("%s")%(label),**kwargs)
        else:
            p = ax.plot(ellipse[:,0], ellipse[:,1], color=col, axes=ax,**kwargs)
    if main_axes: # Add Main axes:
        for idir in range(2):
            l = np.sqrt(d[idir,idir])*v[:,idir].T
            start = gmm.means_[ik,(id[0],id[1])]-l
            endpt = gmm.means_[ik,(id[0],id[1])]+l
            linex = [start[0], endpt[0]]
            liney = [start[1], endpt[1]]
            plt.plot(linex,liney,color=col,axes=ax)
    return p,ax


# In[16]:

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8), dpi=80, facecolor='w', edgecolor='k')
id = np.array([1,2])-1
ax.scatter(Xr[:,id[0]],Xr[:,id[1]],1)
ax.grid(True)
ax.axis('equal')
ax.set_xlabel("Dimension #%i"%(id[0]+1))
ax.set_ylabel("Dimension #%i"%(id[1]+1))
colors = mpl.colors.cnames.items()
colors = iter(plt.cm.rainbow(np.linspace(0, 1, K)))
for ik in np.arange(K):
#    el,ax = plot_GMMellipse0(gmm,id,ik,colors[ik][0],ax,label="Class-%i"%(ik+1))
    el,ax = plot_GMMellipse(gmm,id,ik,next(colors),ax,label="Class-%i"%(ik+1),linewidth=4)
#    el,ax = plot_GMMellipse(gmm,id,ik,next(colors),ax,label="Class-%i"%(ik+1),std=[1],linewidth=4)
ax.legend(loc='upper left')
plt.show()


# In[17]:

def sns_GMMellipse(x,y,gmm=[],id=[],std=[1],main_axes=True,label="?",**kwargs):
    """
        Plot an 1-STD ellipse for a given component (ik) and 2 given dimensions (id) 
        of the GMM model gmm
        This is my routine, simply working with a matplotlib plot method
        I also added the possiblity to plot the main axes of the ellipse
    """
    K = gmm.n_components
#    colors = iter(plt.cm.rainbow(np.linspace(0, 1, K)))
    colors = sns.color_palette("Paired", K)
    for ik in np.arange(K):
#        col = next(colors)
        col = colors[ik]
        covariances = gmm.covariances_[ik][(id[0],id[0],id[1],id[1]),(id[0],id[1],id[0],id[1])].reshape(2,2)
        d, v = np.linalg.eigh(covariances) #  eigenvectors have unit length
        d = np.diag(d)
        theta = np.arange(0,2*np.pi,0.02)
        x = np.sqrt(d[0,0])*np.cos(theta)
        y = np.sqrt(d[1,1])*np.sin(theta)
        xy = np.array((x,y)).T
        ii = 0
        for nstd in np.array(std):
            ii+=1
            ellipse = np.inner(v,xy).T
            ellipse = nstd*ellipse + np.ones((theta.shape[0], 1))*gmm.means_[ik,(id[0],id[1])]
            if ii == 1:
                plt.plot(ellipse[:,0], ellipse[:,1], color=col, label=("%s")%(label),**kwargs)
            else:
                plt.plot(ellipse[:,0], ellipse[:,1], color=col, **kwargs)
        if main_axes: # Add Main axes:
            for idir in range(2):
                l = np.sqrt(d[idir,idir])*v[:,idir].T
                start = gmm.means_[ik,(id[0],id[1])]-l
                endpt = gmm.means_[ik,(id[0],id[1])]+l
                linex = [start[0], endpt[0]]
                liney = [start[1], endpt[1]]
                plt.plot(linex,liney,color=col)


# In[18]:

df = pd.DataFrame(Xr[:,id], columns=["x", "y"])
sns.set(context="notebook",style="whitegrid", color_codes=True)
g = sns.JointGrid(x="x", y="y", data=df, size=7, ratio=4, space=0.1,
                  xlim=vrangec(Xr[:,id]),ylim=vrangec(Xr[:,id]))

g.plot_marginals(sns.distplot, kde=False, color=".5")

g.plot_joint(plt.scatter, c=".5", s=1, linewidth=1, marker="+")

g.plot_joint(sns_GMMellipse, gmm=gmm, id=id, main_axes=False, linewidth=3)

g.set_axis_labels("Dimension #%i"%(id[0]+1), "Dimension #%i"%(id[1]+1))
g.fig.suptitle("PC-%i vs PC-%i"%(id[1]+1,id[0]+1))
plt.show()


# In[19]:

for ii in np.arange(2,4):
    id = np.array([1,ii])-1
    df = pd.DataFrame(Xr[:,id], columns=["x", "y"])
    sns.set(context="notebook",style="whitegrid", color_codes=True)
    g = sns.JointGrid(x="x", y="y", data=df, size=7, ratio=4, space=0.1,
                      xlim=vrangec(Xr[:,id]),ylim=vrangec(Xr[:,id]))
    #g.plot_joint(sns.kdeplot, cmap="Purples_d", kind='hex', linewidth=1, color='k', n_levels=30)
    #g.plot_joint(sns.kdeplot, shade = True, cmap="Purples_d", kind='hex', n_levels=30)
    g.plot_joint(sns.kdeplot, shade = True,                  shade_lowest=False,                 cmap=sns.light_palette("gray",reverse=False,as_cmap=True), kind='hex', n_levels=30)
    #g.plot_joint(sns.kdeplot, shade = False, kind='hex', n_levels=10)
    g.plot_marginals(sns.kdeplot, color="k", shade=True)

    g.plot_joint(sns_GMMellipse, gmm=gmm, id=id, main_axes=False, linewidth=3)

    g.set_axis_labels("Dimension #%i"%(id[0]+1), "Dimension #%i"%(id[1]+1))
    g.fig.suptitle("PC-%i vs PC-%i"%(id[1]+1,id[0]+1))
    plt.show()


# # Classify all the dataset

# ## Load the entire dataset to process

# In[71]:

get_ipython().run_cell_magic(u'time', u'', u'\n# Define where to find the Argo file on STD:\nncfile = \'NATL_HOMOGENEOUS_variables_7subset_*.nc\'\n#ncfile = \'NATL_HOMOGENEOUS_variables_7subset_2.nc\'\n\n# Load Argo SDL data:\ndset = xarray.open_mfdataset(os.path.join(ncroot,ncfile))\n\n# Select the depth layer of analysis:\ndsub = dset.sel(DEPTH=slice(0,-1400))\n\n# And get the array X(Nz,Np): The field to classify with a GMM, \n# Np profiles with Nz depth levels.\nX = dsub[\'TEMP\']\n\n# Size of the training set X:\n[Np, Nz] = X.shape\nprint "Number of raw features (Depth Levels): ", Nz \nprint "Number of samples (N profiles): ", Np')


# ## Apply the PCM to the dataset:

# In[72]:

get_ipython().run_cell_magic(u'time', u'', u'# Normalize data:\nXn = scaler.transform(X)    \n\n# Reduce the dataset (compute the y):\nXr = reducer.transform(Xn) # Here we compute: np.dot(Xn - reducer.mean_, np.transpose(reducer.components_))\n\n# Classify the dataset:\nLABELS = gmm.predict(Xr) # [Np,1]\nPOST   = gmm.predict_proba(Xr) # [Np,Nc]')


# ## Time for figures

# In[73]:

# Plot LABELS map:
lon = dsub['LONGITUDE'].values
lat = dsub['LATITUDE'].values

mp0 = bm(projection='cyl',        llcrnrlat=-1,urcrnrlat=71,        llcrnrlon=-91,urcrnrlon=1,        lat_ts=0,resolution='c')
fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(5,5), dpi=160, facecolor='w', edgecolor='k') 
axes = np.reshape(axes,(1,1))
for ie in range(1):
    mp = mp0
    mp.ax = axes[0,ie]
    plt.jet()
    sc = mp.scatter(lon,lat,1,LABELS)
    mp.drawmeridians(np.arange(0, 360, 10), labels=[0,0,0,1], color='0.8', linewidth=0.5, fontsize=8)
    mp.drawparallels(np.arange(-80, 80, 5), labels=[1,0,0,0], color='0.8', linewidth=0.5, fontsize=8)
    mp.drawlsmask(land_color='0.8', ocean_color='w')
#    mp.drawcoastlines()
    mp.colorbar(sc,ax=axes[0,ie])
#    ax.set_title(tit)
    axes[0,ie].set_title("LABELS (K=%i)"%(K))
plt.tight_layout()
plt.show()


# In[74]:

# Plot POSTERIORS map:
mp0 = bm(projection='cyl',        llcrnrlat=-1,urcrnrlat=71,        llcrnrlon=-91,urcrnrlon=1,        lat_ts=0,resolution='c')
for ie in range(K):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(3,3), dpi=160, facecolor='w', edgecolor='k') 
    mp = mp0
    mp.ax = ax
    mp.drawmeridians(np.arange(0, 360, 10), labels=[0,0,0,1], color='0.8', linewidth=0.5, fontsize=5)
    mp.drawparallels(np.arange(-80, 80, 5), labels=[1,0,0,0], color='0.8', linewidth=0.5, fontsize=5)
#    mp.drawlsmask(ax=ax,land_color='0.8', ocean_color='w')
    mp.drawcoastlines(ax=ax, linestyle='solid', color='0.8', linewidth=0.5)
    sc = mp.scatter(lon,lat,1,POST[:,ie])
#    mp.colorbar(sc,ax=axes[0,ie])
    ax.set_title("POST (K=%i)"%(ie+1),fontsize=8)
    plt.tight_layout()
    plt.show()

