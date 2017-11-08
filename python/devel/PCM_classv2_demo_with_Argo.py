# -*coding: UTF-8 -*-
# Here we show how to use PCM implementation version 2
__author__ = 'gmaze@ifremer.fr'

import sys, os
os.chdir('/Users/gmaze/git/github/obidam/pcm/python/devel/')
execfile('PCM_classv2_demo_with_Argo_startup.py')

# ----------------
# In[Compute a PCM]
#%time
# Number of class in the PCM
K = 8
# K = 3

# Compression level for the dimensionality reduction
maxvar = 99.9 # in %

# Create a PCM instance with same vertical axis as Argo data:

# Usage: pcm(K, DPTmodel, scaling=1, reduction=1, classif='gmm', COVARTYPE='full', maxvar=99.9, verb=False)
m = pcm(K, Z, maxvar=maxvar)
print m
print dir(m)

# Fit the PCM on the training set:
m = m.fit(X, Z)

# Then predict classes:
labels = m.predict(X, Z)

# Or we can do both in one step:
# labels = m.fit_predict(X, Z)

# Compute posterior probability of each class:
post = m.predict_proba(X, Z)

# Compute quantiles (Typical profiles of classes)
Q = m.quant(X, labels=labels, verb=True)
Q.attrs['units'] = Xunit

# Compute quantiles on Normalized data:
xn = m._scaler.transform(m._interpoler.fit_transform(X, Z))
Qn = m.quant(xn, labels=labels, verb=True)
Qn.attrs['units'] = ("Norm[%s]")%(Xunit)

# ----------------
# Create map of labels:
fig,ax,proj = init_map([-82,0,0,65], dxtick=10, dytick=5)

# Add scatter plot
cmap = cmap_discretize(plt.cm.Paired,K)
sc = ax.scatter(lon,lat,s=8,c=labels,edgecolors='face',
            vmin=0,
            vmax=K-1,
            transform=proj,cmap=cmap)
# Finish map:
discrete_colorbar(K, cmap, ticklabels=np.arange(0,K)).set_label("PCM class")
ax.set_title("LABELS (K=%i)"%(K),fontsize=18)
plt.tight_layout()
plt.show()


# Create map of posteriors:
fig,ax,proj = init_map([-82,0,0,65],dxtick=10,dytick=5)

# Add scatter plot
k = np.random.randint(0,K) # The component to plot posteriors for, in range(K)
# cmap = plt.cm.jet
cmap = plt.cm.YlGn
sc = ax.scatter(lon,lat,s=8,c=post[:,k],edgecolors='face',
            vmin=0,
            vmax=1,
            transform=proj,cmap=cmap)
# Finish map:
plt.colorbar(sc, shrink=0.5)
ax.set_title("POSTERIORS FOR COMPONENT k=%i"%(k),fontsize=18)
plt.tight_layout()
plt.show()

# Plot quantiles
fig, ax = plt.subplots(nrows=1, ncols=m.K, figsize=(2*m.K,4), dpi=80, facecolor='w', edgecolor='k',sharey='row')
cmap = cmap_discretize(plt.cm.Paired,m.K)
xlim = np.array([0.9*Q.min(), 1.1*Q.max()])
for k in range(m.K):
    Qk = Q.sel(components=k)
    for q in Qk['quantile']:
        Qkq = Qk.sel(quantile=q)
        ax[k].plot(Qkq.values.T,Z,label=("%0.2f")%(Qkq['quantile']))
    ax[k].set_title(("Component: %i")%(k),color=cmap(k))
    ax[k].legend(loc='lower right')
    ax[k].set_xlim(xlim)
    ax[k].set_ylim(np.array([Z.min(), Z.max()]))
    ax[k].set_xlabel(Q.units)
    if k==0: ax[k].set_ylabel('feature dimension')
    ax[k].grid(True)
plt.tight_layout()


# Plot quantiles
fig, ax = plt.subplots(nrows=1, ncols=m.K, figsize=(2*m.K,4), dpi=80, facecolor='w', edgecolor='k',sharey='row')
cmap = cmap_discretize(plt.cm.Paired,m.K)
xlim = np.array([-4,4])
for k in range(m.K):
    Qk = Qn.sel(components=k)
    for q in Qk['quantile']:
        Qkq = Qk.sel(quantile=q)
        ax[k].plot(Qkq.values.T,Z,label=("%0.2f")%(Qkq['quantile']))
    ax[k].set_title(("Component: %i")%(k),color=cmap(k))
    ax[k].legend(loc='lower right')
    ax[k].set_xlim(xlim)
    ax[k].set_ylim(np.array([Z.min(), Z.max()]))
    ax[k].set_xlabel(Qn.units)
    if k==0: ax[k].set_ylabel('feature dimension')
    ax[k].grid(True)
plt.tight_layout()


fig, ax = plt.subplots(nrows=2, ncols=m.K, figsize=(2*m.K,8), dpi=80, facecolor='w', edgecolor='k',sharey='row')
cmap = cmap_discretize(plt.cm.Paired,m.K)

xlim = np.array([0.9*Q.min(), 1.1*Q.max()])
for k in range(m.K):
    Qk = Q.sel(components=k)
    for q in Qk['quantile']:
        Qkq = Qk.sel(quantile=q)
        ax[0,k].plot(Qkq.values.T,Z,label=("%0.2f")%(Qkq['quantile']))
    ax[0,k].set_title(("Component: %i")%(k),color=cmap(k))
    ax[0,k].legend(loc='lower right')
    ax[0,k].set_xlim(xlim)
    ax[0,k].set_ylim(np.array([Z.min(), Z.max()]))
    ax[0,k].set_xlabel(Q.units)
    if k==0: ax[0,k].set_ylabel('feature dimension')
    ax[0,k].grid(True)

xlim = np.array([-4,4])
for k in range(m.K):
    Qk = Qn.sel(components=k)
    for q in Qk['quantile']:
        Qkq = Qk.sel(quantile=q)
        ax[1,k].plot(Qkq.values.T,Z,label=("%0.2f")%(Qkq['quantile']))
    ax[1,k].legend(loc='lower right')
    ax[1,k].set_xlim(xlim)
    ax[1,k].set_ylim(np.array([Z.min(), Z.max()]))
    ax[1,k].set_xlabel(Qn.units)
    if k==0: ax[1,k].set_ylabel('feature dimension')
    ax[1,k].grid(True)
plt.tight_layout()

