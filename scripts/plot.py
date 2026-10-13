#!/usr/bin/env python3
"""Plot a demo CSV. Plotting runs separately from the C++ numerical library."""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
p=argparse.ArgumentParser()
p.add_argument("csv"); p.add_argument("output")
args=p.parse_args()
data=np.genfromtxt(args.csv,delimiter=",",names=True)
x=np.unique(data["x"]); y=np.unique(data["y"])
u=data["T"].reshape(len(y),len(x))
fig,ax=plt.subplots(figsize=(6,5),layout="constrained")
im=ax.pcolormesh(x,y,u,shading="nearest",cmap="inferno")
ax.set(xlabel="x",ylabel="y",title="Interior temperature",aspect="equal",xlim=(0,1),ylim=(0,1))
fig.colorbar(im,ax=ax,label="Temperature")
fig.savefig(args.output,dpi=180)
