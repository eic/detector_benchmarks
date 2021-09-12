## plot test_matscan.cxx results
## Shujie Li, 08, 2021
## usage: python matscan_plot.py zrange rrange
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math as m
import sys
from matplotlib import colors
from scipy import stats
from matplotlib.colors import LogNorm


plt.rcParams['figure.figsize'] = [10.0, 6.0]
# plt.rcParams['fure.dpi'] = 80
            
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title]


def cart2sph(x,y,z):
    x=np.array(x)
    y=np.array(y)
    z=np.array(z)
    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    elev = np.arctan2(np.sqrt(XsqPlusYsq),z)     # theta (polar to z axis)
    az = np.arctan2(y,x)                           # phi (azimuthal)
    return r, elev, az

## read a single scan table for cylinder volume
def read_table_cylind(zmax=140,rmax=60,indir='./'):
    header  = np.array(['x', 'y', 'z', 'ind', 'material', 'Z', 'density', 'rad_length', 'thickness', 'path_length', 'sum_x0',"end_x","end_y","end_z"])
    df      = pd.read_csv(indir+"all_z%g_r%g.dat" %(zmax,rmax),names=header,delim_whitespace=True)
    df["x0"]=np.array(df["thickness"])/np.array(df["rad_length"])
    df=df.reset_index(drop=True).drop_duplicates()
    
    # calcualte polar angle then eta
    ax = df["x"]
    ay = df["y"]
    az = df["z"]
    r,polar,azi=cart2sph(ax,ay,az)
    eta = -np.log(np.tan(polar/2.0))
    df["eta"] = eta
    df["azi_angle"] = azi
    df["pol_angle"] = polar
    
    print(df.head())
    print(df.material.unique())
    
#     df=df[np.isfinite(df).all(1)] #drop in
    return df

def add_kin(df):
        # calcualte polar angle then eta
    ax = df["x"]
    ay = df["y"]
    az = df["z"]
    r,polar,azi=cart2sph(ax,ay,az)
    eta = -np.log(np.tan(polar/2.0))
    df["eta"] = eta
    df["azi_angle"] = azi
    df["pol_angle"] = polar
    return df

## group by xyz coordinate to get sum x0 per track with selected materials
def get_x0(df0,mat_name=""):
    if len(mat_name)>0:
        df=df0[df0["material"]==mat_name]
    else:
        df=df0
    # df=df[df["material"]!="CarbonFiber_25percent"]
    dfg=df.groupby(["x","y","z"])["x0"].sum()
    dfg=dfg.reset_index()
    dfg=dfg[np.isfinite(dfg).all(1)] #drop inf
    dfg=add_kin(dfg)
    return dfg


## plot mean,min,max X0 of each eta bin
def plot_1d(dfgz,mat_name="",ax="eta"):
    xlow = 0
    values=np.array(dfgz["x0"])*100
    xx = np.array(dfgz[ax])
    if ax=="eta":
        xhi=3.5
        xlow = -xhi
    elif ax=="pol_angle":
        xhi=90
        xx=180/2-xx*180/np.pi
    elif ax=="z": # z projected at barrel radius
        xhi=140
        xx*=(43.23/60)
    ## with bin stats
    nbin = 50
    x0mean,binedge,binnumber=stats.binned_statistic(xx,values, 'mean', bins=nbin,range=[xlow,xhi])
    ## !!! possible bugs in this error bar calculation
    x0max ,binedge,binnumber=stats.binned_statistic(xx,values, 'max',  bins=nbin,range=[xlow,xhi])
    x0min ,binedge,binnumber=stats.binned_statistic(xx,values, 'min',  bins=nbin,range=[xlow,xhi])

    bin_center = (binedge[0:-1]+binedge[1:])/2
    plt.plot(bin_center,x0mean,label=mat_name)
    plt.fill_between(bin_center,x0min,x0max,alpha=0.2)
    plt.xlim(xlow,xhi)
    plt.grid()
    plt.suptitle("total X0")
    plt.xlabel(ax)

    return bin_center, x0mean, x0max, x0min



if __name__ == "__main__":

    ## corresponding to the scan output filenmae from test_matscan.cxx, in cm
    zrange = int(sys.argv[1])
    rrange = int(sys.argv[2])
    outdir = './'
    indir  = './'
    cols = np.array(["x","y","z","material","thickness","path_length"])
    df   = read_table_cylind(zrange,rrange,indir)


    ## -----------------plot side view of geometry, z v.s. R-------------
    plt.figure()
    xe = df["end_x"]
    ye = df["end_y"]
    ze = df["end_z"]
    re = np.sqrt(xe**2+ye**2)
    # c5=(df["end_x"]**2+df["end_y"]**2)<11**2
    # c6=df["material"]=="Aluminum"
    plt.scatter(ze,re,s=0.001,marker=".")

    # plt.xlabel("$\eta$")
    plt.xlabel("z [cm]")
    plt.ylabel("R [cm]")
    # plt.xlim(0,4)
    plt.grid()
    plt.savefig(outdir+"/matscan_geo_z_z%g_r%g.png" %(zrange,rrange))

    ## -----------------plot side view of geometry, eta v.s. R-------------
    plt.figure()
    xe = df["end_x"]
    ye = df["end_y"]
    ze = df["end_z"]
    re = np.sqrt(xe**2+ye**2)
    plt.scatter(df["eta"],re,s=0.001,marker=".")

    plt.xlabel("$\eta$")
    plt.ylabel("R [cm]")
    # plt.xlim(0,4)
    plt.grid()
    plt.savefig(outdir+"/matscan_geo_eta_z%g_r%g.png" %(zrange,rrange))


    ## -----------------plot 1 d matscan z v.s. X0-------------
    df_all=get_x0(df)
    plt.figure()
    ax="eta"
    print(df.columns,df_all.columns)
    print(df_all['eta'])
    _=plot_1d(df_all,"Total",ax)

    mat_list=np.delete(df.material.unique(),0)
    # mat_list.drop("Vacuum")
    # mat_list=['NONE','Air','Beryllium','Aluminum','Silicon','CarbonFiber']
    # mat_list=["Aluminum","Silicon","CarbonFiber"]
    for mat_name in mat_list:
        df_mat = get_x0(df,mat_name)
        _=plot_1d(df_mat,mat_name,ax)

    plt.legend(loc="upper right")
    plt.xlabel("$\eta$")
    plt.ylabel("X/X0 (%)")
    plt.savefig(outdir+"/matscan_1d_z%g_r%g.png" %(zrange,rrange))

    ## -----------------plot 2 d matscan z v.s. phi-------------

    dfgz=df_all
    values=dfgz["x0"]
    xx = dfgz["eta"]
    yy = dfgz["azi_angle"] # in rad
    ## with bin stats
    nbin = 50
    x0mean_2d,x_edge,y_edge,binnumber_2d=stats.binned_statistic_2d(xx,yy, values, 'mean', bins=[nbin,nbin],range=[[-5,5],[-5,5]])
    # x0max_2d ,_,_,_=stats.binned_statistic_2d(xx,yy,values, 'max',  bins=[nbin,nbin],range=[[-5,5],[-5,5]])
    # x0min_2d ,_,_,_=stats.binned_statistic_2d(xx,yy,values, 'min',  bins=[nbin,nbin],range=[[-5,5],[-5,5]])

    x_center = (x_edge[0:-1]+x_edge[1:])/2
    y_center = (y_edge[0:-1]+y_edge[1:])/2 
    # plt.contour(x0mean_2d.transpose(),extent=[x_edge[0],x_edge[-1],y_edge[0],y_edge[-1]],
    #     linewidths=3, cmap = plt.cm.rainbow)
    # plt.colorbar()


    fig=plt.figure(figsize=(7, 6))
    c  = plt.pcolormesh(x_edge,y_edge/np.pi*180,x0mean_2d.transpose(),norm=LogNorm(vmin=0.001, vmax=10))
    fig.colorbar(c)
    # plt.colorbar()
    # plt.ylim(-np.pi,np.pi)
    plt.ylim(-180,180)
    plt.xlim(-5,5)
    plt.xlabel("$\eta$")
    plt.ylabel("Azimuthal angle [degree]")
    plt.suptitle("total X0 [%]")

    plt.savefig(outdir+"/matscan_2d_z%g_r%g.png" %(zrange,rrange))


