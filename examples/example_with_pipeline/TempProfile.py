import pandas as pd
import numpy as np
import scipy.constants
import matplotlib as mp
import matplotlib.pyplot as plt
import os as os

path="/scratch/achau2/Uranus/coll2/Analysis/"
filename = '2ME_b07_cold'
pathinput = '/scratch/achau2/Uranus/coll2/Analysis/'

def readfile(folder,step,ncores):
    frames=[pd.read_csv(path+folder+"/results."+'%05d'%step+"_"+'%05d'%ncores+"_"'%05d'%f+".dat",skiprows=2,header=None,sep='\t',on_bad_lines='skip') for f in range(ncores)]
    df=pd.concat(frames,ignore_index=True).sort_values(by=[0],ignore_index=True)
    df.rename(columns={0: 'id', 1: 'tag', 2: 'm',3:'x',4:'y',5 :'z', 6: 'vx', 7:'vy', 8:'vz', 9:'dens', 10:'int ener', 11:'pres', 12: 'pot ener', 13: 'entr',14:'temp'}, inplace=True)
    df['m']=pd.to_numeric(df['m'])
    df['r']=(df['x']**2+df['y']**2+df['z']**2)**0.5
    df['v']=(df['vx']**2+df['vy']**2+df['vz']**2)**0.5
    return df

def profileplot(lsdf,step,name,title,lscolor,lslabel):        
    fig, axs = plt.subplots(2, 2,figsize=(5.2,5.2),tight_layout=True)
    for i in range(2):
        for j in range(2):
            axs[i, j].grid(True)
            axs[i, j].set_xlabel('r [km]')
    for k in range(len(lsdf)):
        axs[0,0].plot(lsdf[k]['r']/1000,lsdf[k]['dens']/1000,'.',color=lscolor[k])
        axs[0,1].plot(lsdf[k]['r']/1000,lsdf[k]['temp'],'.',color=lscolor[k])
        axs[1,0].plot(lsdf[k]['r']/1000,lsdf[k]['pres']/1000000000,'.',color=lscolor[k])
        axs[1,1].plot(lsdf[k]['r']/1000,lsdf[k]['v'],'.',color=lscolor[k],label=lslabel[k])
    axs[0,0].set_ylabel('dens [g/cm³]')
    axs[0,1].set_ylabel('temp [K]')
    axs[1,0].set_ylabel('pressure [GPa]')
    axs[1,1].set_ylabel('vel [m/s]')
    fig.suptitle(title, fontsize=16)
    plt.figlegend(bbox_to_anchor=(0,0), loc="upper left", bbox_transform=fig.transFigure,ncol=4)
    plt.savefig(path+name +'_%05d'%step+'_profile.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
    plt.close
    
def plotxy(df,step,gi,qtty,b):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(1,1)
    #fig.set_size_inches(5.2,5.2*11/13)
    fig.set_size_inches(5.2,5.2*13/12)
    t=step*100
    dfs=df.sort_values(by=qtty)
    #plot=ax.scatter(x=dfs['x']/1000,y=dfs['y']/1000,c=dfs[qtty],cmap='plasma',s=0.1)
    plot=ax.scatter(x=dfs['x']/1000,y=dfs['y']/1000,c=dfs[qtty],cmap='Dark2',vmin=0,vmax=8,s=1, marker = '.')
    ax.set_ylabel('y [km]')
    ax.set_xlabel('x [km]')
    #ax.set_xlim(-1e5,1e5)
    #ax.set_ylim(-1e5,1e5)
    '''
    ax.annotate('t= '+str(step*100)+' [s]',
            xy=(-1.1e4, -1.8e4), xycoords='data',
            xytext=(20, 0), textcoords='offset points',
            horizontalalignment='center', verticalalignment='bottom')
    ax.annotate('b='+str(b)+', v=v_esc',
            xy=(2.05e4, 1.7e4), xycoords='data',
            xytext=(-20, 0), textcoords='offset points',
            horizontalalignment='center', verticalalignment='bottom')
    '''
    #cbar=plt.colorbar(plot)
    #cbar.set_label('temperature [K]')
    plt.savefig(path+gi+'_'+qtty+'_%05d'%step+'_xy_all.png',dpi=300,bbox_inches='tight', transparent=False)
    plt.close()

def plotxz(df,step,gi,qtty,b):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(1,1)
    #fig.set_size_inches(5.2,5.2*11/13)
    fig.set_size_inches(5.2,5.2*13/12)
    t=step*100
    dfs=df.sort_values(by=qtty)
    #plot=ax.scatter(x=dfs['x']/1000,y=dfs['y']/1000,c=dfs[qtty],cmap='plasma',s=0.1)
    plot=ax.scatter(x=dfs['x']/1000,y=dfs['z']/1000,c=dfs[qtty],cmap='Dark2',vmin=0,vmax=8,s=1, marker = '.')
    ax.set_ylabel('z [km]')
    ax.set_xlabel('x [km]')
    #ax.set_xlim(-1e5,1e5)
    #ax.set_ylim(-1e5,1e5)
    '''
    ax.annotate('t= '+str(step*100)+' [s]',
            xy=(-1.1e4, -1.8e4), xycoords='data',
            xytext=(20, 0), textcoords='offset points',
            horizontalalignment='center', verticalalignment='bottom')
    ax.annotate('b='+str(b)+', v=v_esc',
            xy=(2.05e4, 1.7e4), xycoords='data',
            xytext=(-20, 0), textcoords='offset points',
            horizontalalignment='center', verticalalignment='bottom')
    '''
    #cbar=plt.colorbar(plot)
    #cbar.set_label('temperature [K]')
    plt.savefig(path+gi+'_'+qtty+'_%05d'%step+'_xz_all.png',dpi=300,bbox_inches='tight', transparent=False)
    plt.close()
    
def plotxycom(df,step,gi,qtty):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(4.2,4.2*11/13)
    t=step*100
    dfs=df.sort_values(by=qtty)
    plot=ax.scatter(x=dfs['x']/1000,y=dfs['y']/1000,c=dfs[qtty],cmap='plasma',marker='o',s=0.1,zorder=1)
    ax.scatter(x=[com(df)[0]/1000],y=[com(df)[1]/1000],color='darkgreen',marker='x',s=10,linewidths=1,zorder=3)
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_xlim(-0.8e4,0.8e4)
    ax.set_ylim(-0.8e4,0.8e4)
    ax.annotate('t= '+str(step*100)+' [s]',
            xy=(-0.6e4, -0.75e4), xycoords='data',
            xytext=(20, 0), textcoords='offset points',
            horizontalalignment='center', verticalalignment='bottom')
    cbar=plt.colorbar(plot)
    cbar.set_label('pressure [GPa]')
    plt.grid(linestyle='dotted',color='#a6bddb',zorder=2)
    plt.savefig(path+gi+'_'+qtty+'_%05d'%step+'_centered_movie.png',dpi=300,bbox_inches='tight',transparent=False)
    plt.close()
    
def com(df):
    mtot=df['m'].sum()
    x=(df['x']*df['m']).sum()/mtot
    y=(df['y']*df['m']).sum()/mtot
    z=(df['z']*df['m']).sum()/mtot
    return np.array([x, y, z])

def spec_mom(df):
    return np.array([df['vx'].sum(),df['vy'].sum(),df['vz'].sum()])

def spec_ang_mom(df):
    ag_x = ((df['y'] * df['vz']) - (df['z'] *df['vy'])).sum()
    ag_y = ((df['z'] * df['vx']) - (df['x'] *df['vz'])).sum()
    ag_z = ((df['x'] * df['vy']) - (df['y'] *df['vx'])).sum()
    return np.array([ag_x,ag_y,ag_z])

def findcom(df):
    ir = df.loc[df['tag']==1]
    mcom=ir['m'].sum()
    x=(ir['x']*ir['m']).sum()
    y=(ir['y']*ir['m']).sum()
    z=(ir['z']*ir['m']).sum()
    return np.array([x,y,z])/mcom

def define_parts(df,Rp):
    #recenter the system
    com=findcom(df)
    df['x']=df['x']-com[0]
    df['y']=df['y']-com[1]
    df['z']=df['z']-com[2]
    df['r']=(df['x']**2+df['y']**2+df['z']**2)**0.5
    #compute the specific angular momenta of the particles and eccentricity:
    df['jx']=(df['y'] * df['vz']) - (df['z'] * df['vy'])
    df['jy']=(df['z'] * df['vx']) - (df['x'] * df['vz'])
    df['jz']=(df['x'] * df['vy']) - (df['y'] * df['vx'])
    df['j']=(df['jx']**2+df['jy']**2+df['jz']**2)**0.5
    #give an estimate of the initial planet
    Mp = 4*np.pi/3 *5513.4* Rp**3
    #compute the equivalent semi-major axis
    df['aeq']=df['j'] * df['j']/(scipy.constants.G * Mp)
    df['mu']=scipy.constants.G*(Mp+df['m'])
    df['spec_orb_energy']=0.5*df['v']**2-df['mu']/df['r']
    df['e']=np.sqrt(1+2*df['spec_orb_energy']*df['j']**2/(df['mu']**2))
    mdisk=0
    mesc=0
    list_parts=[]
    list_tag2=[]
    for i in range(len(df)): #faut vectoriser ça en utilisant loc
        if df['r'][i]<Rp or df['aeq'][i]<Rp:
            list_parts.append('planet')
            list_tag2.append(0)
        elif df['e'][i]<1:
            list_parts.append('disk')
            list_tag2.append(1)
            mdisk+=df['m'][i]
        else:
            list_parts.append('escaping')
            mesc+=df['m'][i]
            list_tag2.append(2)
    df['part']=list_parts
    df['tag2']=list_tag2
    return df, mdisk, mesc

def converge_parts(df,Rpguess,precision):
    Mtot=df['m'].sum()
    Rtot=(3*Mtot/(5513.4*4*np.pi))**(1/3)
    err=1
    Rp=Rpguess
    list_df=[]
    errlist=[]
    while err > precision: #en fait la condition doit être sur deux valeurs successives
        element=define_parts(df,Rp)
        Mpafter=Mtot-element[1]-element[2]
        Rpafter=(3*Mpafter/(5513.4*4*np.pi))**(1/3)
        err=abs(Rpafter-Rp)/Rp
        Rp=Rpafter
        list_df.append(element)
        errlist.append(err)
    else:
        ls=list_df[len(list_df)-1]
        Mpafter=Mtot-ls[1]-ls[2]
        Rpafter=(3*Mpafter/(5513.4*4*np.pi))**(1/3)
        print(err)
        print(Rpafter)
        return ls
    
def bin_temperature_profile(df,qtty,dr,dr_th,gi_name): #to use on a already identified planet
    Rp = round(df['r'].max())
    r_min = round(df['r'].min(),-3)
    r_demi = np.arange(r_min-0.5*dr,Rp+0.5*dr,dr)
    r = np.arange(r_min,Rp+dr,dr)
    print(Rp,r_min,len(r_demi))
    T_binned = []
    T_table = df.loc[df['r']<r_demi[0],['r',qtty]]
    if qtty == "m":
        T_binned.append(df.loc[df['r']<=r_demi[0],[qtty]].squeeze().sum())
    else:
        T_binned.append(df.loc[df['r']<=r_demi[0],[qtty]].squeeze().mean())
    #print(df.loc[df['r']<r_demi[0],['temp']].mean())
    for i in range(len(r_demi)-1):
        T_table = pd.concat([T_table,df.loc[(df['r']>r_demi[i]) & (df['r']<r_demi[i+1]),['r',qtty]]],ignore_index=True)
        if qtty == "m":
            T_binned.append(df.loc[(df['r']>r_demi[i]) & (df['r']<=r_demi[i+1]),[qtty]].squeeze().sum())
        else:
            T_binned.append(df.loc[(df['r']>r_demi[i]) & (df['r']<=r_demi[i+1]),[qtty]].squeeze().mean())
    if qtty == "m":
        T_binned[0] = T_binned[0]
        for i in range(1,len(r_demi)):
            T_binned[i] = T_binned[i]+T_binned[i-1]
    else:
        T_binned = T_binned
        """
        T_table = pd.concat([T_table,df.loc[(df['r']>r_demi[i]) & (df['r']<r_demi[i+1]),['r',qtty]]],ignore_index=True)
        if qtty == "m":
            T_binned.append(df.loc[(df['r']>r_demi[i]) & (df['r']<r_demi[i+1]),[qtty]].squeeze().sum())
        else:
            T_binned.append(df.loc[(df['r']>r_demi[i]) & (df['r']<r_demi[i+1]),[qtty]].squeeze().mean())
        """
    dbin = pd.DataFrame()
    dbin['r'] = r
    dbin[qtty] = T_binned
    print(dbin)
    dbin = dbin.dropna()
    dbin = dbin.reset_index(drop=True)
    print(dbin)
    filename = path+gi_name+'_'+qtty+'_binned_'+str(int(dr/1000))+'km_'+str(int(dr_th/1000))+'km.dat'
    if os.path.exists(filename):
        os.remove(filename)
        print(f"File '{filename}' has been deleted.")
    else:
        print(f"The file '{filename}' does not exist.")
    dbin.to_csv(filename,sep='\t',mode='a',index=False)
    print(dbin['r'][0])
    r_smoothed = []
    T_smoothed = []
    for i in range(len(dbin)-1):
        l = np.arange(dbin['r'][i],dbin['r'][i+1],dr_th)
        r_smoothed.extend((l))
        m = np.linspace(dbin[qtty][i],dbin[qtty][i+1],len(l)+1)
        T_smoothed.extend(m[:-1])
    dc = pd.DataFrame()
    dc['r'] = r_smoothed
    dc[qtty] = T_smoothed
    filename = path+gi_name+'_'+qtty+'_smoothed_'+str(int(dr/1000))+'km_'+str(int(dr_th/1000))+'km.dat'
    if os.path.exists(filename):
        os.remove(filename)
        print(f"File '{filename}' has been deleted.")
    else:
        print(f"The file '{filename}' does not exist.")
    dc.to_csv(filename ,sep='\t',mode='a',index=False)
    return r_demi,T_binned,r_smoothed,T_smoothed

def entropy_match_temp_profile(f_ener, f_den, f_T, dr): #to use once the energy, density, temperature are smoothed
    df_ener = pd.read_csv(f_ener,sep='\t')
    df_den = pd.read_csv(f_den,sep='\t')
    df_T = pd.read_csv(f_T,sep='\t')
    df = pd.concat([df_ener, df_den['dens'], df_T['temp']], axis=1)
    print(df)
    V = []
    for i in range(len(dbin)-1):
        l = np.arange(dbin['r'][i],dbin['r'][i+1],dr_th)
        r_smoothed.extend((l))
        m = np.linspace(dbin[qtty][i],dbin[qtty][i+1],len(l)+1)
        T_smoothed.extend(m[:-1])
    return df

gi=pd.read_csv(pathinput+filename+'_identified_new.dat',sep='\t')
RHE = 519734442.6240584
Runknown = 6374242
planet=gi.loc[(gi['tag2']==0) & (gi['r'] <= RHE)]
rest = gi.loc[(gi['tag2']==0) & (gi['r'] > RHE)]
disk=gi.loc[gi['tag2']==1]
escaping=gi.loc[gi['tag2']==2]
fullplanet = gi.loc[gi['tag2']==0]

planet_iron_tar = planet.loc[planet['tag']==1]
planet_sil_tar = planet.loc[planet['tag']==0]
planet_iron_imp = planet.loc[planet['tag']==3]
planet_sil_imp = planet.loc[planet['tag']==2]

planet_iron = pd.concat([planet_iron_tar,planet_iron_imp])
planet_sil = pd.concat([planet_sil_tar, planet_sil_imp])

dr = 50000
dr_sm = 10000

df = bin_temperature_profile(planet,'temp',dr, dr_sm,filename)

df_iron=bin_temperature_profile(planet_iron,'m',dr, dr_sm,filename+'_iron')
df_sil=bin_temperature_profile(planet_sil,'m',dr, dr_sm,filename+'_sil')

df_m=bin_temperature_profile(fullplanet,'m',dr, dr_sm,filename)

bin_temperature_profile(planet,'dens',dr, dr_sm,filename)
bin_temperature_profile(planet,'int ener',dr, dr_sm,filename)
df_p = bin_temperature_profile(planet,'pres',dr, dr_sm,filename)

fig, ax = plt.subplots(1,1)
fig.set_size_inches(5.2,5.2*11/13)
ax.set_xlabel('r [m]')
ax.set_ylabel('m [kg]')
ax.plot(df_iron[0],df_iron[1],'.',color='k')
ax.plot(df_sil[0],df_sil[1],'.',color='goldenrod')
plt.savefig(path+filename+'_cumulativemass_binned_'+str(int(dr/1000))+'km_Rmax_profile.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
plt.close

fig, ax = plt.subplots(1,1)
fig.set_size_inches(5.2,5.2*11/13)
ax.set_xlabel('r [m]')
ax.set_ylabel('m [kg]')
ax.plot(df_iron[2],df_iron[3],'.',color='k',markersize=3)
ax.plot(df_sil[0],df_sil[1],'.',color='goldenrod',markersize=3)
plt.savefig(path+filename+'_cumulativemass_smoothed_'+str(int(dr/1000))+'km_Rmax_profile.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
plt.close

lsdf=[planet_iron_tar, planet_iron_imp,planet_sil_tar, planet_sil_imp]
lscolor=[ '#3182bd', '#fec44f','#a6bddb','#c994c7']
lslabel=['target silicate','impactor silicate','target water','impactor water']
lsd=[planet, disk, escaping]
lsdlabel=['planet', 'disk', 'escaping']

fig, ax = plt.subplots(1,1)
fig.set_size_inches(5.2,5.2*11/13)
axbox = ax.get_position()
for k in range(len(lsd)):
    ax.plot(lsd[k]['r']/1000,lsd[k]['temp'],'.',label=lsdlabel[k])
ax.set_xlabel('Radius [km]')
ax.set_ylabel('Temperature [K]')
ax.plot(df[0]/1000,df[1],color='k')
plt.xscale('log') 
plt.axvline(x = 9730, color = 'k', ls = '--')
plt.axvline(x = 23530, color = 'k', ls = ':')
plt.figlegend(bbox_to_anchor=[axbox.x0+0.5*axbox.width, axbox.y0-0.17], loc="center", bbox_transform=fig.transFigure,ncol=2)
plt.savefig(path+filename+'_temp_average_'+str(int(dr/1000))+'km_disk_profile_Rcmb_Rmax_.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
plt.close

fig, ax = plt.subplots(1,1)
fig.set_size_inches(5.2,5.2*11/13)
for k in range(len(lsdf)):
    ax.plot(lsdf[k]['r']/1000,lsdf[k]['temp'],'.',color=lscolor[k],label=lslabel[k])
ax.set_xlabel('Radius [km]')
ax.set_ylabel('Temperature [K]')
ax.plot(df[0]/1000,df[1],color='k')
ax.set_xlim(0, 25000)
plt.axvline(x = 9730, color = 'k', ls = '--')
plt.axvline(x = 23530, color = 'k', ls = ':')
plt.figlegend(bbox_to_anchor=[axbox.x0+0.5*axbox.width, axbox.y0-0.17], loc="center", bbox_transform=fig.transFigure,ncol=2)
plt.savefig(path+filename+'_temp_average_'+str(int(dr/1000))+'km_profile_Rcmb_Rmax.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
plt.close

fig, ax = plt.subplots(1,1)
fig.set_size_inches(5.2,5.2*11/13)
for k in range(len(lsdf)):
    ax.plot(lsdf[k]['r']/1000,lsdf[k]['pres']/1e9,'.',color=lscolor[k],label=lslabel[k])
ax.set_xlabel('Radius [km]')
ax.set_ylabel('Pressure [GPa]')
ax.plot(df_p[0]/1000,[float(item)/1e9 for item in df_p[1]],color='k')
ax.set_xlim(0, 25000)
plt.axvline(x = 9730, color = 'k', ls = '--')
plt.axvline(x = 23530, color = 'k', ls = ':')
plt.figlegend(bbox_to_anchor=[axbox.x0+0.5*axbox.width, axbox.y0-0.17], loc="center", bbox_transform=fig.transFigure,ncol=2)
plt.savefig(path+filename+'_pres_average_'+str(int(dr/1000))+'km_profile_Rcmb_Rmax.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
plt.close

fig, ax = plt.subplots(1,1)
fig.set_size_inches(5.2,5.2*11/13)
ax.set_xlabel('r [km]')
ax.set_ylabel('cumulative mass ')
#plt.axvline(x = 3.487e6, color = 'k', ls = '--')
plt.axvline(x = 3.865e6, color = 'k', ls = '--')
#plt.axvline(x = 3.45e6, color = 'k', ls = '--')
#plt.axvline(x = 3.83e6, color = 'k', ls = '--')
#plt.axvline(x = 5.08e6, color = 'k', ls = '--')
plt.axhline(y = 0.99, color = 'k', ls = '--')
plt.axhline(y = 0.01, color = 'k', ls = '--')
ax.plot(df_sil[2],df_sil[3]/df_sil[3][len(df_sil[3])-1],'.',color= '#3182bd',markersize=3)
ax.plot(df_iron[2],df_iron[3]/df_iron[3][len(df_iron[3])-1],'.',color='#fec44f',markersize=3)
plt.savefig(path+filename+'_cumulative_mass_sil_iron_smoothed_'+str(int(dr/1000))+'km_Rmax_profile.png',dpi=300,bbox_inches='tight',facecolor='white', transparent=False)
plt.close

omega = ((planet['jx'].sum()**2+planet['jy'].sum()**2+planet['jz'].sum()**2))**0.5/(planet['x']**2+planet['y']**2).sum()
print(omega, 2*np.pi/omega/(3600))

