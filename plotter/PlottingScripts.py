import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
from numpy import mean
from math import sqrt,acos,cos,sin,pi,exp,log,isnan,atan2
from numpy import asarray
from root_pandas import read_root
from matplotlib import gridspec
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import quad
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from textwrap import wrap
import copy

def Tune1(df):
    return df['xsec_tune1_weight'].values

def CV(df):
    wgt = df['xsec_corr_weight'].values
    wgt[wgt == np.inf] = 1
    wgt = np.nan_to_num(wgt)
    wgt[wgt <= 0] = 1
    return wgt

def Pi0Wgt(df):
    wgt0 = CV(df)
    wgt1 = np.nan_to_num(df['pi0weight'].values)
    wgt1[wgt1 <=0] = 1
    return wgt0*wgt1

def Spline(df):
    wgt = df['spline_weight'].values
    wgt[wgt == np.inf] = 1
    wgt = np.nan_to_num(wgt)
    wgt[wgt <= 0] = 1
    return wgt

def poisson_errors(k, CL = 0.6827):
   # 1 Sig = 0.6827
   # 2 Sig = 0.9545
   # 3 Sig = 0.9973
   a = 1.0 - CL
   low, high = (stats.chi2.ppf(a/2, 2*k) / 2, stats.chi2.ppf(1-a/2, 2*k + 2) / 2)
   if k == 0:
       low = 0.0
   return k - low, high - k

 
class sampHist:

    def __init__(self,samp_df,samp_l,samp_c,samp_wind,samp_s):
        self._df = samp_df.copy()
        self._label = samp_l
        self._color = samp_c
        self._scale = samp_s
        self._wi = samp_wind
        if samp_wind == 0:
            self._wgt = np.ones(len(samp_df))
        if samp_wind == 1:
            self._wgt = CV(samp_df)
        if samp_wind == 2:
            self._wgt = Pi0Wgt(samp_df)

    def setweight(self,windex):
        if windex == 0:
            self._wgt = np.ones(len(self._df))
        if windex == 1:
            self._wgt = CV(self._df)
        if windex == 2:
            self._wgt = Pi0Wgt(samp_df)

    def dist(self,s_var):
        return self._df[s_var].values

    def cosdist(self,s_var):
        return np.cos(self._df[s_var].values)

    def applycut(self,s_cut):
        newhist = copy.deepcopy(self)
        newhist._df = self._df.query(s_cut)
        newhist.setweight(self._wi)
        return newhist

class SimpleHisto:

    def __init__ (self,df_df,f_scale,i_wgt,s_color,s_label):
        temp_df = df_df.copy()
        temp_df['myscale'] = f_scale
        self.mydf = temp_df
        self.mycolor = s_color
        self.mylabel = s_label
        self.iwgt = i_wgt

        self.mycut = 'Enu_1m1p > 0'

    def AddCut(self,s_cut):
        self.mycut  = s_cut

    def ClearCut(self):
        self.mycut = 'Enu_1m1p > 0'

    def GetHist(self,s_varname):
        temp_df = self.mydf.query(self.mycut)
        if self.iwgt == 2:
            myweight = Pi0Wgt(temp_df)
        elif self.iwgt == 1:
            myweight = CV(temp_df)
        elif self.iwgt == 0:
            myweight = np.ones(len(temp_df))
        return temp_df[s_varname].values,myweight,temp_df['myscale'].values,self.mycolor,self.mylabel

class StackedHisto:

    def  __init__ (self,a_df_mc,a_df_scale):
        self.mystrata = []
        self.stratalabel = []
        self.stratacolor = []
        self.mylayer = []
        self.layerlabel = []
        self.layercolor = []
        self.layeriwgt = []
        self.stratxweight = []

        temp_a_df = []
        for i in range(len(a_df_mc)):
            temp_df = a_df_mc[i].copy()
            temp_df['myscale'] = a_df_scale[i]
            temp_a_df.append(temp_df)
        if len(temp_a_df) > 0:
            self.mymc = pd.concat(temp_a_df)

        self.mycut = 'Enu_1m1p > 0'

    def ClearCut(self):
        self.mycut = 'Enu_1m1p > 0'

    def AddCut(self,s_cut):
        self.mycut  = s_cut

    def AddLayer(self,df_layer,df_scale,i_wgt,s_label,s_color):
        temp_df = df_layer.copy()
        temp_df['myscale'] = df_scale
        self.mylayer.append(temp_df)
        self.layerlabel.append(s_label)
        self.layercolor.append(s_color)
        self.layeriwgt.append(i_wgt)

    def AddStrata(self,s_strata,s_label,s_color,f_wgt=1.0):
        self.mystrata.append(s_strata)
        self.stratalabel.append(s_label)
        self.stratacolor.append(s_color)
        self.stratxweight.append(f_wgt)

    def GetHists(self,s_varname):
        a_vals = []        # (nxN)
        a_wgts = []
        a_scale = []
        a_cols = []        # (nx1)
        a_labels = []

        # first, run  through strata
        for i in range(len(self.mystrata)):
            # isolate this stratum with current cut
            temp_df = self.mymc.query(self.mystrata[i]+' and '+self.mycut)

            a_vals.append(temp_df[s_varname].values)
            
            a_wgts.append(Pi0Wgt(temp_df) * self.stratxweight[i])   # ah, this is if we want to scale an individual stratum
            a_scale.append(temp_df['myscale'].values)

            a_cols.append(self.stratacolor[i])
            a_labels.append(self.stratalabel[i])

        for i in range(len(self.mylayer)):
            temp_df = self.mylayer[i].query(self.mycut)

            a_vals.append(temp_df[s_varname].values)
            if self.layeriwgt[i] == 2:
                a_wgts.append(Pi0Wgt(temp_df))
            elif self.layeriwgt[i] == 1:
                a_wgts.append(CV(temp_df))
            elif self.layeriwgt[i] == 0:
                a_wgts.append(np.ones(len(temp_df)))
            a_scale.append(temp_df['myscale'].values)

            a_cols.append(self.layercolor[i])
            a_labels.append(self.layerlabel[i])

        return  np.asarray(a_vals),np.asarray(a_wgts),np.asarray(a_scale),a_cols,a_labels

    def ClearLayers(self):
        self.mylayer = []
        self.layerlabel = []
        self.layercolor = []
        self.layerwgt = []

    def ClearStrata(self):
        self.mystrata = []
        self.stratalabel = []
        self.stratacolor = []

    def DumpLayer(self,ind):
        #print('Layer:',self.layerlabel[ind])
        return self.mylayer[ind].query(self.mycut).copy()

    def DumpStratus(self,ind):
        return self.mymc.query(self.mystrata[ind]+' and '+self.mycut).copy()
        #print('Stratus:',self.self.stratalabel[ind])

class distVar:

    def __init__(self,s_name,n_range,s_label='',s_cov=''):
        self.myname = s_name
        self.myrange = n_range
        self.myscov = s_cov
        if s_label == '':
            self.mylabel = s_name
        else:
            self.mylabel = s_label


from scipy.interpolate import interp1d
from scipy.integrate import quad

def PmuGivenX(mu,x):

    # Returns the probability density value that the true mean is mu given that you
    # have observed x events. Actually employs a series approximation because the real formula
    # involves a very very large numerator and denominator and overflows even high precision
    # variables if your bins contain more than a few hundred events
    # Approximation is good to O(0.1%) or better at all values

    pi  = np.pi
    c   = [1.,-1./12,1./288,139./51840,-571./2488320,-163879./209018880]

    if x == 0:
        return np.exp(-mu)
    else:
        poly = sum(c[i]/x**(i+0.5) for i in range(len(c)))
        return 1/np.sqrt(2*pi)*np.exp(x+x*np.log(mu/x)-mu)*poly

def GetErrorsData(xobs,CL=0.6827):
    step    = 0.01
    upperBoundary = int(max(10,xobs+5*np.sqrt(xobs)))
    r = np.arange(0.01,upperBoundary,step)
    s    = PmuGivenX(r,xobs)*step
    PDF1 = interp1d(r,s,bounds_error=False,fill_value=0)
    PPF1 = interp1d(np.cumsum(s),r)
    xobs_low  = float(PPF1((1-CL)/2))
    xobs_high = float(PPF1(1-(1-CL)/2))
    return xobs_low,xobs_high


def distplot_single(myvar,nbins,myrange,stackedhists,stxcoord,m_cov,legpos=0,ymax=-1,fs=(16,11),consFac=0.0):

    #consFac is new and consolidates all MC strata whose contribution is less than a given percent threshold.
    # (my legend was too full of useless stuff)

    vals_mc = np.zeros(nbins)
    vals_mc_raw = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)
    yrr_mc_sys = np.zeros(nbins)
    ndof = 0

    a_labels_evts_draw = []
    a_labels_evts_failThresh = []
    a_cols_draw = []
    a_vals_draw = []
    a_wgts_draw = []
    a_vals_failThresh = []
    a_wgts_failThresh = []
    evts_failThresh = 0.0

    gh_vals,gh_wgts,gh_scale,gh_cols,gh_labels = stackedhists.GetHists(myvar)

    vals_mc_total,binedges = np.histogram(np.hstack(gh_vals),nbins,range=myrange,weights=np.hstack(gh_wgts)*np.hstack(gh_scale))
    mctot = vals_mc_total.sum()
    mcthresh = consFac * mctot

    for i in range(len(gh_vals)):
        h1_raw,_ = np.histogram(gh_vals[i],nbins,range=myrange,weights=gh_wgts[i])     # hist of raw event weights
        h1,_ = np.histogram(gh_vals[i],nbins,range=myrange,weights=gh_wgts[i]*gh_scale[i])

        vals_mc_raw += h1_raw
        vals_mc += h1

        scalesort = gh_scale[i].argsort()
        sorted_vals = gh_vals[i][scalesort[::-1]]
        sorted_scale = gh_scale[i][scalesort[::-1]]
        sorted_wgt = gh_wgts[i][scalesort[::-1]]
        for sc in np.unique(sorted_scale):
            subvals = sorted_vals[sorted_scale==sc]
            subwgts = sorted_wgt[sorted_scale==sc]
            subh1,_ = np.histogram(subvals,nbins,range=myrange,weights=subwgts)
            yerrsq_mc += np.power(np.sqrt(subh1)*sc,2)

        if(h1.sum() > mcthresh):
            a_vals_draw.append(gh_vals[i].copy())
            a_wgts_draw.append((gh_wgts[i]*gh_scale[i]).copy())
            a_labels_evts_draw.append(gh_labels[i]+' (%.2f)'%h1.sum())
            a_cols_draw.append(gh_cols[i])
        else:
            a_vals_failThresh.append(gh_vals[i].copy())
            a_wgts_failThresh.append((gh_wgts[i]*gh_scale[i]).copy())
            a_labels_evts_failThresh.append(gh_labels[i]+' (%.2f)'%h1.sum())
            evts_failThresh += h1.sum()

    if(len(a_vals_failThresh) > 0):
        #Ok. We've separated out strata and layers which fail our 1% threshold.
        # conbine the small contributions together:
        vals_failThresh = np.concatenate(a_vals_failThresh)
        wgts_failThresh = np.concatenate(a_wgts_failThresh)
        a_vals_draw.append(vals_failThresh)
        a_wgts_draw.append(wgts_failThresh)
        a_cols_draw.append('mediumspringgreen')
        a_labels_evts_draw.append('Other: (%.2f)'%evts_failThresh)
        print('Sub  1\% contributions:')
        print(a_labels_evts_failThresh)

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    #m_cov is fractional, so we multiply it by MC
    for i in range(nbins):
        for j in range(nbins):
            if(vals_mc[i] > 0 and vals_mc[j] > 0):
                m_cov[i][j] *= vals_mc[i]*vals_mc[j]
            else:
                m_cov[i][j] = 0

    # Normalization uncertainty:
    fNorm_squared = m_cov.sum() / np.power(vals_mc.sum(),2)
    print('Normalization Uncertainty:',np.sqrt(fNorm_squared))

    yerr_mc_sys = np.sqrt(np.diag(m_cov))
    yerr_mc_total = np.sqrt(np.diag(m_cov) + yerrsq_mc)

    fig,ax = plt.subplots(figsize=fs)
    fig.patch.set_alpha(1)

    ymax = (vals_mc+yerr_mc_sys).max()*1.4

    ax.set_ylim(0,ymax)
    ax.set_xlim(myrange)
    ax.set_ylabel('Events in Sample POT',fontsize=25)
    ax.set_title('MCC9 Data/MC',fontsize=30)

    ax.hist(a_vals_draw,nbins,range=myrange,weights=a_wgts_draw,color=a_cols_draw,stacked=True,linewidth=0,label=a_labels_evts_draw,edgecolor=None)
    ax.hist(np.hstack(a_vals_draw),nbins,range=myrange,weights=np.hstack(a_wgts_draw),histtype='step',zorder=20,color='black',linewidth=2)

    errboxes_tot = []
    errboxes_sys = []

    for i in range(len(bincenters)):
        rect0 = Rectangle((binedges[i],(vals_mc[i]-yerr_mc_sys[i])),binedges[i+1]-binedges[i],yerr_mc_sys[i]*2)
        errboxes_sys.append(rect0)
        rect1 = Rectangle((binedges[i],(vals_mc[i]-yerr_mc_total[i])),binedges[i+1]-binedges[i],yerr_mc_total[i]*2)
        errboxes_tot.append(rect1)
    pc_sys = PatchCollection(errboxes_tot,facecolor='red', alpha=.3,hatch='X',edgecolor='white')
    pc_sys_outline = PatchCollection(errboxes_tot,facecolor='none', alpha=.9,hatch='X',edgecolor='white',zorder=11)
    pc_tot = PatchCollection(errboxes_sys,facecolor=None,alpha=.1,hatch='/',zorder=12)
    ax.add_collection(pc_sys)
    ax.add_collection(pc_sys_outline)
    ax.add_collection(pc_tot)
    ax.hist(np.zeros(1),(1,2),facecolor=None,alpha=.1,hatch='//',zorder=0,label='MC Systematic Error')
    ax.hist(np.zeros(1),(1,2),facecolor='red', alpha=.3,hatch='X',edgecolor='white',label='MC Sys+Stat Error')
    ax.legend(loc='upper right',fontsize=20,frameon=False,ncol=3)

    return fig,ax


def distplot_wratio_dvar(dvar,nbins,StackedHists,datahist,m_cov=-1,legpos=0,ymax=-1,normshift=1,fs=(16,11),consFac=0.0,debug=False):
    return distplot_wratio(dvar.myname,nbins,dvar.myrange,StackedHists,datahist,dvar.mylabel,m_cov,legpos,ymax,normshift,fs,consFac,debug)

def distplot_wratio(myvar,nbins,myrange,stackedhists,datahist,stxcoord,m_cov=-1,legpos=0,ymax=-1,normshift=1,fs=(16,11),consFac=0.0,debug=False):

    #consFac is new and consolidates all MC strata whose contribution is less than a given percent threshold.
    # (my legend was too full of useless stuff)
    statsonly = False
    if np.array_equal(np.asarray(m_cov),np.asarray((-1))):
        statsonly = True
        m_cov = np.zeros((nbins,nbins))
        

    vals_mc = np.zeros(nbins)
    vals_mc_raw = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)
    yrr_mc_sys = np.zeros(nbins)
    ndof = 0

    a_labels_evts_draw = []
    a_labels_evts_failThresh = []
    a_cols_draw = []
    a_vals_draw = []
    a_wgts_draw = []
    a_vals_failThresh = []
    a_wgts_failThresh = []
    evts_failThresh = 0.0

    gh_vals,gh_wgts,gh_scale,gh_cols,gh_labels = stackedhists.GetHists(myvar)
    data_vals,data_wgt,data_scale,_,data_label = datahist.GetHist(myvar)

    vals_mc_total,binedges = np.histogram(np.hstack(gh_vals),nbins,range=myrange,weights=np.hstack(gh_wgts)*np.hstack(gh_scale))
    mctot = vals_mc_total.sum()
    mcthresh = consFac * mctot


    for i in range(len(gh_vals)):
        h1_raw,_ = np.histogram(gh_vals[i],nbins,range=myrange,weights=gh_wgts[i])     # hist of raw event weights
        h1,_ = np.histogram(gh_vals[i],nbins,range=myrange,weights=gh_wgts[i]*gh_scale[i])

        vals_mc_raw += h1_raw
        vals_mc += h1

        scalesort = gh_scale[i].argsort()
        sorted_vals = gh_vals[i][scalesort[::-1]]
        sorted_scale = gh_scale[i][scalesort[::-1]]
        sorted_wgt = gh_wgts[i][scalesort[::-1]]
        for sc in np.unique(sorted_scale):
            subvals = sorted_vals[sorted_scale==sc]
            subwgts = sorted_wgt[sorted_scale==sc]
            subh1,_ = np.histogram(subvals,nbins,range=myrange,weights=subwgts)
            yerrsq_mc += np.power(np.sqrt(subh1)*sc,2)

        if(h1.sum() > mcthresh):
            a_vals_draw.append(gh_vals[i].copy())
            a_wgts_draw.append((gh_wgts[i]*gh_scale[i]).copy())
            a_labels_evts_draw.append(gh_labels[i]+' (%.2f)'%h1.sum())
            a_cols_draw.append(gh_cols[i])
        else:
            a_vals_failThresh.append(gh_vals[i].copy())
            a_wgts_failThresh.append((gh_wgts[i]*gh_scale[i]).copy())
            a_labels_evts_failThresh.append(gh_labels[i]+' (%.2f)'%h1.sum())
            evts_failThresh += h1.sum()

    if(len(a_vals_failThresh) > 0):
        #Ok. We've separated out strata and layers which fail our 1% threshold.
        # conbine the small contributions together:
        vals_failThresh = np.concatenate(a_vals_failThresh)
        wgts_failThresh = np.concatenate(a_wgts_failThresh)
        a_vals_draw.append(vals_failThresh)
        a_wgts_draw.append(wgts_failThresh)
        a_cols_draw.append('mediumspringgreen')
        a_labels_evts_draw.append('Other: (%.2f)'%evts_failThresh)
        print('Sub  1\% contributions:')
        print(a_labels_evts_failThresh)

    vals_data_raw,_ = np.histogram(data_vals,nbins,range=myrange,weights=data_wgt)
    vals_data,_ = np.histogram(data_vals,nbins,range=myrange,weights=np.multiply(data_wgt,data_scale))
    bincenters = np.diff(binedges)/2 + binedges[:-1]

    #jarretbars
    a_obslo = []
    a_obshi = []
    for i in range(nbins):
        obslo,obshi = poisson_errors(vals_data[i])
        a_obshi.append(obshi)
        a_obslo.append(obslo)

    #m_cov is fractional, so we multiply it by MC
    for i in range(nbins):
        for j in range(nbins):
            if(vals_mc[i] > 0 and vals_mc[j] > 0):
                m_cov[i][j] *= vals_mc[i]*vals_mc[j]
            else:
                m_cov[i][j] = 0

    # Normalization uncertainty:
    fNorm_squared = m_cov.sum() / np.power(vals_mc.sum(),2)
    print('Normalization Uncertainty:',np.sqrt(fNorm_squared))

    yerr_mc_sys = np.sqrt(np.diag(m_cov))
    yerr_mc_total = np.sqrt(np.diag(m_cov) + yerrsq_mc)

    yerrsq_data = np.zeros(nbins)
    for i in range(nbins):
        if vals_mc_raw[i] > 0:
            ndof += 1
            if vals_data[i] > 0:
                yerrsq_data[i] += (3.0*vals_data[i]*vals_mc[i]*normshift)/(vals_mc[i]*normshift+2.0*vals_data[i])
            else:
                yerrsq_data[i] += vals_mc[i]*normshift/2.0
            m_cov[i][i] += yerrsq_data[i] + yerrsq_mc[i]
        else:
            m_cov[i][i] += 999

    yerr_data = np.sqrt(yerrsq_data)

    if debug:
        print('mc stat error',np.divide(np.sqrt(yerrsq_mc),vals_mc))
        print('data stat error CNP',np.divide(np.sqrt(yerrsq_data),vals_data))
        print('mc_sys_error',np.divide(np.sqrt(yerr_mc_sys),vals_mc))
        for i in range(20):
            print('cov_collapsed:',m_cov[i][:])
    
    er_rat_line = np.zeros(nbins)
    er_rat_line_sys = np.zeros(nbins)
    er_rat_dotshi = np.zeros(nbins)
    er_rat_dotslo = np.zeros(nbins)

    for i in range(nbins):
        if vals_mc[i] > 0:
            er_rat_line[i] = yerr_mc_total[i]/float(vals_mc[i])
            er_rat_line_sys[i] = yerr_mc_sys[i]/float(vals_mc[i])
            er_rat_dotshi[i] = a_obshi[i]/float(vals_mc[i])
            er_rat_dotslo[i] = a_obslo[i]/float(vals_mc[i])

    chisq = 0.0
    invcov = np.linalg.inv(m_cov)

    # calc chi2
    for i in range(nbins):
        for j in range(nbins):
            if vals_mc_raw[i] > 0 and vals_mc_raw[j] > 0:
                chisq += (vals_mc[i] - vals_data[i]) * (vals_mc[j] - vals_data[j]) * invcov[i][j]
    pval = 1 - stats.chi2.cdf(chisq, ndof)
    print(chisq/float(ndof),pval)

    fig,ax = plt.subplots(figsize=fs)
    fig.patch.set_alpha(1)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, .75])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    ymax = max((vals_data+np.asarray(a_obshi)).max(),(vals_mc+yerr_mc_sys).max())*1.4

    ax0.set_ylim(0,ymax)
    ax0.set_xlim(myrange)
    ax1.set_ylim(0,2)
    ax1.set_xlim(myrange)
    ax1.set_xlabel(stxcoord,fontsize=30)
    ax0.set_ylabel('Events in Full Sample',fontsize=30)
    ax1.set_ylabel('Data/MC',fontsize=30)
    ax0.set_title('MCC9 Data/MC',fontsize=35)
    ax0.tick_params(axis='x', labelsize=25)
    ax0.tick_params(axis='y', labelsize=25)
    ax1.tick_params(axis='x', labelsize=25)
    ax1.tick_params(axis='y', labelsize=25)
    ax0.hist(a_vals_draw,nbins,range=myrange,weights=a_wgts_draw,color=a_cols_draw,stacked=True,linewidth=0,label=a_labels_evts_draw,edgecolor=None)
    ax0.hist(np.hstack(a_vals_draw),nbins,range=myrange,weights=np.hstack(a_wgts_draw),histtype='step',zorder=20,color='black',linewidth=2)

    ax0.errorbar(bincenters,vals_data,fmt='o',yerr=(a_obslo,a_obshi),color='black',capsize=5,label=data_label+' (%i)'%vals_data.sum(),markersize=8,zorder=20,elinewidth=2)
    ax0.errorbar(bincenters,vals_data,fmt='o',color='white',zorder=19,markersize=16)

    errboxes_tot = []
    errboxes_sys = []

    for i in range(len(bincenters)):
        rect0 = Rectangle((binedges[i],(vals_mc[i]-yerr_mc_sys[i])),binedges[i+1]-binedges[i],yerr_mc_sys[i]*2)
        errboxes_sys.append(rect0)
        rect1 = Rectangle((binedges[i],(vals_mc[i]-yerr_mc_total[i])),binedges[i+1]-binedges[i],yerr_mc_total[i]*2)
        errboxes_tot.append(rect1)
    pc_sys = PatchCollection(errboxes_tot,facecolor='red', alpha=.3,hatch='X',edgecolor='white')
    pc_sys_outline = PatchCollection(errboxes_tot,facecolor='none', alpha=.9,hatch='X',edgecolor='white',zorder=11)
    pc_tot = PatchCollection(errboxes_sys,facecolor=None,alpha=.1,hatch='/',zorder=12)
    ax0.add_collection(pc_sys)
    ax0.add_collection(pc_sys_outline)
    ax0.add_collection(pc_tot)
    if not statsonly:
        ax0.hist(np.zeros(1),(1,2),facecolor=None,alpha=.1,hatch='//',zorder=0,label='MC Systematic Error')
        ax0.hist(np.zeros(1),(1,2),facecolor='red', alpha=.3,hatch='X',edgecolor='white',label='MC Sys+Stat Error')
    else:
        ax0.hist(np.zeros(1),(1,2),facecolor='red', alpha=.3,hatch='X',edgecolor='white',label='MC Stat Error')
    ax0.legend(loc='upper right',fontsize=20,frameon=False,ncol=3)

    errboxes_rat_tot = []
    errboxes_rat_sys = []
    for i in range(len(er_rat_dotshi)):
        rect0 = Rectangle((binedges[i],1-er_rat_line[i]),binedges[i+1]-binedges[i],er_rat_line[i]*2)
        errboxes_rat_tot.append(rect0)
        rect1 = Rectangle((binedges[i],1-er_rat_line_sys[i]),binedges[i+1]-binedges[i],er_rat_line_sys[i]*2)
        errboxes_rat_sys.append(rect1)

    pc_rat_tot = PatchCollection(errboxes_rat_tot, facecolor='red', alpha=.3)
    pc_rat_sys = PatchCollection(errboxes_rat_sys, facecolor=None, alpha=.1,hatch='/',zorder=12)
    ax1.add_collection(pc_rat_tot)
    ax1.add_collection(pc_rat_sys)
    ax1.hist(np.zeros(1),(1,2),facecolor=None,alpha=.1,hatch='//',zorder=0)

    ax1.hist(np.zeros(1),(1,2),facecolor='red',alpha=.3,zorder=0)
    ax1.errorbar(bincenters,np.true_divide(vals_data,vals_mc),yerr=(er_rat_dotshi,er_rat_dotslo),fmt='o',color='maroon',capsize=0,markersize=8,elinewidth=2)

    #ax1.legend(loc='lower right',fontsize=15,frameon=False)
    if debug:
        print('Data/MC')
        print(np.true_divide(vals_data,vals_mc))
        print(bincenters)
        print('data:',vals_data)
        print('pred:',vals_mc)
        print('error:',yerr_mc_total)
    
    ax1.axhline(1,color='black',linestyle=':')
    ax0.annotate(r'$\sum$data/$\sum$MC = %.2f $\pm$ %.2f'%(vals_data.sum()/float(vals_mc.sum()),np.sqrt(fNorm_squared)),xy=(.01,.92),xycoords='axes fraction',fontsize=20,bbox=dict(boxstyle="square", fc="ghostwhite",alpha=.8))

    ax1.annotate('p-value: %.3f\n'%pval+r'$\chi^2_{CNP}/%i  (dof)$: %.3f'%(ndof,chisq/float(ndof)),xy=(.85,.75), xycoords='axes fraction',fontsize=20,bbox=dict(boxstyle="round4", fc="w",alpha=.9),zorder=30)

    plt.tight_layout()

    return fig,ax0,ax1,pval


def truthplot(myvar,nbins,myrange,predhists,stxcoord,legpos=0,ymax=-1,normshift=1):

    fig.patch.set_alpha(1)
    vals_mc = np.zeros(nbins)
    vals_troof = np.zeros(nbins)
    yerrsq_mc = np.zeros(nbins)

    a_labels_evts = []

    for i in range(len(predhists)):
        h1_raw,binedges = np.histogram(predhists[i].dist(myvar),nbins,range=myrange,weights=predhists[i]._wgt)     # hist of raw event weights
        h1 = h1_raw * predhists[i]._scale
        a_scaledweight.append(predhists[i]._wgt*predhists[i]._scale)

        yerrsq_mc += np.power(np.sqrt(h1_raw)*predhists[i]._scale,2)

        vals_mc += h1
        a_labels_evts.append(predhists[i]._label+' (%.2f)'%h1.sum())

        htroof,_ = np.histogram(predhists[i].dist('MC_energyInit'),nbins,range=myrange,weights=predhists[i]._wgt)
        vals_troof += htroof* predhists[i]._scale


    yerr_mc = np.sqrt(yerrsq_mc)

    bincenters = np.diff(binedges)/2 + binedges[:-1]

    if ymax < 0: ymax = vals_mc.max()*1.35

    ax.set_ylim(0,ymax)
    ax.set_xlim(myrange)
    ax.set_ylabel('Events in Sample',fontsize=30)
    ax.set_title('MCC9 MC/Truth',fontsize=35)

    ax.hist(gh_vals,nbins,range=myrange,weights=a_scaledweight,color=gh_cols,stacked=True,linewidth=0,label=a_labels_evts,edgecolor=None)
    plt.scatter(bincenters,vals_troof,color='black',label='Truth',zorder=10)

    ax.legend(loc='upper right',fontsize=15,frameon=False,ncol=3)

    return ax

def mcplot(_vartest,_myrange,_mybins,_predhists,_obshist):
    _pvals,_,_,_pwgts,_ = _predhists.GetHists(_vartest)
    _ovals,_,_,_owgts,_ = _obshist.GetHist(_vartest)

    phist = np.zeros(_mybins)
    for i in range(len(_pvals)):
        _ph,_ = np.histogram(_pvals[i],bins=_mybins,range=_myrange,weights=_pwgts[i])
        phist += _ph

    ohist,binedges = np.histogram(_ovals,bins=_mybins,range=_myrange,weights=_owgts)
    bincenters = np.diff(binedges)/2 + binedges[:-1]

    fig,ax = plt.subplots(figsize=(16,11))
    plt.scatter(bincenters,np.true_divide(phist,ohist,where=ohist!=0))
    ax.set_ylabel('Raw MC / Data',fontsize=25)
    ax.set_xlabel(_vartest)
    print('MC raw',phist.sum())
    print('Data raw',ohist.sum())


def DrawMatrix(cov,nbins,myrange):

    binedges = np.linspace(myrange,nbins)
    X, Y = np.meshgrid(binedges,binedges)
    fig,ax = plt.subplots(figsize=(10,10))
    crat = ax.pcolormesh(X, Y,cov,cmap='cool')
    cbar = fig.colorbar(crat)


def bless_MC_labels(row):
    mclabel = ''
    intlabel = ''
    parentlabel = ''
    pizero = [1090,1086,1090,1080,1015,1013,1011,1008,1006,1004]
    piplusminus = [1085,1079,1032,1017,1014,1007,1005,1003,1028,1021,1016,1012,1010,1009]

    if abs(row['nu_pdg']) == 12:
        intlabel = 'nue'
    elif abs(row['nu_pdg']) == 14:
        intlabel = 'numu'

    if not (row['MC_nproton']==1 and row['MC_nlepton']==1):
        return 'nLmP'
    elif not 0 < row['MC_scedr'] <= 5.0:
        return 'offvtx'
    elif not abs((row['MC_energyInit']-row['Enu_1m1p'])/row['MC_energyInit']) < 0.2:
        return 'badreco'
    else:
        if row['nu_interaction_type'] == 1001:
            mclabel = 'CCQE'
        elif row['nu_interaction_type'] == 1000:
            mclabel = 'MEC'
        elif row['nu_interaction_type'] in pizero:
            mclabel = 'pizero'
        elif row['nu_interaction_type'] in piplusminus:
            mclabel = 'piplusminus'
        else:
            mclabel = 'other'

    return '%s_%s'%(intlabel,mclabel)

def bless_int_labels(row):
    mclabel = ''
    intlabel = ''
    parentlabel = ''
    pizero = [1090,1086,1090,1080,1015,1013,1011,1008,1006,1004]
    piplusminus = [1085,1079,1032,1017,1014,1007,1005,1003,1028,1021,1016,1012,1010,1009]

    if abs(row['nu_pdg']) == 12:
        intlabel = 'nue'
    elif abs(row['nu_pdg']) == 14:
        intlabel = 'numu'

    if row['nu_interaction_type'] == 1001:
        mclabel = 'CCQE'
    elif row['nu_interaction_type'] == 1000:
        mclabel = 'MEC'
    elif row['nu_interaction_type'] in pizero:
        mclabel = 'pizero'
    elif row['nu_interaction_type'] in piplusminus:
        mclabel = 'piplusminus'
    else:
        mclabel = 'other'

    return '%s_%s'%(intlabel,mclabel)


def dist2d_statsonly(dvar1,dvar2,nbins1,nbins2,obsHist,mcHist):

    vals_obs_var1,wgt_obs,scale_obs,_,_ = obsHist.GetHist(dvar1.myname)
    vals_obs_var2,_,_,_,_ = obsHist.GetHist(dvar2.myname)

    Hobs, xedges, yedges = np.histogram2d(vals_obs_var1,vals_obs_var2,[nbins1,nbins2],[dvar1.myrange,dvar2.myrange],weights=wgt_obs*scale_obs)

    a_vals_pred_var1,a_wgt_pred,a_scale_pred,_,_ = mcHist.GetHists(dvar1.myname)
    a_vals_pred_var2,_,_,_,_ = mcHist.GetHists(dvar2.myname)

    vals_pred_var1 = np.concatenate(a_vals_pred_var1)
    vals_pred_var2 = np.concatenate(a_vals_pred_var2)
    wgt_pred = np.concatenate(a_wgt_pred)
    scale_pred = np.concatenate(a_scale_pred)

    Hpred,xedges, yedges = np.histogram2d(vals_pred_var1,vals_pred_var2,[nbins1,nbins2],[dvar1.myrange,dvar2.myrange],weights=wgt_pred*scale_pred)

    Hobs = Hobs.T
    Hpred = Hpred.T
    X, Y = np.meshgrid(xedges, yedges)

    #Now we'll make a matrix of chi2 contribution
    Hchi2 = np.zeros((nbins1,nbins2))

    totchi2 = 0.0
    for i in range(nbins1):
        for j in range(nbins2):
            if Hobs[i][j] > 0:
                err_cnp = (3.0*Hobs[i][j]*Hpred[i][j])/(Hpred[i][j]+2.0*Hobs[i][j])
            else:
                err_cnp = Hpred[i][j]/2.0
            Hchi2[i][j] = np.power(Hobs[i][j]-Hpred[i][j],2)/err_cnp
            totchi2 += Hchi2[i][j]

    fig,ax = plt.subplots(figsize=(10,10))
    ax.set_xlabel(dvar1.mylabel,fontsize=20)
    ax.set_ylabel(dvar2.mylabel,fontsize=20)

    crat = ax.pcolormesh(X, Y,Hchi2,cmap='cool')
    cbar = fig.colorbar(crat)
    cbar.ax.set_ylabel(r'$\chi^2$ Contribution', rotation=270,fontsize=15)

    plt.title(r'$\chi^2$',fontsize=25)
    plt.xlim(dvar1.myrange)
    plt.ylim(dvar2.myrange)

    return Hobs,Hpred,totchi2

# Cov matrices from detsys pickle pack. all tested and good ass of July 2
class Cov:
    def __init__(self,s_detsyspickle):

        with open(s_detsyspickle,'rb') as handle: (self.a_overlap_sys,self.a_cv_sys,self.s_detsyslist) = pickle.load(handle)
        self.s_cuts = 'Enu_1m1p > 0'
        self.s_cuts_cv = 'Enu_1m1p_cv > 0'

    def SetCuts(self,s_cuts,s_cuts_cv):
        self.s_cuts = s_cuts
        self.s_cuts_cv = s_cuts_cv

    def Nominal(self,dvar,nbins,draw=False):

        cov_tru = np.zeros((nbins,nbins))
        for sysi in range(len(self.a_overlap_sys)):
            myvardf = self.a_overlap_sys[sysi].query(self.s_cuts)
            myvarcv = self.a_cv_sys[sysi].query(self.s_cuts_cv)

            var_sys = myvardf[dvar.myname]
            var_cv = myvarcv[dvar.myname+'_cv']

            hCV,binedges = np.histogram(var_cv,bins=nbins,range=dvar.myrange)
            h0,_ = np.histogram(var_sys,bins=nbins,range=dvar.myrange)

            for i in range(nbins):
                for j in range(nbins):
                    cov_tru[i][j] += (h0[i]-hCV[i])*(h0[j]-hCV[j])/(hCV[i]*hCV[j])

        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(binedges,binedges)
            crat_tru = ax.pcolormesh(X, Y,cov_tru,cmap='cool')
            cbar = fig.colorbar(crat_tru)
            ax.set_xlabel(dvar.mylabel,fontsize=20)
            ax.set_ylabel(dvar.mylabel,fontsize=20)

        return cov_tru

    def Flat(self,dvar,nbins,draw=False):

        flatsys = 0;
        for sysi in range(len(self.a_overlap_sys)):
            myvardf = self.a_overlap_sys[sysi].query(self.s_cuts)
            myvarcv = self.a_cv_sys[sysi].query(self.s_cuts_cv)

            var_sys = myvardf[dvar.myname]
            var_cv = myvarcv[dvar.myname+'_cv']

            hCV,binedges = np.histogram(var_cv,bins=1,range=dvar.myrange)
            h0,_ = np.histogram(var_sys,bins=1,range=dvar.myrange)

            flatsys += (h0.sum()-hCV.sum())*(h0.sum()-hCV.sum())/(hCV.sum()*hCV.sum())

        if(draw):
            print('Flatsys:',flatsys,' - (',np.sqrt(flatsys),' frac error)')
        return np.ones((nbins,nbins))*flatsys

    def Polyfit(self,dvar,nbins,draw=False):

        cov_poly = np.zeros((nbins,nbins))
        for sysi in range(len(self.a_overlap_sys)):
            myvardf = self.a_overlap_sys[sysi].query(self.s_cuts)
            myvarcv = self.a_cv_sys[sysi].query(self.s_cuts_cv)

            var_sys = myvardf[dvar.myname]
            var_cv = myvarcv[dvar.myname+'_cv']

            hCV,binedges = np.histogram(var_cv,bins=nbins,range=dvar.myrange)
            h0,_ = np.histogram(var_sys,bins=nbins,range=dvar.myrange)
            bincenters = np.diff(binedges)/2 + binedges[:-1]

            truRat = np.true_divide(h0,hCV,out=np.ones_like(bincenters),where=hCV!=0)

            # get polyfit degs
            aics = []
            degs = []
            for deg in range(min(nbins-2,int(nbins/2))):
                params = deg + 1
                polyRat = np.polyfit(bincenters, truRat, deg)
                fRat = np.poly1d(polyRat)

                # now calculate chi2 for fit
                yerr_rat = np.true_divide(np.sqrt(fRat(bincenters)*hCV),hCV,out=np.zeros_like(bincenters),where=hCV!=0)
                chi2_fit = np.power(np.true_divide(fRat(bincenters)-truRat,yerr_rat),2).sum()
                aic = chi2_fit + 2*params + 2*params*(params+1)/float(nbins-params-1)
                aics.append(aic)
                degs.append(deg)

            polyterms = degs[np.argmin(aics)]
            print(self.s_detsyslist[sysi],'Polyfit Degrees:',polyterms,aics[np.argmin(aics)])
            polyRat = np.polyfit(bincenters, np.true_divide(h0,hCV,where=hCV!=0), polyterms)
            fRat = np.poly1d(polyRat)
            h0_fit = fRat(bincenters)*hCV

            for i in range(nbins):
                for j in range(nbins):
                    cov_poly[i][j] += (h0_fit[i]-hCV[i])*(h0_fit[j]-hCV[j])/(hCV[i]*hCV[j])

        if(draw):
            fig,ax = plt.subplots(figsize=(9,9))
            X, Y = np.meshgrid(binedges,binedges)
            crat_poly = ax.pcolormesh(X, Y,cov_poly,cmap='cool')
            cbar = fig.colorbar(crat_poly)
            ax.set_xlabel(dvar.mylabel,fontsize=20)
            ax.set_ylabel(dvar.mylabel,fontsize=20)

        return cov_poly




def datadataplot_wratio_dvar(dvar,nbins,datahist1,datahist2,legpos=0,ymax=-1,normshift=1,fs=(16,11)):
    return datadataplot_wratio(dvar.myname,nbins,dvar.myrange,datahist1,datahist2,dvar.mylabel,legpos,ymax,normshift,fs)

def datadataplot_wratio(myvar,nbins,myrange,datahist1,datahist2,stxcoord,legpos=0,ymax=-1,normshift=1,fs=(16,11)):

    ndof = 0

    data1_vals,data1_wgt,data1_scale,_,data1_label = datahist1.GetHist(myvar)
    data2_vals,data2_wgt,data2_scale,_,data2_label = datahist2.GetHist(myvar)

    vals_data1_raw,binedges = np.histogram(data1_vals,nbins,range=myrange,weights=data1_wgt)
    vals_data1,_ = np.histogram(data1_vals,nbins,range=myrange,weights=np.multiply(data1_wgt,data1_scale))
    vals_data2_raw,_ = np.histogram(data2_vals,nbins,range=myrange,weights=data2_wgt)
    vals_data2,_ = np.histogram(data2_vals,nbins,range=myrange,weights=np.multiply(data2_wgt,data2_scale))
    bincenters = np.diff(binedges)/2 + binedges[:-1]

    #jarretbars
    a_obslo1 = []
    a_obshi1 = []
    a_obslo2 = []
    a_obshi2 = []
    for i in range(nbins):
        obslo,obshi = poisson_errors(vals_data1[i])
        a_obshi1.append(obshi-vals_data[i])
        a_obslo1.append(vals_data[i]-obslo)
        obslo,obshi = poisson_errors(vals_data2[i])
        a_obshi2.append(obshi)
        a_obslo2.append(obslo)

    #m_cov is fractional, so we multiply it by MC
    m_cov = np.zeros((nbins,nbins))

    yerrsq_data = np.zeros(nbins)
    for i in range(nbins):
        if vals_data1[i] > 0 or  vals_data2[i] > 0:
            ndof += 1
            m_cov[i][i] += vals_data1[i] + vals_data2[i]
        else:
            m_cov[i][i] += 999

    er_rat_linehi = np.zeros(nbins)
    er_rat_linelo = np.zeros(nbins)
    er_rat_dotshi = np.zeros(nbins)
    er_rat_dotslo = np.zeros(nbins)

    for i in range(nbins):
        if vals_data2[i] > 0:
            er_rat_linehi[i] = a_obshi2/float(vals_data2[i])
            er_rat_linelo[i] = a_obslo2/float(vals_data2[i])
            er_rat_dotshi[i] = a_obshi1[i]/float(vals_data2[i])
            er_rat_dotslo[i] = a_obslo1[i]/float(vals_data2[i])

    chisq = 0.0
    invcov = np.linalg.inv(m_cov)

    # calc chi2
    for i in range(nbins):
        for j in range(nbins):
            if vals_data2[i] > 0 or vals_data1[i] > 0:
                chisq += (vals_data1[i] - vals_data2[i]) * (vals_data1[j] - vals_data2[j]) * invcov[i][j]
    pval = 1 - stats.chi2.cdf(chisq, ndof)
    print(chisq/float(ndof),pval)

    fig,ax = plt.subplots(figsize=fs)
    fig.patch.set_alpha(1)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, .75])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    ymax = max((vals_data1+np.asarray(a_obshi1)).max(),(vals_data2+a_obshi2).max())*1.4

    ax0.set_ylim(0,ymax)
    ax0.set_xlim(myrange)
    ax1.set_ylim(0,2)
    ax1.set_xlim(myrange)
    ax1.set_xlabel(stxcoord,fontsize=25)
    ax0.set_ylabel('Events in Sample',fontsize=25)
    ax1.set_ylabel('Data1/Data2',fontsize=25)
    ax0.set_title('MCC9 Data/Data',fontsize=30)

    ax0.errorbar(bincenters,vals_data1,fmt='o',yerr=(a_obslo1,a_obshi1),color='orange',capsize=5,label=data1_label+' (%i)'%vals_data1.sum(),markersize=8,zorder=20,elinewidth=2)
    ax0.errorbar(bincenters,vals_data2,fmt='o',yerr=(a_obslo2,a_obshi2),color='blue',capsize=10,label=data1_label+' (%i)'%vals_data1.sum(),markersize=16,zorder=19,elinewidth=4)

    ax0.legend(loc='upper right',fontsize=20,frameon=False,ncol=3)

    errboxes_rat_tot = []
    for i in range(len(er_rat_dotshi)):
        rect0 = Rectangle((binedges[i],1-er_rat_linelo[i]),binedges[i+1]-binedges[i],er_rat_linehi[i]+er_rat_linelo[i])
        errboxes_rat_tot.append(rect0)

    pc_rat_tot = PatchCollection(errboxes_rat_tot, facecolor='blue', alpha=.3)
    ax1.add_collection(pc_rat_tot)
    ax1.errorbar(bincenters,np.true_divide(vals_data1,vals_data2),yerr=(er_rat_dotshi,er_rat_dotslo),fmt='o',color='orange',capsize=0,markersize=8,elinewidth=2)

    ax1.axhline(1,color='black',linestyle=':')
    ax0.annotate(r'$\sum$data1/$\sum$data2 = %.2f'%(vals_data1.sum()/float(vals_data2.sum())),xy=(.01,.92),xycoords='axes fraction',fontsize=20,bbox=dict(boxstyle="square", fc="ghostwhite",alpha=.8))

    ax1.annotate('p-value: %.3f\n'%pval+r'$\chi^2_{CNP}/%i  (dof)$: %.3f'%(ndof,chisq/float(ndof)),xy=(.85,.75), xycoords='axes fraction',fontsize=20,bbox=dict(boxstyle="round4", fc="w",alpha=.9),zorder=30)

    plt.tight_layout()

    return fig,ax0,ax1,pval
