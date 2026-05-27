# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 18:00:52 2024

@author: emil
"""
import matplotlib.pyplot as plt
class PrettyPlot:
    def setTitles(ax,title,xtitle,ytitle):
        '''
        Method for setting title and axis labels for a single axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis object to modify.
        title : str
            Figure title.
        xtitle : str
            Label for x-axis.
        ytitle : str
            Label for y-axis.

        Returns
        -------
        None
        '''
        ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        plt.tight_layout()


    def setTitlesV(ax,title,xtitle,ytitle1,ytitle2):
        '''
        Method for setting titles and labels for vertically stacked subplots.

        Parameters
        ----------
        ax : array-like of matplotlib.axes.Axes
            Array containing subplot axes.
        title : str
            Title for the top subplot.
        xtitle : str
            Label for shared x-axis.
        ytitle1 : str
            Label for upper subplot y-axis.
        ytitle2 : str
            Label for lower subplot y-axis.

        Returns
        -------
        None
        '''
        ax[0].set_title(title)
        ax[1].set_xlabel(xtitle)
        ax[0].set_ylabel(ytitle1)
        ax[1].set_ylabel(ytitle2)
        plt.tight_layout()
        
    def setTitlesH(ax,title,title2,xtitle,xtitle2,ytitle):
        '''
        Method for setting titles and labels for horizontally stacked subplots.

        Parameters
        ----------
        ax : array-like of matplotlib.axes.Axes
            Array containing subplot axes.
        title : str
            Title for left subplot.
        title2 : str
            Title for right subplot.
        xtitle : str
            Label for left subplot x-axis.
        xtitle2 : str
            Label for right subplot x-axis.
        ytitle : str
            Shared y-axis label.

        Returns
        -------
        None
        '''
        ax[0].set_title(title)
        ax[1].set_title(title2)
        ax[0].set_xlabel(xtitle)
        ax[1].set_xlabel(xtitle2)
        ax[0].set_ylabel(ytitle)
        plt.tight_layout()

    def makePretty(figSize,a=10,b=5):
        '''
        Method for applying plotting style settings to matplotlib figures.

        Parameters
        ----------
        figSize : float
            Scaling factor for plot appearance.
        a : float
            Figure width scaling factor.
        b : float
            Figure height scaling factor.

        Returns
        -------
        None
        '''
        plt.rc("font", family=["DejaVu Sans"])
        plt.rc("axes", labelsize=18*figSize, titlesize=22*figSize)
        plt.rc("xtick", labelsize=16*figSize, top=True, direction="in")
        plt.rc("ytick", labelsize=16*figSize, right=True, direction="in")
        plt.rc("legend", fontsize=12*figSize)
        plt.rc("figure", figsize=(a*figSize, b*figSize))
        plt.rc("lines", linewidth=2*figSize, markersize=8*figSize)
        plt.rc("errorbar",capsize=6*figSize)
    def shareXAxis(axes):
        '''
        Method for sharing the x-axis among multiple axes.

        Parameters
        ----------
        axes : array-like of matplotlib.axes.Axes
            Collection of axes objects.

        Returns
        -------
        None
        '''
        for ax in axes[1:]:
            ax.sharex(axes[0])
        
    def shareYAxis(axes):
        '''
        Method for sharing the y-axis among multiple axes.

        Parameters
        ----------
        axes : array-like of matplotlib.axes.Axes
            Collection of axes objects.

        Returns
        -------
        None
        '''
        for ax in axes[1:]:
            ax.sharey(axes[0])

    def shareAxis(axes):
        '''
        Method for sharing both x-axis and y-axis among multiple axes.

        Parameters
        ----------
        axes : array-like of matplotlib.axes.Axes
            Collection of axes objects.

        Returns
        -------
        None
        '''
        for ax in axes[1:]:
            ax.sharex(axes[0])
            ax.sharey(axes[0])
        
    def shareAxiss(axes):
        '''
        Method for sharing axes in a 2D grid of subplots.

        Parameters
        ----------
        axes : 2D array-like of matplotlib.axes.Axes
            Grid of subplot axes.

        Returns
        -------
        None
        '''
        for ax in axes:
            for a in ax:
                a.sharex(axes[0][0])
                a.sharey(axes[0][0])

    def plot_residuals(xs,ys,func,func_xs,popt,label="data",xlabel='x',ylabel='y',title='Residuals plot'):
        '''
        Method for plotting fitted data together with residuals.

        Parameters
        ----------
        xs : array-like
            x-values of data.
        ys : array-like
            y-values of data.
        func : callable
            Function used for fitting.
        popt : array-like
            Optimized fit parameters.
        label : str
            Label for plotted data.
        xlabel : str
            Label for x-axis.
        ylabel : str
            Label for y-axis.
        title : str
            Plot title.

        Returns
        -------
        None
        '''
        resi_fig, (ax_main, ax_resid) = plt.subplots(2, 1, sharex=True, figsize=(6, 7),
        gridspec_kw={"height_ratios": [2, 1],  # top plot is 3x taller
                     "hspace": 0.05})
        
        ax_main.plot(xs,ys,'x-',label=label)
        ax_main.plot(func_xs,func(func_xs,*popt), label='fit')
        ax_resid.plot(xs,ys-func(xs,*popt),'x')
        ax_resid.axhline(0,color='k')
    
        ax_resid.set_xlabel(xlabel)
        ax_resid.set_ylabel(r'$\Delta$')
        ax_main.set_ylabel(ylabel)
        ax_main.set_title(title)

    def savefig(filename):
        plt.savefig(f'/root/SR_Module/Plots/{filename}.svg', format="svg",bbox_inches='tight')