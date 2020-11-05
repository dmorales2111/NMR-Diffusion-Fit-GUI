#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diffusion curve-fitting tool, version 1.0.0

@author: daniel
"""

import tkinter as tk
import tkinter.font as tkfont
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import nmrglue as ng
from scipy.optimize import curve_fit
from tkinter import filedialog as fd
import os

class MainApplication():
    def __init__(self,master):
        self.master = master
        master.title('Welcome to Dan\'s diffusion fitting tool! Ver 1.0.0')
        pad=3
        master.geometry("{0}x{1}+0+0".format(
            master.winfo_screenwidth()-pad, master.winfo_screenheight()-pad))

        master.resizable(True, True)
        self.textbox = tk.Text(master)

        # colors for the text box, depending on the output
        self.textbox.tag_config('open', background = 'violet',
                                foreground = 'black')
        self.textbox.tag_config('fit', background = 'purple',
                                foreground = 'white')
        self.textbox.tag_config('clear', background = 'black', foreground = 'white')


        #sub windows in program
        self.welcome = WelcomeFrame(master, tk.RAISED, 'grey', self.textbox)
        self.data = DataWindow(master, 200, 200, tk.RAISED, 'grey', self.textbox)


        self.welcome.grid(row = 0, column = 0)
        self.data.grid(row = 0, column = 1, padx = 20)
        self.textbox.grid(row = 1, column = 0, columnspan = 2, pady = 10,
                          sticky = 'NSEW'
                          )


class WelcomeFrame(tk.Frame):
    def __init__(self, master, relief, bg, box):
        super(WelcomeFrame,self).__init__(master)
        self.master = master
        self.relief = relief
        self.box = box
        self.hintlabel = tk.Label(self,
                                  text = "Please input your directory with " +
                                  'your data (Directory should contain ' +
                                  'procpar and integ_series from SpinWorks):',
                                  font = tkfont.Font(size = 15)
                                  )
        self.hintlabel.grid(row = 1, column = 0, pady = 20, columnspan = 2)
        self.dir_sub_button = tk.Button(self,text = 'Open Folder',
                                        command = self.open
                                        )
        self.dir_sub_button.grid(row = 2, column = 0, sticky = 'NSEW')
        self.quitbutton = tk.Button(self,
                                    text = 'Quit',
                                    command = self.master.destroy
                                    )
        self.quitbutton.grid(row =2, column = 1, sticky = 'NSEW')
    def open(self):
        global direct
        direct = fd.askdirectory()
        self.box.insert(tk.END, 'Using files found in ' + direct + '\n' + '\n',
                        'open'
                        )

class DataWindow(tk.Frame):
    def __init__(self,master, width, height, relief, bg, box):
        super(DataWindow,self).__init__(master)
        self.master = master
        self.width = width
        self.height = height
        self.relief = relief
        self.box = box

        self.f = Figure(figsize=(4,4),dpi=150)
        self.ax = self.f.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row = 0, column = 0, columnspan=2,
                                         sticky = 'NSEW')

        self.v = tk.IntVar()

        #choose 1 or 2 component fits with radio buttons
        self.oneComponent = tk.Radiobutton(self,
                                           text = '1-component Fit',
                                           variable = self.v,
                                           value = 0)
        self.twoComponent = tk.Radiobutton(self,
                                           text = '2-component Fit',
                                           variable = self.v,
                                           value = 1)
        #button for plotting data
        self.dataButton = tk.Button(self,
                                    text = 'Plot Data',
                                    command = self.plotData
                                    )

        self.oneComponent.grid(row = 1, column = 0)
        self.twoComponent.grid(row = 2, column = 0)
        self.dataButton.grid(row = 1, column = 1)

        #fitting button, uses input from radio buttons
        self.fitButton = tk.Button(self,
                                   text = 'Fit That Bad Boy',
                                   command = self.curvefit,
                                   )
        self.fitButton.grid(row = 2, column = 1)

        #clear data from existing plots, necessary for fitting multiple sets of
        #data

        self.clearButton = tk.Button(self,
                                     text = 'Clear existing data',
                                     command = self.wipe)
        self.clearButton.grid(row = 3, column = 1)


    def wipe(self):
        self.ax.set_title('')

        if len(self.ax.lines) > 1:
            self.ax.get_legend().remove()

        while len(self.ax.lines) > 0:
            self.ax.lines[0].remove()

        self.canvas.draw()
        self.box.insert(tk.END, 'Data cleared from graph.' + '\n' + '\n', 'clear')


    def plotData(self):
        global direct
        global delta
        global Delta
        global gyro
        global xdata
        global ydata
        os.chdir(direct)
        parameters = ng.varian.read_procpar('procpar')



        #parses procpar file to produce gradient values (in G/cm),
        #big del and little del (in s)
        gzlvls = np.array(parameters['gzlvl1']['values']).astype(float)
        gradcal = float(parameters['gcal_']['values'][0])
        grad_values = gradcal * gzlvls

        Delta = float(parameters['del']['values'][0])
        delta = float(parameters['gt1']['values'][0])

        #list of gyromagnetic ratios (in Hz/T)

        gyro_list = {
            'H1':42577478.51818,
            'F19':40052000,
            'Li7':16546000,
            'Na23':11262000
            }

        gyro = gyro_list.get(parameters['tn']['values'][0])
        with open(r"integ_series.txt",  'r') as f:
            lines = f.readlines()
            int_messy = lines[-1]
            int_messy = int_messy.strip()
            integrals_almost = int_messy.split()
            integrals = np.array(integrals_almost[3:]).astype(float)

        #normalize
        normint = integrals / np.max(integrals)

        xdata = grad_values
        ydata = normint

        plt.clf()

        self.ax.plot(xdata,ydata,'p')
        self.ax.set_ylabel('Normalized Inteisity')
        self.ax.set_xlabel('Gradient (G/cm)')
        self.ax.set_title('%s Diffusion' % parameters['tn']['values'][0])
        self.canvas.draw()


    def curvefit(self):
        global delta
        global Delta
        global gyro
        global xdata
        global ydata
        global ax


        #self.plotData()

        #define the fitting functions, Stejskal-Tanner (1 component)
        def fitfunc1(x,I0,D):
            return I0 * np.exp(-D *
                               (gyro * 2 * np.pi * delta * x/100)**2 *
                               (Delta - (delta/3)))

        def fitfunc2(x,I1,D1,I2,D2):
            return I1 * np.exp(-D1 *
                               (gyro * 2 * np.pi * delta * x/100)**2 *
                               (Delta - (delta/3))) +   I2 * np.exp(-D2 *
                               (gyro * 2 * np.pi * delta * x/100)**2 *
                               (Delta - (delta/3)))

        funcparam = [fitfunc1, fitfunc2]
        initparam = [np.array([0.99,1e-11]), np.array([0.5,1e-11,0.5,1e-11])]


        popt, pcov = curve_fit(funcparam[self.v.get()], xdata, ydata, p0 = initparam[self.v.get()])

        #graph 1-component fit
        if self.v.get() == 0:
            #popt, pcov = curve_fit(funcparam[0], xdata, ydata, p0 = initparam[0])
            self.ax.plot(xdata, fitfunc1(xdata, *popt))

            self.ax.legend(labels = ('Experimental data',
                                     'Fit, I0 = %5.3f, D = %5.2e m^2/s' % tuple(popt)
                                     )
                           )
            stdDiff = np.sqrt(np.diag(pcov))[1]
            self.canvas.draw()
            self.box.insert(tk.END, 'Fit complete. Fitted Diffusion coefficient' +
                            ' is D = %5.2e m^2/s, Standard deviation' % popt[1] +
                            ' sigma = %5.1e m^2/s' % stdDiff +
                            '\n' + '\n',
                            'fit')
        else:
          # popt, pcov = curve_fit(funcparam[1], xdata, ydata,
           #                        p0 = initparam[1])

           self.ax.plot(xdata, fitfunc2(xdata,*popt))
           self.ax.legend(labels = ('Experimental Data',
                          'Fit, I1 = %5.3f, D1 = %5.2e m^2/s, I2 = %5.3f, D2 = %5.2e m^2/s' % tuple(popt)),
                          fontsize = 'xx-small'
                          )
           stdDiff2 = np.sqrt(np.diag(pcov))[[1,3]]
           self.canvas.draw()
           self.box.insert(tk.END,
                           'Fit Complete. Fitted Diffusion coefficients: ' +
                           'D1 = %5.2e m^2/s, D2 = %5.2e m^2s '
                           % tuple(popt[[1,3]]) + '\n' +
                           'Standard Deviations: sigma1 = %5.1e m^2/s, '
                           % stdDiff2[0] + 'sigma2 = %5.1e m^2/s' % stdDiff2[1]
                           + '\n' + '\n')


        self.box.insert(tk.END,
                        'Fit graph saved in ' + direct + '\n' + '\n',
                        'fit'
                        )





root = tk.Tk()
MainApplication(root)
root.mainloop()
