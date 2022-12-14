# This program was written in its entirety by Jacob Ivanov, Undergraduate Research Assistant under Dr. Georges Pavlidis for the Nanoscale Imaging & Transport Laboratory at the University of Connecticut. Contact me at (860) 999-3105 or jacob.ivanov@uconn.edu for any issues.

# filepath must be either relative to current program or a global filepath
global_filepath = "/Users/jacobivanov/Library/CloudStorage/OneDrive-UniversityofConnecticut/NITL/Data Analysis/Data Files/Bivek Dec7 20 um Cooling Data.xlsx" 
data_filename = global_filepath.split('/')[-1]

# Below are the "toggles" for each time constant regression
N1 = True
N2 = True
N3 = False
N4 = False
N5 = False

# Below are the default bounds, though they can be changed below:
a_min, a_max = 0, 1000

# Standard Heating Bounds:
# b_min, b_max = -1000, 0
# Standard Cooling Bounds:
b_min, b_max = 0, 1000
tau_min, tau_max = 0, 1000

# Basic Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import random
from datetime import datetime

# Defines General Temperature Function
def general_temp_function(t, a, *params):
   # Seperates params arguments in respective lists, i.e.
   # general_temp_regress(t, a, b1, tau1, b2, tau2...)
   # b_list = [b1, b2, ...]
   # tau_list = [tau1, tau2, ...]
   b_list = []
   tau_list = []
   for i in range(0, len(params) // 2):
      pair = [params[i * 2], params[i * 2 + 1]]
      b_list.append(pair[0])
      tau_list.append(pair[1])

   # Calculates
   N = len(params) // 2
   sum_ = 0
   for term in range(0, N):
      sum_ += b_list[term] * np.exp(-t / tau_list[term])
   return a + sum_ 

# Defines Passing Functions
def pass1(t, a, b1, tau1):
   return general_temp_function(t, a, b1, tau1)
def pass2(t, a, b1, tau1, b2, tau2):
   return general_temp_function(t, a, b1, tau1, b2, tau2)
def pass3(t, a, b1, tau1, b2, tau2, b3, tau3):
   return general_temp_function(t, a, b1, tau1, b2, tau2, b3, tau3)
def pass4(t, a, b1, tau1, b2, tau2, b3, tau3, b4, tau4):
   return general_temp_function(t, a, b1, tau1, b2, tau2, b3, tau3, b4, tau4)
def pass5(t, a, b1, tau1, b2, tau2, b3, tau3, b4, tau4, b5, tau5):
   return general_temp_function(t, a, b1, tau1, b2, tau2, b3, tau3, b4, tau4, b5, tau5)

# Reads Data From Given File
df = pd.read_excel(global_filepath)

'''
EDIT THE BELOW DEPENDING ON THE DATA FILE
'''
# HARD CODED VERSION
# t = np.array([0, 20, 40, 60, 80, 100])
t = []
for i in range(0, df.shape[0]):
   # t.append(df.at[i, "Time"])
   t.append(df.at[i, "Time (Cooling)"])
t = np.array(t)
# print(t)

# HARD CODED VERSION
# T_data = np.array([0, 42.068, 58.668, 63.468, 65.668, 67.468])
T_data = []
for i in range(0, df.shape[0]):
   # T_data.append(df.at[i, "Temperature"])
   T_data.append(df.at[i, "Temperature (Cooling)"])
T_data = np.array(T_data)
# print(T_data)

# Bounds Generation
T_fit_N1bounds = (
   np.array([a_min, b_min, tau_min]), 
   np.array([a_max, b_max, tau_max]))

T_fit_N2bounds = (
   np.array([a_min, b_min, tau_min, b_min, tau_min]), 
   np.array([a_max, b_max, tau_max, b_max, tau_max]))

T_fit_N3bounds = (
   np.array([a_min, b_min, tau_min, b_min, tau_min, b_min, tau_min]), 
   np.array([a_max, b_max, tau_max, b_max, tau_max, b_max, tau_max]))

T_fit_N4bounds = (
   np.array([a_min, b_min, tau_min, b_min, tau_min, b_min, tau_min, b_min, tau_min]), 
   np.array([a_max, b_max, tau_max, b_max, tau_max, b_max, tau_max, b_max, tau_max]))

T_fit_N5bounds = (
   np.array([a_min, b_min, tau_min, b_min, tau_min, b_min, tau_min, b_min, tau_min, b_min, tau_min]), 
   np.array([a_max, b_max, tau_max, b_max, tau_max, b_max, tau_max, b_max, tau_max, b_max, tau_max]))

# Bounds and Curve Fit Calculation
if N1 == True:
   T_fit_N1params, T_fit_N1covars = curve_fit(pass1, t, T_data, bounds = T_fit_N1bounds)

   T_fit_N1uncers = []
   for coeff in range(0, 2 + 1):
      T_fit_N1uncers.append(np.sqrt(T_fit_N1covars[coeff][coeff]))
   T_fit_N1uncers = np.array(T_fit_N1uncers)


   T_fit_N1 = []
   for ti in t:
      T_fit_N1.append(general_temp_function(ti, 
      T_fit_N1params[0], 
      T_fit_N1params[1], 
      T_fit_N1params[2]))

if N2 == True:
   T_fit_N2params, T_fit_N2covars = curve_fit(pass2, t, T_data, bounds = T_fit_N2bounds)

   T_fit_N2uncers = []
   for coeff in range(0, 4 + 1):
      T_fit_N2uncers.append(np.sqrt(T_fit_N2covars[coeff][coeff]))
   T_fit_N2uncers = np.array(T_fit_N2uncers)

   T_fit_N2 = []
   for ti in t:
      T_fit_N2.append(general_temp_function(ti, 
      T_fit_N2params[0], 
      T_fit_N2params[1], 
      T_fit_N2params[2], 
      T_fit_N2params[3], 
      T_fit_N2params[4]))

if N3 == True:
   T_fit_N3params, T_fit_N3covars = curve_fit(pass3, t, T_data, bounds = T_fit_N3bounds)

   T_fit_N3uncers = []
   for coeff in range(0, 6 + 1):
      T_fit_N3uncers.append(np.sqrt(T_fit_N3covars[coeff][coeff]))
   T_fit_N3uncers = np.array(T_fit_N3uncers)

   T_fit_N3 = []
   for ti in t:
      T_fit_N3.append(general_temp_function(ti, 
      T_fit_N3params[0], 
      T_fit_N3params[1], 
      T_fit_N3params[2], 
      T_fit_N3params[3], 
      T_fit_N3params[4], 
      T_fit_N3params[5], 
      T_fit_N3params[6]))

if N4 == True:
   T_fit_N4params, T_fit_N4covars = curve_fit(pass4, t, T_data, bounds = T_fit_N4bounds)

   T_fit_N4uncers = []
   for coeff in range(0, 8 + 1):
      T_fit_N4uncers.append(np.sqrt(T_fit_N4covars[coeff][coeff]))
   T_fit_N4uncers = np.array(T_fit_N4uncers)

   T_fit_N4 = []
   for ti in t:
      T_fit_N4.append(general_temp_function(ti, 
      T_fit_N4params[0], 
      T_fit_N4params[1], 
      T_fit_N4params[2], 
      T_fit_N4params[3], 
      T_fit_N4params[4], 
      T_fit_N4params[5], 
      T_fit_N4params[6],
      T_fit_N4params[7], 
      T_fit_N4params[8]))

if N5 == True:
   T_fit_N5params, T_fit_N5covars = curve_fit(pass5, t, T_data, bounds = T_fit_N5bounds)

   T_fit_N5uncers = []
   for coeff in range(0, 10 + 1):
      T_fit_N5uncers.append(np.sqrt(T_fit_N5covars[coeff][coeff]))
   T_fit_N5uncers = np.array(T_fit_N5uncers)

   T_fit_N5 = []
   for ti in t:
      T_fit_N5.append(general_temp_function(ti, 
      T_fit_N5params[0], 
      T_fit_N5params[1], 
      T_fit_N5params[2], 
      T_fit_N5params[3], 
      T_fit_N5params[4], 
      T_fit_N5params[5], 
      T_fit_N5params[6],
      T_fit_N5params[7], 
      T_fit_N5params[8],
      T_fit_N5params[9], 
      T_fit_N5params[10]))

# Defines Coefficient of Determination (R^2) Function
def R2(data, fit):
   residuals = data - fit
   ss_res = np.sum(residuals ** 2)
   ss_tot = np.sum((data - np.mean(data)) ** 2)
   return 1 - (ss_res / ss_tot)

# Output .text file
# /Users/jacobivanov/Library/CloudStorage/OneDrive-UniversityofConnecticut/NITL/Data Analysis/Completed Regressions/data3.xlsx Temperature Regression Parameters.txt
with open("/Users/jacobivanov/Library/CloudStorage/OneDrive-UniversityofConnecticut/NITL/Data Analysis/Completed Regressions/{0} Temperature Regression Parameters.txt".format(data_filename), mode = 'w+') as output:
   output.write("The following text file was generated by General Temperature Regression.py program written by Jacob Ivanov, Undergraduate Researcher for the University of Connecticut Nanoscale Imaging & Transport Laboratory.\n")
   output.write("This file lists all parameters for each regression.\n")
   output.write("\n")
   output.write("Data File: {0}\n".format(global_filepath))
   output.write("Timestamp: {0}\n".format(datetime.now()))
   output.write("\n")

   if N1 == True:
      output.write("N = 1 Regression:\n")
      output.write("a    = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N1params[0], T_fit_N1uncers[0]))
      output.write("b1   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N1params[1], T_fit_N1uncers[1]))
      output.write("tau1 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N1params[2], T_fit_N1uncers[2]))
      output.write("R2   = {0:.6f}\n".format(R2(T_data, T_fit_N1)))
      output.write('\n')
   
   if N2 == True:
      output.write("N = 2 Regression:\n")
      output.write("a    = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N2params[0], T_fit_N2uncers[0]))
      output.write("b1   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N2params[1], T_fit_N2uncers[1]))
      output.write("tau1 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N2params[2], T_fit_N2uncers[2]))
      output.write("b2   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N2params[3], T_fit_N2uncers[3]))
      output.write("tau2 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N2params[4], T_fit_N2uncers[4]))
      output.write("R2   = {0:.6f}\n".format(R2(T_data, T_fit_N2)))
      output.write('\n')
   
   if N3 == True:
      output.write("N = 3 Regression:\n")
      output.write("a    = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[0], T_fit_N3uncers[0]))
      output.write("b1   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[1], T_fit_N3uncers[1]))
      output.write("tau1 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[2], T_fit_N3uncers[2]))
      output.write("b2   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[3], T_fit_N3uncers[3]))
      output.write("tau2 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[4], T_fit_N3uncers[4]))
      output.write("b3   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[5], T_fit_N3uncers[5]))
      output.write("tau3 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N3params[6], T_fit_N3uncers[6]))
      output.write("R2   = {0:.6f}\n".format(R2(T_data, T_fit_N3)))
      output.write('\n')

   if N4 == True:
      output.write("N = 4 Regression:\n")
      output.write("a    = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[0], T_fit_N4uncers[0]))
      output.write("b1   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[1], T_fit_N4uncers[1]))
      output.write("tau1 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[2], T_fit_N4uncers[2]))
      output.write("b2   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[3], T_fit_N4uncers[3]))
      output.write("tau2 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[4], T_fit_N4uncers[4]))
      output.write("b3   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[5], T_fit_N4uncers[5]))
      output.write("tau3 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[6], T_fit_N4uncers[6]))
      output.write("b4   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[7], T_fit_N4uncers[7]))
      output.write("tau4 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N4params[8], T_fit_N4uncers[8]))
      output.write("R2   = {0:.6f}\n".format(R2(T_data, T_fit_N4)))
      output.write('\n')

   if N5 == True:
      output.write("N = 5 Regression:\n")
      output.write("a    = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[0], T_fit_N5uncers[0]))
      output.write("b1   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[1], T_fit_N5uncers[1]))
      output.write("tau1 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[2], T_fit_N5uncers[2]))
      output.write("b2   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[3], T_fit_N5uncers[3]))
      output.write("tau2 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[4], T_fit_N5uncers[4]))
      output.write("b3   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[5], T_fit_N5uncers[5]))
      output.write("tau3 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[6], T_fit_N5uncers[6]))
      output.write("b4   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[7], T_fit_N5uncers[7]))
      output.write("tau4 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[8], T_fit_N5uncers[8]))
      output.write("b5   = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[9], T_fit_N5uncers[9]))
      output.write("tau5 = {0:+011.6f} ± {1:011.6f}\n".format(T_fit_N5params[10], T_fit_N5uncers[10]))
      output.write("R2   = {0:.6f}\n".format(R2(T_data, T_fit_N5)))
      output.write('\n')

# Output MatPlotLib Settings
if N2 == False:
   a, b = 1, 1
if N3 == False:
   a, b = 1, 2
if N4 == False:
   a, b = 1, 3
if N5 == False:
   a, b = 2, 2

# fig, ax = plt.subplots(2, 3, figsize = (17, 11), dpi = 80)
fig, ax = plt.subplots(2, 3, figsize = (17, 11), dpi = 80)

fig.suptitle("Temperature over Time Regression Comparison, Data File: {0}".format(data_filename))
plt.setp(ax[-1, :], xlabel= "Time (μs)")
plt.setp(ax[:, 0], ylabel= "Temperature (°C)")

# Raw Data
ax[0, 0].scatter(t, T_data, color = '#cccccc', label = 'Data')
ax[0, 0].set_title("Raw Data")
ax[0, 0].legend(loc = "upper right")

# N = 1
if N1 == True:
   ax[0, 1].scatter(t, T_data, color = '#cccccc')
   ax[0, 1].plot(t, T_fit_N1, color = '#3d9cbf', linestyle = 'dashed', label = "tau1 = {0:.3f},\nR2 = {1:.4f}".format(T_fit_N1params[2], R2(T_data, T_fit_N1)))
   ax[0, 1].set_title("N = 1")
   ax[0, 1].legend(loc = "upper right")

# N = 2
if N2 == True:
   ax[0, 2].scatter(t, T_data, color = '#cccccc')
   ax[0, 2].plot(t, T_fit_N2, color = '#6b9caa', linestyle = 'dashed', label = "tau1 = {0:.3f},\ntau2 = {1:.3f},\nR2 = {2:.4f}".format(T_fit_N2params[2], T_fit_N2params[4], R2(T_data, T_fit_N2)))
   ax[0, 2].set_title("N = 2")
   ax[0, 2].legend(loc = "upper right")

# N = 3
if N3 == True:
   ax[1, 0].scatter(t, T_data, color = '#cccccc')
   ax[1, 0].plot(t, T_fit_N3, color = '#507f80', linestyle = 'dashed', label = "tau1 = {0:.3f},\ntau2 = {1:.3f},\ntau3 = {2:.3f},\nR2 = {3:.4f}".format(T_fit_N3params[2], T_fit_N3params[4], T_fit_N3params[6], R2(T_data, T_fit_N3)))
   ax[1, 0].set_title("N = 3")
   ax[1, 0].legend(loc = "upper right")

# N = 4
if N4 == True:
   ax[1, 1].scatter(t, T_data, color = '#cccccc')
   ax[1, 1].plot(t, T_fit_N4, color = '#2f3942', linestyle = 'dashed', label = "tau1 = {0:.3f},\ntau2 = {1:.3f},\ntau3 = {2:.3f},\ntau4 = {3:.3f},\nR2 = {4:.4f}".format(T_fit_N4params[2], T_fit_N4params[4], T_fit_N4params[6], T_fit_N4params[8], R2(T_data, T_fit_N4)))
   ax[1, 1].set_title("N = 4")
   ax[1, 1].legend(loc = "upper right")

# N = 5
if N5 == True:
   ax[1, 2].scatter(t, T_data, color = '#cccccc')
   ax[1, 2].plot(t, T_fit_N5, color = '#153b42', linestyle = 'dashed', label = "tau1 = {0:.3f},\ntau2 = {1:.3f},\ntau3 = {2:.3f},\ntau4 = {3:.3f},\ntau5 = {4:.3f},\nR2 = {5:.4f}".format(T_fit_N5params[2], T_fit_N5params[4], T_fit_N5params[6], T_fit_N5params[8], T_fit_N5params[10], R2(T_data, T_fit_N5)))
   ax[1, 2].set_title("N = 5")
   ax[1, 2].legend(loc = "upper right")

plt.savefig("/Users/jacobivanov/Library/CloudStorage/OneDrive-UniversityofConnecticut/NITL/Data Analysis/Completed Regressions/{0} Regressions.png".format(data_filename))
plt.show()