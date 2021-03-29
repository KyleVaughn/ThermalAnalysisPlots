import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np
import sys 
from scipy.signal import savgol_filter

# plots
ExportData_bool = False
plotTempvsMass_bool = False
plotTempvsdmdt_bool = True
plotTempvsdmdt_normalized_bool = False
plotinvKvsdmdt_bool = False
plotinvKvsdmdt_stepwise_normalized_bool = False
plotStepMassvsTime_bool = True
plotStepdmdtvsTime_bool = True

#Title
custom_title = False
custom_title_names = ['DrugName']

# Legend
custom_legend = False
custom_legend_names = ['TGA pan, step','DSC pan, step','TGA pan, ramp','PE pan, ramp','DSC pan, ramp']
#'Stepwise: 100-200 C','Stepwise: 180-330 C']
#, '26.3 mg','35.7 mg'
turn_legend_linefit_off = False

# Linear fit for enthalpy calc
linear = False # Perform the linear fit
linear_ranges = [(230,270)] # Temperature range in C for the linear fit

# Smoothing using Savgol filter
smooth = False # plot smoothed data
smooth_ramp_only_and_I_DEMAND_THE_SAME_COLORS = False # If only plotting smoothed ramp data, and u want the lines to be the original data color
smooth_window = 31 # points in smoothing window. must be an odd number
smooth_order = 1  # order of polynomial fit

# Scatter plot
scatter = True # Plot the full data as points
scatter_smoothed_plt4 = True # Plot the smoothed ramp data for the K^-1 vs -dmdt plot as point data
pt_size = 1 # point size

# Isotherm
# Assumes monotonic heating. Isotherm starts when |T - T_iso| < match tolerance
# Isotherm ends when |T - T_iso| > deviation tolerance
# Since temperature will fluctuate, and is output at discrete intervals, we need tolerances.
isotherm_match_tolerance = 2e-1 # if within this many degrees of isotherm, start isotherm
isotherm_deviation_tolerance = 3.5 # if isotherm deviates this many degrees, consider it stopped

# Fonts
fontname = "Arial"
fontsize = 12
fontsize_title = 14

################################################################################################
class data:
    def __init__(self, title, isRamp, time, temp, weight, iso_temps, iso_timeSteps):
        self.title = title
        self.isRamp = isRamp
        self.time = time
        self.temp = temp
        self.temp_K = temp + 273.15
        self.weight = weight
        self.iso_temps = np.array(iso_temps)
        if not isRamp:
            self.iso_temps_K = self.iso_temps + 273.15
        self.iso_timeSteps = np.array(iso_timeSteps)
        self.iso_dmdt = []
        self.iso_start_weight = []
        self.iso_stop_weight = []

        ref_weight = weight[0]
        self.weight_percent = 100.0 + 100.0*(weight-ref_weight)/ref_weight

        self.dmdt =  np.gradient(self.weight, self.time)
        self.dmdt_smooth = savgol_filter(self.dmdt, 31, 1)
        print(f"Reference weight for {title} is {ref_weight:.5f} mg")

    def printAll(self):
        print("Title: ", self.title)
        print("isRamp: ", self.isRamp)
        print("Time: ", self.time)
        print("Temp in C: ", self.temp)
        print("Temp in K: ", self.temp_K)
        print("Weight: ", self.weight)
        print("Weight %: ", self.weight_percent)
        print("Iso Temps: ", self.iso_temps)
        print("ISo Timesteps: ", self.iso_timeSteps)



###############################################################################################        
def CtoK(T):
    return T + 273.15

def KtoC(T):
    return T - 273.15

def inv2C(Tinv):
    return 1.0/Tinv - 273.15
 
def C2inv(T):
    return 1/(T + 273.15)


###############################################################################################
def plotTempvsMass(datalist):
    fig = plt.figure(1)
    for d in datalist:
        plt.plot(d.temp, d.weight_percent, label=d.title)
    plt.xlabel(r'Temperature ($\degree$C)', fontname=fontname, fontsize=fontsize)
    plt.ylabel('% Mass', fontname=fontname, fontsize=fontsize)
    plt.title("Percent Mass vs. Time", fontname=fontname, fontsize=fontsize_title)
    plt.tick_params(labelright=True, right=True, labeltop=True, top=True, which='both', direction='in')
    plt.minorticks_on()
    plt.legend()
    fig.tight_layout()


###############################################################################################
def plotTempvsdmdt(datalist):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    color_idx = 0

    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    secax = ax.secondary_xaxis('top', functions=(KtoC, CtoK))
    for d in datalist:
        if d.isRamp:
            if scatter:
                ax.scatter(d.temp_K[:-1], d.dmdt, s=pt_size, color=colors[color_idx], label=d.title)
        else:
            ax.plot(d.iso_temps_K, d.iso_dmdt, '-o', color=colors[color_idx], label=d.title)
            if scatter:
                ax.scatter(d.temp_K, d.dmdt, s=pt_size, color=colors[color_idx+1])
                color_idx = color_idx + 1
        color_idx = color_idx + 1
    if smooth_ramp_only_and_I_DEMAND_THE_SAME_COLORS:
        color_idx = 0
    for d in datalist:
        if d.isRamp:
            if smooth:
                ax.plot(d.temp_K[:-1], d.dmdt_smooth, color=colors[color_idx], label=d.title + ' smoothed')
                color_idx = color_idx + 1
    ax.set_xlabel('Temperature (K)', fontname=fontname, fontsize=fontsize)
    secax.set_xlabel(r'Temperature ($\degree$C)', fontname=fontname, fontsize=fontsize)
    ax.set_ylabel(r'$\frac{dm}{dt}$ (mg/min)', fontname=fontname, fontsize=fontsize)
    plt.title(r"$\frac{dm}{dt}$ vs. Temperature", fontname=fontname, fontsize=fontsize_title)
    ax.tick_params(labelright=True, right=True, which='both', direction='in')
    secax.tick_params(which='both', direction='in')
    ax.minorticks_on()
    secax.minorticks_on()
    ax.legend()
    fig.tight_layout()


###############################################################################################
def plotTempvsdmdt_normalized(datalist):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    color_idx = 0

    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    secax = ax.secondary_xaxis('top', functions=(KtoC, CtoK))
    for d in datalist:
        if d.isRamp:
            if scatter:
                ax.scatter(d.temp_K[:-1], 100.0*d.dmdt/d.weight[0], s=pt_size, color=colors[color_idx], label=d.title)
        else:
            ax.plot(d.iso_temps_K, 100.0*d.iso_dmdt/d.weight[0], '-o', color=colors[color_idx], label=d.title)
        color_idx = color_idx + 1
    if smooth_ramp_only_and_I_DEMAND_THE_SAME_COLORS:
        color_idx = 0
    for d in datalist:
        if d.isRamp:
            if smooth:
                ax.plot(d.temp_K[:-1], 100.0*d.dmdt_smooth/d.weight[0], color=colors[color_idx], label=d.title + ' smoothed')
                color_idx = color_idx + 1
    ax.set_xlabel('Temperature (K)', fontname=fontname, fontsize=fontsize)
    secax.set_xlabel(r'Temperature ($\degree$C)', fontname=fontname, fontsize=fontsize)
    ax.set_ylabel(r'$\frac{100}{m_0}\frac{dm}{dt}$ (%/min)', fontname=fontname, fontsize=fontsize)
    plt.title(r"Mass Loss Rate per mg Loaded vs. Temperature", fontname=fontname, fontsize=fontsize_title)
    ax.tick_params(labelright=True, right=True, which='both', direction='in')
    secax.tick_params(which='both', direction='in')
    ax.minorticks_on()
    secax.minorticks_on()
    ax.legend()
    fig.tight_layout()


########################################################################################################
def plotinvKvsdmdt(datalist):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    colors = colors + colors
    color_idx = 0

    fig = plt.figure(4)
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    secax = ax.secondary_xaxis('top', functions=(inv2C, C2inv))
    for d in datalist:
        if d.isRamp:
            if scatter:
                ax.plot(1.0/np.array(d.temp_K[:-1]), -d.dmdt, 'o', ms=pt_size, color=colors[color_idx], label=d.title)
        else:
            ax.plot(1.0/np.array(d.iso_temps_K), -d.iso_dmdt, '-o', color=colors[color_idx], label=d.title)
        color_idx = color_idx + 1
    if smooth_ramp_only_and_I_DEMAND_THE_SAME_COLORS:
        color_idx = 0
    for d in datalist:
        if d.isRamp:
            if smooth:
                if scatter_smoothed_plt4:
# Note that matplotlib has a problem with log scale scatter plots, hence using plot instead of scatter
# for figure 4 only.
                    ax.plot(1.0/np.array(d.temp_K[:-1]), -d.dmdt_smooth, 'o', ms=pt_size, color=colors[color_idx], label=d.title + ' smoothed')
                else:
                    ax.plot(1.0/np.array(d.temp_K[:-1]), -d.dmdt_smooth, color=colors[color_idx], label=d.title + ' smoothed')

            color_idx = color_idx + 1

    # Linear fit
    if linear:
        if len(linear_ranges) < 1:
            print("Linear range requires at least one tuple or temperatures: [(min temp0, max temp0), (min temp1, max temp1), etc.]")
            print("If the linear fit is not desired, set 'linear = False' in the options")
            sys.exit()
        for tup in linear_ranges:
            if len(tup) != 2:
                print("Linear range requires exactly two arguments per tuple: [(min temp, max temp)]")
                print("If the linear fit is not desired, set 'linear = False' in the options")
                sys.exit()
            if tup[0] >= tup[1]:
                print("Linear range requires [(min temp, max temp)], with min temp < max temp")
                sys.exit()
        if len(datalist) != len(linear_ranges):
            print(f"Length of datalist ({len(datalist)}) does not equal length of linear ranges ({len(linear_ranges)})")
            print("Linear range requires: [(min temp0, max temp0), (min temp1, max temp1), etc.]")
            print("If the linear fit is not desired, set 'linear = False' in the options")
            sys.exit()
        for i,d in enumerate(datalist):
            linear_range = linear_ranges[i]
            print(f'The linear range is {linear_range[0]}-{linear_range[1]}C')
            if d.isRamp:
                lin_idx = (linear_range[0] < d.temp) & (d.temp < linear_range[1])
                # Check that this range has 2 valid points 
                valid_pts= 0
                for tf in lin_idx:
                    if tf:
                        valid_pts = valid_pts + 1
                if valid_pts < 2:
                    print('Cannot find the necessary 2 points in the linear data range for a line fit.' 
                            + ' Using entire data set for line fit.')
                    lin_temp_inv0 = 1.0/np.array(d.temp_K)
                    lin_dmdt = -d.dmdt
                else:
                    lin_temp_inv0 = 1.0/np.array(d.temp_K[lin_idx])
                    lin_dmdt = -d.dmdt[lin_idx[:-1]]
            else:
                lin_idx = (linear_range[0] < d.iso_temps) & (d.iso_temps < linear_range[1])
                # Check that this range has 2 valid points 
                valid_pts= 0
                for tf in lin_idx:
                    if tf:
                        valid_pts = valid_pts + 1
                if valid_pts < 2:
                    print('Cannot find the necessary 2 points in the linear data range for a line fit.' \
                            + ' Using entire data set for line fit.')
                    lin_temp_inv0 = 1.0/d.iso_temps_K
                    lin_dmdt = -d.iso_dmdt
                else:
                    lin_temp_inv0 = 1.0/d.iso_temps_K[lin_idx]
                    lin_dmdt = -d.iso_dmdt[lin_idx]

            # Handle negative values
            pos_idx = (lin_dmdt > 0.0) 
            neg_ctr = 0
            for tf in pos_idx:
                if not tf:
                    neg_ctr = neg_ctr + 1
            if neg_ctr > 1:
                print(f"{neg_ctr} negative values encountered in -dm/dt." \
                        + "Cannot take log of these values, so I'm ignoring them")
            
            lin_temp_inv = lin_temp_inv0[pos_idx]
            lin_logdmdt = np.log(lin_dmdt[pos_idx])

            TT = np.linspace(lin_temp_inv[0], lin_temp_inv[-1], 100)
            results = np.polyfit(lin_temp_inv, lin_logdmdt, 1, full=True)
            m = results[0][0]
            b = results[0][1]
            if not results[1].size > 0:
                RSS = 0.0
            else:
                RSS = results[1][0]
            y = lin_logdmdt
            ybar = np.sum(y)/len(y)
            TSS = np.sum((y - ybar)**2)
            Rsquared = 1 - RSS/TSS
            print(f"For: {d.title}")
            if b > 0:
                print(f"Line:   y={m}x + {b}")
            else:
                print(f"Line:   y={m}x - {-b}")
            print(f"RSS:    {RSS}")
            print(f"R^2:    {Rsquared}")
            print(f"Enthalpy = {m*(-8.3145e-3):.3f} kJ/mol")
            line = np.poly1d(results[0])
            if b > 0:
                ax.plot(TT, np.exp(line(TT)), color=colors[color_idx],
                        label=r"$\ln\left(-\frac{dm}{dt}\right)$" \
                        + f"={m:.1f}" \
                        + r"$T^{-1}$+" \
                        + f"{b:.1f}," \
                        + f" $R^2$={Rsquared:.4f}")
            else:
                ax.plot(TT, np.exp(line(TT)), color=colors[color_idx],
                        label=r"$\ln\left(-\frac{dm}{dt}\right)$" \
                        + f"={m:.1f}" \
                        + r"$T^{-1}$-" \
                        + f"{-b:.1f}," \
                        + f" $R^2$={Rsquared:.4f}")
            color_idx = color_idx + 1





    ax.set_xlabel('1/T (K$^{-1}$)', fontname=fontname, fontsize=fontsize)
    secax.set_xlabel(r'Temperature ($\degree$C)', fontname=fontname, fontsize=fontsize)
    ax.set_ylabel(r'$-\frac{dm}{dt}$ (mg/min)', fontname=fontname, fontsize=fontsize)
    plt.title(r"$-\frac{dm}{dt}$ vs. 1/T", fontname=fontname, fontsize=fontsize_title)
    ax.tick_params(labelright=True, right=True, which='both', direction='in')
    secax.tick_params(which='both', direction='in')
    ax.minorticks_on()
    secax.minorticks_on()
    ax.legend()
    fig.tight_layout()






###############################################################################################
def ExportData(datalist):
    
    # Temp, inv_t, wt, wt%, dmdt, 
    for d in datalist:
        f = open(d.title + "_export.txt", "w")
        if d.isRamp:
            # write the header
            f.write("Title: " + d.title + "\n")
            f.write("Temp (K), Inverse Temp (1/K), Weight (mg), Weight %, dm/dt (mg/s) \n")
            # For each value in the range of dmdt (temp - 1), write the line
            for i in range(len(d.temp)-1):
                line = [d.temp_K[i], 1.0/d.temp_K[i], d.weight[i], d.weight_percent[i], d.dmdt[i]]
                str_line = str(line[0])
                for item in line[1:]:
                    str_line = str_line + ", " + str(item)
                f.write(str_line + "\n")
            # Write the last value for everything but dmdt, which is already exhausted
            i = len(d.temp) - 1
            line = [d.temp_K[i], 1.0/d.temp_K[i], d.weight[i], d.weight_percent[i]]
            str_line = str(line[0])
            for item in line[1:]:
                str_line = str_line + ", " + str(item)
            f.write(str_line + "\n")
        else: # it's not a ramp file
            print("Title: ", d.title)
            print("Temp (C), Temp (K), Inverse Temp (1/K), Start Weight (mg), Stop Weight (mg), dm/dt (mg/s)")
            for i in range(len(d.iso_temps_K) - 1):
                print(d.iso_temps[i], 
                    d.iso_temps_K[i], 
                    1.0/d.iso_temps_K[i], 
                    d.iso_start_weight[i], 
                    d.iso_stop_weight[i],
                    d.iso_dmdt[i]
                )
            i = len(d.iso_temps_K) - 1
            print(d.iso_temps[i], 
                d.iso_temps_K[i], 
                1.0/d.iso_temps_K[i], 
                d.iso_start_weight[i], 
                d.iso_stop_weight[i],
            )
        f.close()
        
        
        
###############################################################################################



def importData(custom_legend_names):
    datalist = []
    filenames = sys.argv[1:]

    # Did they give any files?
    if not filenames:
        print("No files given. Please specify files to process with 'python3 processTGA.py file1.txt file2.txt etc.'")
        sys.exit()

    # Are custom legend names being used?
    if custom_legend:
        if( len(custom_legend_names) != len(filenames) ):
                print("Custom legend names must correspond to same number of files")
                sys.exit()

    # Read data from each file
    for index, filename in enumerate(filenames):
        try:
            f = open(filename, 'r', encoding='utf-16')
            lines = f.readlines()
        except IOError:
            print("Could not open/read file:", filename)
            sys.exit()

        # Determine if ramp or step data
        isRamp = False
        for ctr, line in enumerate(lines):
            if 'ProcName' in line:
                words = line.split()
                if words[1] == 'Ramp':
                    isRamp = True
    
        # Import data
        if isRamp:
            iso_temps = []
            iso_timeSteps = []
            for ctr, line in enumerate(lines):
                if 'StartOfData' in line:
                   startOfDataLine = ctr
                elif 'Sample' in line:
                    words = line.split()
                    if words[0] == 'Sample':
                        title = words[1]
        else:
            iso_temps = []
            iso_timeSteps = []
            for ctr, line in enumerate(lines):
                if 'StartOfData' in line:
                   startOfDataLine = ctr
                elif 'Sample' in line:
                    words = line.split()
                    if words[0] == 'Sample':
                        title = words[1]
                elif 'OrgMethod' and 'Ramp' in line:
                    if 'Isothermal' in lines[ctr+1]:
                        # temp
                        words = line.split()
                        iso_temps.append(float(words[-2]))
                        # time step
                        words = lines[ctr+1].split()
                        iso_timeSteps.append(float(words[-2]))
    
        filedata = np.loadtxt(filename,
                          skiprows = startOfDataLine + 1,  
                          delimiter = '\t', 
                          encoding = 'utf-16',
                          )

        time = filedata[:,0]
        temp = filedata[:,1]
        weight = filedata[:,2]
        if custom_legend:
            title = custom_legend_names[index]

        # Initialize data objects
        datalist.append( data(title, isRamp, time, temp, weight, iso_temps, iso_timeSteps) )
    return datalist



###############################################################################################
def processStepData(datalist):
    # Assumed monotonic heating
    for d in datalist:
        if not d.isRamp:
            assert(len(d.iso_temps) == len(d.iso_timeSteps))
            startStopIndices = []
            for i in range(len(d.iso_temps)):
                iso_temp = d.iso_temps[i]
                iso_timeStep = d.iso_timeSteps[i]
                iso_start, iso_stop = getIsothermStartStopIndices(d.temp, d.time, iso_temp, iso_timeStep)
                startStopIndices.append((iso_start, iso_stop))

            iso_dmdt = np.zeros(len(d.iso_temps))
            temp_inv = np.zeros(len(d.iso_temps))
            for i in range(len(d.iso_temps)):
                mass_start = d.weight[startStopIndices[i][0]]
                mass_stop  = d.weight[startStopIndices[i][1]]
                d.iso_start_weight.append(mass_start)
                d.iso_stop_weight.append(mass_stop)
                time_start = d.time[startStopIndices[i][0]]
                time_stop  = d.time[startStopIndices[i][1]]
                iso_dmdt[i] = (mass_stop - mass_start)/(time_stop - time_start)
            d.iso_dmdt = iso_dmdt




###############################################################################################
def getIsothermStartStopIndices(temp_data, time_data, iso_temp, iso_timeStep):
    # Return first value within tolerance.
    # If no values in tolerance, return closest match
    error = abs(temp_data - iso_temp)
    start = len(temp_data)-1
    stop = len(temp_data)-1
    # start index
    for index, value in enumerate(error):
        if value < isotherm_match_tolerance:
            start = index
            break
    if start == len(temp_data)-1:
        start = error.argmin()
    
    # stop index
    time_stop = len(temp_data)-1
    temp_stop = len(temp_data)-1
    time_start = time_data[start]
    # temp stop
    for index, value in enumerate(error[start:]):
        if value > isotherm_deviation_tolerance:
            temp_stop = index + start
            break
    
    # time stop
    elapsed_time = time_data - time_start
    for index, value in enumerate(elapsed_time[start:]):
        if value >= iso_timeStep:
            time_stop = index + start
            break
    
    stop = min(time_stop, temp_stop)
    return start, stop




########################################################################################################
def plotinvKvsdmdt_stepwise_normalized(datalist):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    colors = colors + colors
    color_idx = 0

    fig = plt.figure(5)
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    secax = ax.secondary_xaxis('top', functions=(inv2C, C2inv))
    for d in datalist:
        if d.isRamp:
            if scatter:
                ax.plot(1.0/np.array(d.temp_K[:-1]), np.divide(-d.dmdt, d.weight[:-1]), 'o', ms=pt_size, color=colors[color_idx], label=d.title)
        else:
            ax.plot(1.0/np.array(d.iso_temps_K), np.divide(-d.iso_dmdt, d.weight[:-1]), '-o', color=colors[color_idx], label=d.title)
        color_idx = color_idx + 1
    if smooth_ramp_only_and_I_DEMAND_THE_SAME_COLORS:
        color_idx = 0
    for d in datalist:
        if d.isRamp:
            if smooth:
                if scatter_smoothed_plt4:
# Note that matplotlib has a problem with log scale scatter plots, hence using plot instead of scatter
# for figure 4 only.
                    ax.plot(1.0/np.array(d.temp_K[:-1]), np.divide(-d.dmdt_smooth, d.weight[:-1]), 'o', ms=pt_size, color=colors[color_idx], label=d.title + ' smoothed')
                else:
                    ax.plot(1.0/np.array(d.temp_K[:-1]), np.divide(-d.dmdt_smooth, d.weight[:-1]), color=colors[color_idx], label=d.title + ' smoothed')

            color_idx = color_idx + 1

    # Linear fit
    if linear:
        if len(linear_ranges) < 1:
            print("Linear range requires at least one tuple or temperatures: [(min temp0, max temp0), (min temp1, max temp1), etc.]")
            print("If the linear fit is not desired, set 'linear = False' in the options")
            sys.exit()
        for tup in linear_ranges:
            if len(tup) != 2:
                print("Linear range requires exactly two arguments per tuple: [(min temp, max temp)]")
                print("If the linear fit is not desired, set 'linear = False' in the options")
                sys.exit()
            if tup[0] >= tup[1]:
                print("Linear range requires [(min temp, max temp)], with min temp < max temp")
                sys.exit()
        if len(datalist) != len(linear_ranges):
            print(f"Length of datalist ({len(datalist)}) does not equal length of linear ranges ({len(linear_ranges)})")
            print("Linear range requires: [(min temp0, max temp0), (min temp1, max temp1), etc.]")
            print("If the linear fit is not desired, set 'linear = False' in the options")
            sys.exit()
        for i,d in enumerate(datalist):
            linear_range = linear_ranges[i]
            print(f'The linear range is {linear_range[0]}-{linear_range[1]}C')
            if d.isRamp:
                lin_idx = (linear_range[0] < d.temp) & (d.temp < linear_range[1])
                # Check that this range has 2 valid points 
                valid_pts= 0
                for tf in lin_idx:
                    if tf:
                        valid_pts = valid_pts + 1
                if valid_pts < 2:
                    print('Cannot find the necessary 2 points in the linear data range for a line fit.' 
                            + ' Using entire data set for line fit.')
                    lin_temp_inv0 = 1.0/np.array(d.temp_K)
                    lin_dmdt = np.divide(-d.dmdt, d.weight)
                else:
                    lin_temp_inv0 = 1.0/np.array(d.temp_K[lin_idx])
                    lin_dmdt = np.divide(-d.dmdt[lin_idx[:-1]], d.weight[lin_idx])
            else:
                lin_idx = (linear_range[0] < d.iso_temps) & (d.iso_temps < linear_range[1])
                # Check that this range has 2 valid points 
                valid_pts= 0
                for tf in lin_idx:
                    if tf:
                        valid_pts = valid_pts + 1
                if valid_pts < 2:
                    print('Cannot find the necessary 2 points in the linear data range for a line fit.' \
                            + ' Using entire data set for line fit.')
                    lin_temp_inv0 = 1.0/d.iso_temps_K
                    lin_dmdt = np.divide(-d.iso_dmdt, d.iso_start_weight)
                else:
                    lin_temp_inv0 = 1.0/d.iso_temps_K[lin_idx]
                    lin_dmdt = np.divide(-d.iso_dmdt[lin_idx], d.iso_start_weight[lin_idx])

            # Handle negative values
            pos_idx = (lin_dmdt > 0.0) 
            neg_ctr = 0
            for tf in pos_idx:
                if not tf:
                    neg_ctr = neg_ctr + 1
            if neg_ctr > 1:
                print(f"{neg_ctr} negative values encountered in -dm/dt." \
                        + "Cannot take log of these values, so I'm ignoring them")
            
            lin_temp_inv = lin_temp_inv0[pos_idx]
            lin_logdmdt = np.log(lin_dmdt[pos_idx])

            TT = np.linspace(lin_temp_inv[0], lin_temp_inv[-1], 100)
            results = np.polyfit(lin_temp_inv, lin_logdmdt, 1, full=True)
            m = results[0][0]
            b = results[0][1]
            if not results[1].size > 0:
                RSS = 0.0
            else:
                RSS = results[1][0]
            y = lin_logdmdt
            ybar = np.sum(y)/len(y)
            TSS = np.sum((y - ybar)**2)
            Rsquared = 1 - RSS/TSS
            print(f"For: {d.title}")
            if b > 0:
                print(f"Line:   y={m}x + {b}")
            else:
                print(f"Line:   y={m}x - {-b}")
            print(f"RSS:    {RSS}")
            print(f"R^2:    {Rsquared}")
            print(f"Enthalpy = {m*(-8.3145e-3):.3f} kJ/mol")
            line = np.poly1d(results[0])
            if b > 0:
                ax.plot(TT, np.exp(line(TT)), color=colors[color_idx],
                        label=r"$\ln\left(-\frac{dm}{dt}\right)$" \
                        + f"={m:.1f}" \
                        + r"$T^{-1}$+" \
                        + f"{b:.1f}," \
                        + f" $R^2$={Rsquared:.4f}")
            else:
                ax.plot(TT, np.exp(line(TT)), color=colors[color_idx],
                        label=r"$\ln\left(-\frac{dm}{dt}\right)$" \
                        + f"={m:.1f}" \
                        + r"$T^{-1}$-" \
                        + f"{-b:.1f}," \
                        + f" $R^2$={Rsquared:.4f}")
            color_idx = color_idx + 1





    ax.set_xlabel('1/T (K$^{-1}$)', fontname=fontname, fontsize=fontsize)
    secax.set_xlabel(r'Temperature ($\degree$C)', fontname=fontname, fontsize=fontsize)
    ax.set_ylabel(r'$-\frac{dm}{dt}$ (mg/min)', fontname=fontname, fontsize=fontsize)
    plt.title(r"$-\frac{dm}{dt}$ vs. 1/T", fontname=fontname, fontsize=fontsize_title)
    ax.tick_params(labelright=True, right=True, which='both', direction='in')
    secax.tick_params(which='both', direction='in')
    ax.minorticks_on()
    secax.minorticks_on()
    ax.legend()
    fig.tight_layout()


###############################################################################################
def plotStepMassvsTime(datalist):
    fig = plt.figure(6)
    for d in datalist:
        if d.isRamp:
            print("PlotStepMassvsTime is only valid for step files!")
            sys.exit()
        else:
            assert(len(d.iso_temps) == len(d.iso_timeSteps))
            startStopIndices = []
            for i in range(len(d.iso_temps)):
                iso_temp = d.iso_temps[i]
                iso_timeStep = d.iso_timeSteps[i]
                iso_start, iso_stop = getIsothermStartStopIndices(d.temp, d.time, iso_temp, iso_timeStep)
                startStopIndices.append((iso_start, iso_stop))

            iso_mass = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            iso_time = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            for i in range(len(d.iso_temps)):
                iso_mass = d.weight[startStopIndices[i][0] : startStopIndices[i][1] + 1]
                iso_time = d.time[startStopIndices[i][0] : startStopIndices[i][1] + 1]
                plt.plot(iso_time, iso_mass, label=r"" + d.title + " " + str(d.iso_temps[i]) + "$\degree$C") 

    plt.xlabel(r'Time (min)', fontname=fontname, fontsize=fontsize)
    plt.ylabel('Mass (g)', fontname=fontname, fontsize=fontsize)
    plt.title("Mass vs. Time", fontname=fontname, fontsize=fontsize_title)
    plt.tick_params(labelright=True, right=True, labeltop=True, top=True, which='both', direction='in')
    plt.minorticks_on()
    plt.legend()
    fig.tight_layout()


    fig = plt.figure(7)
    for j, d in enumerate(datalist):
        if d.isRamp:
            print("PlotStepMassvsTime is only valid for step files!")
            sys.exit()
        else:
            assert(len(d.iso_temps) == len(d.iso_timeSteps))
            startStopIndices = []
            for i in range(len(d.iso_temps)):
                iso_temp = d.iso_temps[i]
                iso_timeStep = d.iso_timeSteps[i]
                iso_start, iso_stop = getIsothermStartStopIndices(d.temp, d.time, iso_temp, iso_timeStep)
                startStopIndices.append((iso_start, iso_stop))
                                                                                                              
            iso_mass = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            iso_time = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            for i in range(len(d.iso_temps)):
                if linear_ranges[j][0] <= d.iso_temps[i] and d.iso_temps[i] <= linear_ranges[j][1]:
                    iso_mass = d.weight[startStopIndices[i][0] : startStopIndices[i][1] + 1]
                    iso_time = d.time[startStopIndices[i][0] : startStopIndices[i][1] + 1] - d.time[startStopIndices[i][0]]
                    plt.plot(iso_time, iso_mass, label=r"" + d.title + " " + str(d.iso_temps[i]) + "$\degree$C") 
                                                                                                              
    plt.xlabel(r'Time (min)', fontname=fontname, fontsize=fontsize)
    plt.ylabel('Mass (g)', fontname=fontname, fontsize=fontsize)
    plt.title("Mass vs. Time", fontname=fontname, fontsize=fontsize_title)
    plt.tick_params(labelright=True, right=True, labeltop=True, top=True, which='both', direction='in')
    plt.minorticks_on()
    plt.legend()
    fig.tight_layout()


###############################################################################################
def plotStepdmdtvsTime(datalist):
    fig = plt.figure(8)
    for j, d in enumerate(datalist):
        if d.isRamp:
            print("PlotStepMassvsTime is only valid for step files!")
            sys.exit()
        else:
            assert(len(d.iso_temps) == len(d.iso_timeSteps))
            startStopIndices = []
            for i in range(len(d.iso_temps)):
                iso_temp = d.iso_temps[i]
                iso_timeStep = d.iso_timeSteps[i]
                iso_start, iso_stop = getIsothermStartStopIndices(d.temp, d.time, iso_temp, iso_timeStep)
                startStopIndices.append((iso_start, iso_stop))
                                                                                                              
            iso_dmdt = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            iso_time = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            for i in range(len(d.iso_temps)):
                if linear_ranges[j][0] <= d.iso_temps[i] and d.iso_temps[i] <= linear_ranges[j][1]:
                    iso_dmdt = d.dmdt[startStopIndices[i][0] : startStopIndices[i][1] + 1]
                    iso_time = d.time[startStopIndices[i][0] : startStopIndices[i][1] + 1] - d.time[startStopIndices[i][0]]
                    plt.plot(iso_time, iso_dmdt, label=r"" + d.title + " " + str(d.iso_temps[i]) + "$\degree$C") 
            for i in range(len(d.iso_temps)):
                if linear_ranges[j][0] <= d.iso_temps[i] and d.iso_temps[i] <= linear_ranges[j][1]:
                    iso_dmdt = d.dmdt_smooth[startStopIndices[i][0] : startStopIndices[i][1] + 1]
                    iso_time = d.time[startStopIndices[i][0] : startStopIndices[i][1] + 1] - d.time[startStopIndices[i][0]]
                    plt.plot(iso_time, iso_dmdt, label=r"" + d.title + " " + str(d.iso_temps[i]) + "$\degree$C smooth") 
                                                                                                              
    plt.xlabel(r'Time (min)', fontname=fontname, fontsize=fontsize)
    plt.ylabel(r'$\frac{dm}{dt}$ (mg/min)', fontname=fontname, fontsize=fontsize)
    plt.title(r"$\frac{dm}{dt}$ vs. Time", fontname=fontname, fontsize=fontsize_title)
    plt.tick_params(labelright=True, right=True, labeltop=True, top=True, which='both', direction='in')
    plt.minorticks_on()
    plt.legend()
    fig.tight_layout()


    fig = plt.figure(9)
    for j, d in enumerate(datalist):
        if d.isRamp:
            print("PlotStepMassvsTime is only valid for step files!")
            sys.exit()
        else:
            assert(len(d.iso_temps) == len(d.iso_timeSteps))
            startStopIndices = []
            for i in range(len(d.iso_temps)):
                iso_temp = d.iso_temps[i]
                iso_timeStep = d.iso_timeSteps[i]
                iso_start, iso_stop = getIsothermStartStopIndices(d.temp, d.time, iso_temp, iso_timeStep)
                startStopIndices.append((iso_start, iso_stop))
                                                                                                              
            iso_dmdt = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            iso_time = np.zeros(startStopIndices[i][1] - startStopIndices[i][0] + 1)
            for i in range(len(d.iso_temps)):
                if linear_ranges[j][0] <= d.iso_temps[i] and d.iso_temps[i] <= linear_ranges[j][1]:
                    iso_dmdt = d.dmdt_smooth[startStopIndices[i][0] : startStopIndices[i][1] + 1]
                    iso_time = d.time[startStopIndices[i][0] : startStopIndices[i][1] + 1] - d.time[startStopIndices[i][0]]
                    plt.plot(iso_time, np.gradient(iso_dmdt, iso_time), 
                            label=r"" + d.title + " " + str(d.iso_temps[i]) + "$\degree$C") 
                                                                                                              
    plt.xlabel(r'Time (min)', fontname=fontname, fontsize=fontsize)
    plt.ylabel(r'$\frac{d^2m}{dt^2}$ (mg/min$^2$)', fontname=fontname, fontsize=fontsize)
    plt.title(r"$\frac{d^2m}{dt^2}$ vs. Time", fontname=fontname, fontsize=fontsize_title)
    plt.tick_params(labelright=True, right=True, labeltop=True, top=True, which='both', direction='in')
    plt.minorticks_on()
    plt.legend()
    fig.tight_layout()

###############################################################################################
if __name__ == "__main__":
    datalist = importData(custom_legend_names)
    processStepData(datalist)
    if ExportData_bool:
        ExportData(datalist)
    if plotTempvsMass_bool:
        plotTempvsMass(datalist)
    if plotTempvsdmdt_bool:
        plotTempvsdmdt(datalist)
    if plotTempvsdmdt_normalized_bool:
        plotTempvsdmdt_normalized(datalist)
    if plotinvKvsdmdt_bool:
        plotinvKvsdmdt(datalist)
    if plotinvKvsdmdt_stepwise_normalized_bool:
        plotinvKvsdmdt_stepwise_normalized(datalist)
    if plotStepMassvsTime_bool:
        plotStepMassvsTime(datalist)
    if plotStepdmdtvsTime_bool:
        plotStepdmdtvsTime(datalist)
    plt.show()
