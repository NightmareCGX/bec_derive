#!/home/sijie.pan/tools/anaconda3/bin/python
import os
import gc
import psutil
#from pympler import tracker, muppy, summary
import netCDF4
import numpy as np
import datetime as dt
import sys



fcstdir = "/work/wof/realtime/FCST/2019"
postdir = "/work/skinnerp/WoFS_summary_files/2019"
dates = ( "20190430", "20190501", "20190502", "20190503", "20190506", \
          "20190507", "20190508", "20190509", "20190510", "20190513", \
          "20190514", "20190515", "20190516", "20190517", "20190518", \
          "20190520", "20190521", "20190522", "20190523", "20190524", \
          "20190525", "20190526", "20190528", "20190529", "20190530" )

varnames = ["U", "V", "W", "PSFC", "T", "QVAPOR", "QCLOUD", \
            "QRAIN", "QSNOW", "QICE", "QGRAUP", "QHAIL"]
times = ("2000", "2030", "2100", "2130", "2200", "2230", "2300", "2330", \
         "0000", "0030", "0100", "0130", "0200", "0230", "0300")
fcstinterval = 30 #minutes
fcsthours = 3 #hours
fcstindex = 36 #the fcstindex-th(fcsthours*12) 5-min post output
fcstindex_trylimit = 31
#wrfout for 2018 events, wrfwof for 2019 events
if fcstdir[-4:] == "2018":
    fcstprefix = "wrfout"
elif fcstdir[-4:] == "2019":
    fcstprefix = "wrfwof"
ref_threshs = ["50", "40", "30", "20", "10", "0", "-999"]
#load memory tracking package
#tr = tracker.SummaryTracker()



class mosaic():
    

    def __init__(self, dir, time, timestr):
        self.__file = dir + "/news-e_SWT_" + time + ".nc"
        self.__pmmmosaic = np.empty(0)
        self.__posttime = time
        self.__fcsttime = timestr
        self.__TimeUpdate = False
        
        try:
            self.__index = fcstindex
            while self.__index >= fcstindex_trylimit:
                try:
                    mosaicIn = netCDF4.Dataset(self.__file, "r")
                    print("Opening pmm mosaic file {}".format(self.__file))
                    break
                except:
                    self.__index = self.__index - 1
                    self.__posttime = str(self.__index) + "_" + self.__posttime[3:17] + (dt.datetime.strptime(self.__posttime[17:21],"%H%M") - dt.timedelta(minutes = 5)).strftime("%H%M")
                    self.__file = dir + "/news-e_SWT_" + self.__posttime + ".nc"
                    self.__fcsttime = (dt.datetime.strptime(self.__fcsttime, "%Y-%m-%d_%H:%M:%S") - dt.timedelta(minutes = 5)).strftime("%Y-%m-%d_%H:%M:%S")
        except:
            pass
        else:
            if self.__index < fcstindex_trylimit:
                print("Forecast timestamp is {}! Exceed try limit {}!".format(fcstindex-self.__index+1, fcstindex-fcstindex_trylimit+1))
                sys.exit(1)
        
        if self.__posttime != time or self.__fcsttime != timestr:
            self.__TimeUpdate = True
       
        self.__pmmmosaic = np.array(mosaicIn.variables["comp_dz_pmm"][15:-15,15:-15])
        
        mosaicIn.close()
        del self.__file, mosaicIn
    
    
    def mosaicReturn(self):
        return self.__pmmmosaic


    def ifTimeUpdate(self):
        return self.__TimeUpdate


    def timeReturn(self):
        if hasattr(self,'_mosaic__posttime') and hasattr(self,'_mosaic__fcsttime'):
            return self.__fcsttime
        else:
            print("The mosaic file for 3-h forecast exists, no time update is necessary.")
            print("Please check the code to make sure it is running properly.")
            sys.exit(1)

         
    def memRelease(self):
        del self.__pmmmosaic
        gc.collect()




class wrfnc():
    
    
    __gravity = 9.81
    __p0 = 100000
    __RdCp = 287 / 1004
    
    
    def __init__(self, dir, timestr):
        self.__file = dir + "/" + fcstprefix + "_d01_" + timestr
        self.__coord = np.empty(0)
        self.__varValue = dict.fromkeys(varnames)
        
        try:
            wrfIn = netCDF4.Dataset(self.__file, "r")
            print("Opening ensemble output {}".format(self.__file))
        except:
            print("Ensemble file {} does not exist!".format(self.__file))
            sys.exit(1)
        
        #define vertical coord
        self.__coord = np.array(wrfIn.variables["P"][0,:,15:-15,15:-15] + wrfIn.variables["PB"][0,:,15:-15,15:-15])
        #tempPH = wrfIn.variables["PH"][0,:,:,:]
        #for nz in range(0,np.shape(tempPH)[0]):
        #    self.__coord = np.append(self.__coord, np.mean(wrfIn.variables["PHB"][0,nz,:,:] + tempPH[nz,:,:]))
        #del tempPH
        
        for var in varnames:
            if var == "PSFC":
                self.__varValue[var] = np.array(wrfIn.variables[var][0,15:-15,15:-15])
            elif var == "T":
#               self.__varValue[var] = (np.array(wrfIn.variables[var][0,:,15:-15,15:-15]) + 300) * (self.__coord / self.__p0) ** self.__RdCp
                self.__varValue[var] = np.array(wrfIn.variables[var][0,:,15:-15,15:-15])
            elif var == "U":
                self.__varValue[var] = np.array(wrfIn.variables[var][0,:,15:-15,15:-16])
            elif var == "V":
                self.__varValue[var] = np.array(wrfIn.variables[var][0,:,15:-16,15:-15])
            elif var == "W":
                self.__varValue[var] = np.array(wrfIn.variables[var][0,1:,15:-15,15:-15])
            else:
                self.__varValue[var] = np.array(wrfIn.variables[var][0,:,15:-15,15:-15])
                
        wrfIn.close()
        del self.__file, wrfIn, var


    def varReturn(self, var):
        if len(var) > 1 and var != "PSFC":
            return self.__varValue[var] * 1000
        else:
            return self.__varValue[var]


    def coordReturn(self):
        return self.__coord
    
    
    def memRelease(self):
        del self.__varValue, self.__coord
        gc.collect()




class case():
    
    
    def __init__(self, date):
        datetimes = []
        fcstdirs = []
        postdirs = []
        self.__press = dict.fromkeys(ref_threshs)
        self.__stds = dict.fromkeys(ref_threshs)
        self.__array = dict.fromkeys(ref_threshs)

        for ikey in ref_threshs:
            self.__stds[ikey] = dict.fromkeys(varnames)
        
        for index, time in enumerate(times):
            if index == 0:
                datetime = date + time
            else:
                datetime = (dt.datetime.strptime(datetime,"%Y%m%d%H%M") + dt.timedelta(minutes = fcstinterval)).strftime("%Y%m%d%H%M")
            
            fcst_dir = fcstdir + "/" + date + "/" + time
            post_dir = postdir + "/" + date + "/" + time
            datetimes.append(datetime)
            postdirs.append(post_dir)
            fcstdirs.append(fcst_dir)

        dir_tuple = tuple(zip(datetimes, postdirs, fcstdirs))
        del datetime, fcst_dir, post_dir, datetimes, postdirs, fcstdirs
        
        for index, (inittime, post_dir, fcst_dir) in enumerate(dir_tuple):
            self.__masks = dict.fromkeys(ref_threshs)
            self.__vars = dict.fromkeys(varnames)
            self.__coord = np.empty(0)
            print("Processing post-processed file for time: {}".format(inittime))
            
            #prepare time sting for read-in use
            fcsttime = (dt.datetime.strptime(inittime,"%Y%m%d%H%M") + dt.timedelta(hours = fcsthours)).strftime("%Y-%m-%d_%H:%M:%S")
            posttime = (dt.datetime.strptime(inittime,"%Y%m%d%H%M") + dt.timedelta(hours = fcsthours)).strftime("%H%M")
            posttime = str(fcstindex) + "_" + inittime[:8] + "_" + inittime[8:] + "_" + posttime
            
            #post file read-in
            postIn  = mosaic(post_dir, posttime, fcsttime)
            # check if time string needs to be updated
            if postIn.ifTimeUpdate():
                print("The original time string for WRF outputs is {}.".format(fcsttime))
                fcsttime = postIn.timeReturn()
                print('The new time string for WRF outputs is {}.'.format(fcsttime))
            for ikey in ref_threshs:
                # grid points with reflectivity > {thresh} dBZ are masked (i.e., true for >60dBZ, >50 and <60 ...)
                current_thresh  = float(ikey)
#               if ikey == "60":
#                   self.__masks[ikey] = np.ma.masked_greater(postIn.mosaicReturn(), current_thresh)
                if ikey == "50":
                    self.__masks[ikey] = np.ma.masked_greater_equal(postIn.mosaicReturn(), current_thresh)
                elif ikey == "-999":
                    current_thresh = 0.0
                    self.__masks[ikey] = np.ma.masked_less_equal(postIn.mosaicReturn(), current_thresh)
                elif current_thresh % 20 == 10:
                    self.__masks[ikey] = np.ma.masked_less_equal(abs(postIn.mosaicReturn() - current_thresh - 5.), 5.)
                else:
                    self.__masks[ikey] = np.ma.masked_less(abs(postIn.mosaicReturn() - current_thresh - 5.), 5.)

            # Release unused memory
            postIn.memRelease()
            
            #ensemble files read and process
            for i_ens in range(1, 19):
                memDir = fcst_dir + "/" + "ENS_MEM_" + str(i_ens)
                memIn = wrfnc(memDir, fcsttime)
                
                if i_ens == 1:
                    self.__coord = memIn.coordReturn()
                else:
                    self.__coord += memIn.coordReturn()
                    
                for varname in varnames:
                    if i_ens == 1:
                        self.__vars[varname] = memIn.varReturn(varname)
                    elif i_ens == 2:
                        self.__vars[varname] = np.stack((self.__vars[varname], memIn.varReturn(varname)), axis = 0)
                    else:
                        self.__vars[varname] = np.concatenate((self.__vars[varname], np.array([memIn.varReturn(varname)])), axis = 0)
                
                memIn.memRelease()

            del fcsttime, posttime, memDir, postIn, memIn
            
            #Calculating spread beween 6-h fcst members at a certain time
            # mask pressure based on post-process results
            self.__coord /= 18
            for ikey in ref_threshs:
                self.__devmasked = np.ma.masked_array(self.__coord, mask = np.broadcast_to(self.__masks[ikey].mask, self.__coord.shape))
                if index == 0:
                    self.__press[ikey] = list(self.__coord[self.__devmasked.mask])
                else:
                    self.__press[ikey] += list(self.__coord[self.__devmasked.mask])
            # mask dev based on post-process results
            for varname in varnames:
                self.__dev = np.zeros(np.shape(self.__vars[varname][0]))
                for i_ens in range(0, 18):
                    self.__dev += self.__vars[varname][i_ens]
                self.__vars[varname] = np.concatenate((self.__vars[varname], np.array([self.__dev / 18])), axis = 0)
                
                self.__dev = np.zeros(np.shape(self.__vars[varname][0]))
                for i_ens in range(0, 18):
                    self.__dev += np.square(self.__vars[varname][i_ens]-self.__vars[varname][18])
                #devide by 17 for unbiased estimation for 18 members
                self.__dev = np.sqrt(self.__dev / 17)
                
                #mask based on post-process results
                for ikey in ref_threshs:
                    self.__devmasked = np.ma.masked_array(self.__dev, mask = np.broadcast_to(self.__masks[ikey].mask, self.__dev.shape))
                    if index == 0:
                        self.__stds[ikey][varname] = list(self.__devmasked.data[self.__devmasked.mask])
                    else:
                        self.__stds[ikey][varname] += list(self.__devmasked.data[self.__devmasked.mask])
                    
            del i_ens, self.__vars, self.__coord, self.__masks, self.__devmasked, self.__dev
            gc.collect()
            print("Total CPU time for date ", inittime, ": ", psutil.cpu_times_percent())
            print("Total memory used for date ", inittime, ": ", psutil.Process(os.getpid()).memory_info().rss)
        
        #rearrangement and sorting
        for ikey in ref_threshs:
            self.__array[ikey] = [self.__press[ikey]]
            for varname in varnames:
                if varname != "PSFC":
                    self.__array[ikey].append(self.__stds[ikey][varname])
            self.__array[ikey] = np.array(self.__array[ikey])
            self.__array[ikey] = self.__array[ikey][:,self.__array[ikey][0,:].argsort()]
            self.__press[ikey] = self.__array[ikey][0]
            for index, varname in enumerate(varnames):
                if index == 0:
                    varind = index
                if varname != "PSFC":
                    varind = varind + 1          
                    self.__stds[ikey][varname] = self.__array[ikey][varind]
#       for index, ikey in enumerate(ref_threshs):
#          if index == 0:
#              self.__array_size = [len(self.__press[ikey])]
#          else:
#              self.__array_size.append(len(self.__press[ikey]))
#       print("Array size after sort: ")
#       printfmt = "{:^8s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}"
#       print(printfmt.format("varname","60 dBZ", "50 dBZ", "40 dBZ", "30 dBZ", "20 dBZ", "10 dBZ", "0 dBZ", "clear sky"))
#       printfmt = "{:^8s}{:>12d}{:>12d}{:>12d}{:>12d}{:>12d}{:>12d}{:>12d}{:>12d}"
#       print(printfmt.format("Pressure", self.__array_size[0], self.__array_size[1], self.__array_size[2], self.__array_size[3], self.__array_size[4], self.__array_size[5], self.__array_size[6], self.__array_size[7]))
#       for varname in varnames:
#           self.__array_size = []
#           for ikey in ref_threshs:
#               self.__array_size.append(len(self.__stds[ikey][varname]))
#           print(printfmt.format(varname, self.__array_size[0], self.__array_size[1], self.__array_size[2], self.__array_size[3], self.__array_size[4], self.__array_size[5], self.__array_size[6], self.__array_size[7]))
#       del self.__array_size
        del self.__array
        del index, inittime, post_dir, fcst_dir, dir_tuple, varname, varind
        gc.collect()
    
    
    def dataReturn(self, varname, thresh):
        if thresh not in ref_threshs:
            print("Error occurred in case.dataReturn.")
            print("{} is not a valid threshold.".format(thresh))
            sys.exit(1)
        if varname != "coord" and varname not in varnames:
            print("Error occurred in case.dataReturn.")
            print("{} is not a valid variable name.".format(varname))
            sys.exit(1)
        if varname == "coord":
            return self.__press[thresh]
        else:
            return self.__stds[thresh][varname]
           
        
    def memRelease(self):
        del self.__press, self.__stds
        gc.collect()




def main():
    

    for ind_date, date in enumerate(dates):
        caseStats = case(date)
        
        for thresh in ref_threshs:
            if thresh == "-999":
                statFile_name = "clearsky_"+ date + ".txt"
            else:
                statFile_name = "ref" + thresh + "dBZ_" + date + ".txt"
            with open(statFile_name, "w") as output:
                for n in range(len(caseStats.dataReturn("coord",thresh))):
                    print(caseStats.dataReturn("coord", thresh)[n], end = ' ', file = output)
                    for index, varname in enumerate(varnames):
                        if varname != "PSFC" and index != len(varnames) - 1:
                            print(caseStats.dataReturn(varname, thresh)[n], end = ' ', file = output)
                        elif varname != "PSFC":
                            print(caseStats.dataReturn(varname, thresh)[n], file = output)

        caseStats.memRelease()
        del caseStats
        gc.collect()
        print("Total CPU time after processing event ", date, ": ", psutil.cpu_times_percent())
        print("Total memory used after processing event ", date, ": ", psutil.Process(os.getpid()).memory_info().rss)




if __name__ == "__main__":
    main()
    #del tr
