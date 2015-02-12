#!/usr/bin/env python

""" var_info.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri Feb  6 11:02:11 AEDT 2015

"""


class VariablesInfo(object):
  
  RAW_SOURCE = {
                'DTR' :{'units': 'degC','long_name':'Diurnal temperature range'},
                'FD' :{'units': 'days','long_name':'Frosting days'},
                'ID' :{'units': 'days','long_name':'Icing days'},
                'SU' :{'units': 'days','long_name':'Summer days'},
                'TR' :{'units': 'days','long_name':'Tropical nights'},
                'TXx' :{'units': 'degC','long_name':'Maximum daily maximum temperature'},
                'TXn' :{'units': 'degC','long_name':'Minimum daily maximum temperature'},
                'TNx' :{'units': 'degC','long_name':'Maximum daily minimum temperature'},
                'TNn' :{'units': 'degC','long_name':'Minimum daily minimum temperature'},
                'R10mm' :{'units': 'days','long_name':'Days with rainfall larger than 10mm'},
                'R20mm' :{'units': 'days','long_name':'Days with rainfall larger than 20mm'},
                'Rnnmm' :{'units': 'days','long_name':''},
                'SDII' :{'units': 'mm/day','long_name':'Simple precipitaiton intensity index'},
				'R95p' :{'units': 'mm','long_name':'Accumulated preciptiation from events above the 95th percentile'},
				'R99p' :{'units': 'mm','long_name':'Accumulated preciptiation from events above the 99th percentile'},
				'PRCPtot' :{'units': 'mm','long_name':'Total accumulated preciptiation'},
				'Rx5day':{'units': 'mm','long_name':'5-day maximum accumulated preciptiation'},
				'Rx1day':{'units': 'mm','long_name':'Daily maximum preciptiation'},
				'TX10p':{'units': 'days','long_name':'No. days tmax below 10th percentile'},
				'TX50p':{'units': 'days','long_name':'No. days tmax above 50th percentile'},
				'TX90p':{'units': 'days','long_name':'No. days tmax above 90th percentile'},
				'TN10p':{'units': 'days','long_name':'No. days tmin below 10th percentile'},
				'TN50p':{'units': 'days','long_name':'No. days tmin above 50th percentile'},
				'TN90p':{'units': 'days','long_name':'No. days tmin above 90th percentile'},
				
                }
                
  def get_var_names(self):
    return self.RAW_SOURCE.keys()
  
  def get_units(self,varname):
    return self.RAW_SOURCE[varname]['units'][:]
    
  def get_longname(self,varname):
    return self.RAW_SOURCE[varname]['long_name'][:]