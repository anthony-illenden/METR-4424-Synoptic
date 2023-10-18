#noahbrauer, edited by Stacey Hitchcock & Michael Hosek 2023
############# Function to convert to isentropic coordinates is courtesy of MetPy; Modifications made by Noah Brauer


import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import xarray as xr
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, add_timestamp
from metpy.units import units

#List the files; YOU WILL NEED TO MODIFY THIS...

hgt_file = 'hgt.201102.nc'
temp_file = 'air.201102.nc'
q_file = 'shum.201102.nc'
u_file = 'uwnd.201102.nc'
v_file = 'vwnd.201102.nc'

#Read in each file 

nc_hgt = Dataset(hgt_file, 'r')
nc_temp = Dataset(temp_file, 'r')
nc_q = Dataset(q_file, 'r')
nc_u = Dataset(u_file, 'r')
nc_v = Dataset(v_file, 'r')


#Read in file attributes and extract times (Feb 19,2019 is what I'm pulling here)

lat = nc_hgt.variables['lat'][:]
lon = nc_hgt.variables['lon'][:] 

narr = {}
time = nc_hgt.variables['time'][:]
timeUnits = nc_hgt.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])

#Assign time index for our time of interest; YOU WILL HAVE TO CHANGE THIS LINE

time_index = np.where(narr['day']==2)[0]

#Now read in meteorological data for this day
level = nc_hgt.variables['level'][:] #In hPa
z = nc_hgt.variables['hgt'][time_index,:,:,:] #In meters
temp = nc_temp.variables['air'][time_index,:,:,:] #In Kelvin
q = nc_q.variables['shum'][time_index,:,:,:]
u = nc_u.variables['uwnd'][time_index,:,:,:]
v = nc_v.variables['vwnd'][time_index,:,:,:]

#Select your time of choice (03 UTC on 2/19 is the default now; YOU MAY NEED TO CHANGE THESE. Change the 1 to 2 if you want 06 UTC, etc.; times are in 3 hour increments)

temp = temp[0,:,:,:]
z = z[0,:,:,:]
q = q[0,:,:,:]
u = u[0,:,:,:]
v = v[0,:,:,:]

#Assign proper units to each variables
#Assign proper units to each variables
temp = units.Quantity(temp, "kelvin")
level = units.Quantity(level,"hectopascal")


#Define isentropic levels

isentlevs = [296.] * units.kelvin

#Now convert to isentropic coordinates

isentropic = mpcalc.isentropic_interpolation(isentlevs, level, temp, q, u, v, z, temperature_out=True)

isentprs, isenttmp, isentspech, isentu, isentv, isenthgt = isentropic

isentspec = isentspech*1000

def ms_to_knot(wind):
    knots = wind/0.514
    return knots

u_knots = ms_to_knot(isentu)
v_knots = ms_to_knot(isentv)

isentspech=units.Quantity(isentspech, "g/kg")

isentrh = 100 * mpcalc.relative_humidity_from_specific_humidity(isentprs, isenttmp, isentspech)

isenthgt=units.Quantity(isenthgt, "m")

msf = mpcalc.montgomery_streamfunction(isenthgt, isenttmp)

rh_nan = np.ones((1,277,349))*np.nan
u_knots_nan = np.ones((1,277,349))*np.nan
v_knots_nan = np.ones((1,277,349))*np.nan
isentprs_nan = np.ones((1,277,349))*np.nan
isentspec_nan = np.ones((1,277,349))*np.nan #in kg/kg
msf_nan = np.ones((1,277,349))*np.nan

for i in range(rh_nan.shape[1]):
    for j in range(rh_nan.shape[2]):
        
        if lon[i,j]>=0:
            rh_nan[:,i,j] = np.nan
            u_knots_nan[:,i,j] = np.nan
            v_knots_nan[:,i,j] = np.nan
            isentprs_nan[:,i,j] = np.nan
            isentspec_nan[:,i,j] = np.nan
            msf_nan[:,i,j] = np.nan
            
        else:
            rh_nan[:,i,j] = isentrh[:,i,j]
            u_knots_nan[:,i,j] = u_knots[:,i,j]
            v_knots_nan[:,i,j] = v_knots[:,i,j]
            isentprs_nan[:,i,j] = isentprs[:,i,j]
            isentspec_nan[:,i,j] = isentspec[:,i,j]
            msf_nan[:,i,j] = msf[:,i,j]



    
print(np.nanmax(isentspec_nan))

msf_2d = msf[0, :, :]

# Set up our projection
crs = ccrs.LambertConformal(central_longitude=-100.0, central_latitude=45.0)

# Coordinates to limit map area
bounds = [(-120., -75., 25., 50.)]
# Choose a level to plot, in this case 296 K
level = 0

#Can tweak figure size if you would like
fig = plt.figure(figsize=(15., 8.))

ax = fig.add_subplot(1, 1, 1, projection=crs)
ax.set_extent(*bounds, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.75)
ax.add_feature(cfeature.STATES, linewidth=0.5)

# Plot the surface
clevmsf = np.arange(200, 400, .5)
clevmsf = np.arange(0, 4000, 5)
cs = ax.contour(lon, lat, msf_nan[level, :, :]*10, levels=clevmsf,
                colors='k', linewidths=2.0, linestyles='solid', transform=ccrs.PlateCarree())
cs.clabel(fontsize=10, inline=1, inline_spacing=7, fmt='%i', rightside_up=True,
          use_clabeltext=True)

clevisent = np.arange(0, 1000, 25)
cs = ax.contour(lon, lat, isentprs_nan[level, :, :], clevisent,
                 colors='blue', linewidths=1.0, linestyles='solid', transform=ccrs.PlateCarree())
ax.clabel(cs, fontsize=10, inline=1, inline_spacing=7,
           fmt='%i', rightside_up=True, use_clabeltext=True)

#cf = ax.contourf(lon, lat, rh_nan[level, :, :], range(10, 106, 5), cmap=plt.cm.gist_earth_r, transform=ccrs.PlateCarree())
cf = ax.contourf(lon, lat, isentspec_nan[level, :, :], range(0,20, 1), cmap=plt.cm.gist_earth_r, transform=ccrs.PlateCarree())


#Add a colorbar
cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05,
                  extendrect='True')
cb.set_label('Specific Humidity (g/kg)', size='x-large')


#Plot wind barbs on the isentropic surface
ax.barbs(lon, lat, u_knots_nan[level, :, :], v_knots_nan[level, :, :], length=6,regrid_shape=20, transform=ccrs.PlateCarree())
plt.title('Montgomery Streamfunction, Specific Humidity, Wind 2/2/2011 0000 UTC', size = 20)
fig.tight_layout()

hgt_file = 'hgt.201601.nc'
temp_file = 'air.201601.nc'
q_file = 'shum.201601.nc'
u_file = 'uwnd.201601.nc'
v_file = 'vwnd.201601.nc'

#Read in each file 

nc_hgt = Dataset(hgt_file, 'r')
nc_temp = Dataset(temp_file, 'r')
nc_q = Dataset(q_file, 'r')
nc_u = Dataset(u_file, 'r')
nc_v = Dataset(v_file, 'r')


#Read in file attributes and extract times (Feb 19,2019 is what I'm pulling here)

lat = nc_hgt.variables['lat'][:]
lon = nc_hgt.variables['lon'][:] 

narr = {}
time = nc_hgt.variables['time'][:]
timeUnits = nc_hgt.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])

#Assign time index for our time of interest; YOU WILL HAVE TO CHANGE THIS LINE

time_index = np.where(narr['day']==22)[0]

#Now read in meteorological data for this day
level = nc_hgt.variables['level'][:] #In hPa
z = nc_hgt.variables['hgt'][time_index,:,:,:] #In meters
temp = nc_temp.variables['air'][time_index,:,:,:] #In Kelvin
q = nc_q.variables['shum'][time_index,:,:,:]
u = nc_u.variables['uwnd'][time_index,:,:,:]
v = nc_v.variables['vwnd'][time_index,:,:,:]

temp = temp[4,:,:,:]
z = z[4,:,:,:]
q = q[4,:,:,:]
u = u[4,:,:,:]
v = v[4,:,:,:]

#Assign proper units to each variables
#Assign proper units to each variables
temp = units.Quantity(temp, "kelvin")
level = units.Quantity(level,"hectopascal")


#Define isentropic levels

isentlevs = [296.] * units.kelvin

#Now convert to isentropic coordinates

isentropic = mpcalc.isentropic_interpolation(isentlevs, level, temp, q, u, v, z, temperature_out=True)

isentprs, isenttmp, isentspech, isentu, isentv, isenthgt = isentropic

isentspec = isentspech*1000

def ms_to_knot(wind):
    knots = wind/0.514
    return knots

u_knots = ms_to_knot(isentu)
v_knots = ms_to_knot(isentv)

isentspech=units.Quantity(isentspech, "g/kg")

isentrh = 100 * mpcalc.relative_humidity_from_specific_humidity(isentprs, isenttmp, isentspech)

isenthgt=units.Quantity(isenthgt, "m")

msf = mpcalc.montgomery_streamfunction(isenthgt, isenttmp)

rh_nan = np.ones((1,277,349))*np.nan
u_knots_nan = np.ones((1,277,349))*np.nan
v_knots_nan = np.ones((1,277,349))*np.nan
isentprs_nan = np.ones((1,277,349))*np.nan
isentspec_nan = np.ones((1,277,349))*np.nan #in kg/kg
msf_nan = np.ones((1,277,349))*np.nan

for i in range(rh_nan.shape[1]):
    for j in range(rh_nan.shape[2]):
        
        if lon[i,j]>=0:
            rh_nan[:,i,j] = np.nan
            u_knots_nan[:,i,j] = np.nan
            v_knots_nan[:,i,j] = np.nan
            isentprs_nan[:,i,j] = np.nan
            isentspec_nan[:,i,j] = np.nan
            msf_nan[:,i,j] = np.nan
            
        else:
            rh_nan[:,i,j] = isentrh[:,i,j]
            u_knots_nan[:,i,j] = u_knots[:,i,j]
            v_knots_nan[:,i,j] = v_knots[:,i,j]
            isentprs_nan[:,i,j] = isentprs[:,i,j]
            isentspec_nan[:,i,j] = isentspec[:,i,j]
            msf_nan[:,i,j] = msf[:,i,j]



    
print(np.nanmax(isentspec_nan))

msf_2d = msf[0, :, :]

# Set up our projection
crs = ccrs.LambertConformal(central_longitude=-100.0, central_latitude=45.0)

# Coordinates to limit map area
bounds = [(-120., -75., 25., 50.)]
# Choose a level to plot, in this case 296 K
level = 0

#Can tweak figure size if you would like
fig = plt.figure(figsize=(15., 8.))

ax = fig.add_subplot(1, 1, 1, projection=crs)
ax.set_extent(*bounds, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.75)
ax.add_feature(cfeature.STATES, linewidth=0.5)

# Plot the surface
clevmsf = np.arange(200, 400, .5)
clevmsf = np.arange(0, 4000, 5)
cs = ax.contour(lon, lat, msf_nan[level, :, :]*10, levels=clevmsf,
                colors='k', linewidths=2.0, linestyles='solid', transform=ccrs.PlateCarree())
cs.clabel(fontsize=10, inline=1, inline_spacing=7, fmt='%i', rightside_up=True,
          use_clabeltext=True)

clevisent = np.arange(0, 1000, 25)
cs = ax.contour(lon, lat, isentprs_nan[level, :, :], clevisent,
                 colors='blue', linewidths=1.0, linestyles='solid', transform=ccrs.PlateCarree())
ax.clabel(cs, fontsize=10, inline=1, inline_spacing=7,
           fmt='%i', rightside_up=True, use_clabeltext=True)

#cf = ax.contourf(lon, lat, rh_nan[level, :, :], range(10, 106, 5), cmap=plt.cm.gist_earth_r, transform=ccrs.PlateCarree())
cf = ax.contourf(lon, lat, isentspec_nan[level, :, :], range(0,20, 1), cmap=plt.cm.gist_earth_r, transform=ccrs.PlateCarree())


#Add a colorbar
cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05,
                  extendrect='True')
cb.set_label('Specific Humidity (g/kg)', size='x-large')


#Plot wind barbs on the isentropic surface
ax.barbs(lon, lat, u_knots_nan[level, :, :], v_knots_nan[level, :, :], length=6,regrid_shape=20, transform=ccrs.PlateCarree())
plt.title('Montgomery Streamfunction, Specific Humidity, Wind 1/22/2016 1200 UTC', size = 20)
fig.tight_layout()