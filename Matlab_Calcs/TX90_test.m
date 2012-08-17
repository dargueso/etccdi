syear=2040
eyear=2059

ncid=netcdf.open('/Users/daniel/Analyses/Sydney2km/Bias_correction/CCRC_DAM_2040-2059_TX_GAUSM.nc','NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'tasmax_bc');
Tmax=netcdf.getVar(ncid,varid);
ncid=netcdf.open('/Users/daniel/scratch/ETCCDI/fclimdex_RCM/thresholds.nc','NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'thresax90');
thx90=netcdf.getVar(ncid,varid);
thx90= permute(thx90,[3,1,2]);
thx90=repmat(thx90,[20,1,1]);

leapdays=localize_leapdays(syear,eyear);
ind_leapdays=find(leapdays==1)

n=1
for i=1:size(ind_leapdays)
   
thx90=insertrows(thx90,thx90(ind_leapdays(n)-1),ind_leapdays(n));
n=n+1
end

thx90=permute(thx90,[2,3,1]);


dx90=Tmax>=thx90;
year_d=years_day(syear,eyear);

for yr=syear:eyear
    
    dx90_y(:,:,yr-syear+1)=sum(dx90(:,:,year_d==yr),3);
end