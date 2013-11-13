%Heat Waves indices for Sydney 2km project
%It calculates:  
%   the number of days with Tmax > threshold
%   the number of heat waves (independent consecutive periods of at least 5 days with Tmax > threshold)
%   the number of heat waves days. Days within heat waves.
%   the same for Tmin using a different threshold.

%+==============================
% Daniel Argueso
% University of New South Wales
% 13 Aug 2012
% Last version: 13 Aug 2012
% email: d.argueso@unsw.edu.au
%==============================

addpath('~/Scripts/ccrc-scripts/Matlab/');





syear=1950;
eyear=2009;
nyears=eyear-syear+1;

m_clim=monthsyears(syear,eyear);
year_index= years_day(syear,eyear);
nmonths=monthsforave(syear,eyear);




%WRF
diriwrf='/home/z3393020/Analyses/NARCliM/NNRP/Bias_corrected/';
[diriwrf,'CCRC_NARCliM_R3_tasmin_bc.nc']
simid='NNRP'
ncid=netcdf.open([diriwrf,'CCRC_NARCliM_R3_tasmax_bc.nc'],'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'tasmax_bc');
tasmax_WRF=netcdf.getVar(ncid,varid);

ncid=netcdf.open([diriwrf,'CCRC_NARCliM_R3_tasmin_bc.nc'],'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'tasmin_bc');
tasmin_WRF=netcdf.getVar(ncid,varid);


filename='~/Analyses/share/geo_em_files/geo_em.d01.narclim.AUS_ONLY.nc';
ncid=netcdf.open(filename,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'XLAT_M');
lat=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'XLONG_M');
lon=netcdf.getVar(ncid,varid);;
d_sizes=size(tasmax_WRF);


%%AWAP
%simid='AWAP'
%diriwrfx='/home/z3393020/Analyses/AWAP/NARCliM/'
%diriwrfn='/home/z3393020/Analyses/AWAP/NARCliM/'
%ncid=netcdf.open([diriwrfx,'AWAP_tmax_1950-2009_tmax.nc'],'NC_NOWRITE');
%varid=netcdf.inqVarID(ncid,'tmax');
%tasmax_WRF=netcdf.getVar(ncid,varid);
%d_sizes=size(tasmax_WRF);
%
%print size
%varid=netcdf.inqVarID(ncid,'lon');
%lon=netcdf.getVar(ncid,varid);
%
%lon=repmat(lon,[1,d_sizes(2)]);
%
%varid=netcdf.inqVarID(ncid,'lat');
%lat=netcdf.getVar(ncid,varid);
%
%lat=repmat(lat,[1,d_sizes(1)])';
% 
%ncid=netcdf.open([diriwrfn,'AWAP_tmin_1950-2009_tmin.nc'],'NC_NOWRITE');
%varid=netcdf.inqVarID(ncid,'tmin');
%tasmin_WRF=netcdf.getVar(ncid,varid);


%##################




ncid=netcdf.create(['/home/z3393020/Analyses/NARCliM/heatwaves/BoM/CCRC_NARCliM_R3_heatwavesBoM_',simid,'.nc'],'NC_CLOBBER');
lat_dimID = netcdf.defDim(ncid,'y',d_sizes(2));
lon_dimID = netcdf.defDim(ncid,'x',d_sizes(1));
%time_dimID = netcdf.defDim(ncid,'time',max(nmonths));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'comment','File containing information of heat waves for both maximum and minimum temperature - as defined by BoM')
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'institution','Climate Change Research Centre (CCRC), University of New South Wales, Sydney, Australia')
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'author', 'NetCDF author: Daniel Argueso. d.argueso@unsw.edu.au')
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'date', date)

% variddx = netcdf.defVar(ncid,['TX90_days'],'float',[lon_dimID lat_dimID time_dimID])
% netcdf.putAtt(ncid,variddx,'long_name',['Number of days with Tmax > 90th percentile']);
% netcdf.putAtt(ncid,variddx,'units','days/year');
varidhwx= netcdf.defVar(ncid,['TX90_heatwaves'],'float',[lon_dimID lat_dimID])
netcdf.putAtt(ncid,varidhwx,'long_name',['Number of times Tmax > 90th percentile, for at least 3 consecutive days (',num2str(syear),'-',num2str(eyear),')']);
netcdf.putAtt(ncid,varidhwx,'units','');
varidhwdx= netcdf.defVar(ncid,['TX90_hw_days'],'float',[lon_dimID lat_dimID])
netcdf.putAtt(ncid,varidhwdx,'long_name',['Total number of heat-wave days (',num2str(syear),'-',num2str(eyear),')']);
netcdf.putAtt(ncid,varidhwx,'units','days');

% variddn = netcdf.defVar(ncid,['TN90_days'],'float',[lon_dimID lat_dimID time_dimID])
% netcdf.putAtt(ncid,variddn,'long_name',['Number of days with Tmin > 90th percentile']);
% netcdf.putAtt(ncid,variddn,'units','days/year');
varidhwn= netcdf.defVar(ncid,['TN90_heatwaves'],'float',[lon_dimID lat_dimID])
netcdf.putAtt(ncid,varidhwn,'long_name',['Number of times Tmin > 90th percentile, for at least 3 consecutive days (',num2str(syear),'-',num2str(eyear),')']);
netcdf.putAtt(ncid,varidhwn,'units','');
varidhwdn= netcdf.defVar(ncid,['TN90_hw_days'],'float',[lon_dimID lat_dimID])
netcdf.putAtt(ncid,varidhwdn,'long_name',['Total number of heat-wave days (',num2str(syear),'-',num2str(eyear),')']);
netcdf.putAtt(ncid,varidhwn,'units','days');


varidlat=netcdf.defVar(ncid,'lat','float',[lon_dimID lat_dimID]);
netcdf.putAtt(ncid,varidlat,'long_name','latitude');
netcdf.putAtt(ncid,varidlat,'units','degrees_north');
varidlon=netcdf.defVar(ncid,'lon','float',[lon_dimID lat_dimID]);
netcdf.putAtt(ncid,varidlon,'long_name','longitude');
netcdf.putAtt(ncid,varidlon,'units','degrees_east');
%varidtime=netcdf.defVar(ncid,'time','float',[time_dimID]);
%netcdf.putAtt(ncid,varidlon,'long_name','time');
%netcdf.putAtt(ncid,varidlon,'units',['days since',num2str(syear),'-01-01 00:00:00']);

%################ END OF USER MODIFICATION ########################




for mtc=1:12
p90x(:,:,mtc)=prctile(tasmax_WRF(:,:,m_clim==mtc),90,3);
p90n(:,:,mtc)=prctile(tasmin_WRF(:,:,m_clim==mtc),90,3);
end


    


for mth=1:max(nmonths)
    pmth=mod(mth,12);
    if pmth==0 
        pmth=12;
    end
    
    TX_days(:,:,nmonths==mth)=tasmax_WRF(:,:,nmonths==mth)>repmat(p90x(:,:,pmth),[1,1,sum(nmonths==mth)]);
    TN_days(:,:,nmonths==mth)=tasmin_WRF(:,:,nmonths==mth)>repmat(p90n(:,:,pmth),[1,1,sum(nmonths==mth)]);
end
    
    
    auxx=consecdays(TX_days);
    auxn=consecdays(TN_days);
    
%for mth=1:max(nmonths)
    TX_heatwaves(:,:)=sum(auxx(:,:,:)==3,3);
    TN_heatwaves(:,:)=sum(auxn(:,:,:)==3,3);
    TX_hwaves_days(:,:)=sum(auxx(:,:,:)>=3,3);
    TN_hwaves_days(:,:)=sum(auxn(:,:,:)>=3,3);
    
%end
    

   



netcdf.endDef(ncid)
%netcdf.putVar(ncid,variddx,single(TX_days))
netcdf.putVar(ncid,varidhwx,single(TX_heatwaves))
netcdf.putVar(ncid,varidhwdx,single(TX_hwaves_days))


%netcdf.putVar(ncid,variddn,single(TN_days))
netcdf.putVar(ncid,varidhwn,single(TN_heatwaves))
netcdf.putVar(ncid,varidhwdn,single(TN_hwaves_days))


netcdf.putVar(ncid,varidlon,lon)
netcdf.putVar(ncid,varidlat,lat)
%netcdf.putVar(ncid,varidtime,[1:max(nmonths)])

netcdf.close(ncid)



% varname=genvarname(['TX_',num2str(TX1)]);
%  eval([varname ' = zeros(d_sizes(1),d_sizes(2));'])
% varname=genvarname(['TX_',num2str(TX2)]);
%  eval([varname ' = zeros(d_sizes(1),d_sizes(2));'])
% varname=genvarname(['TN_',num2str(TN1)]);
%  eval([varname ' = zeros(d_sizes(1),d_sizes(2));'])
% varname=genvarname(['TN_',num2str(TN2)]);
%  eval([varname ' = zeros(d_sizes(1),d_sizes(2));'])



