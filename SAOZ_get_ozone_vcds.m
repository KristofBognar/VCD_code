function [dscd_S, rcd_S, avg_vcd, dscd_S2, rcd_S2, avg_vcd2] = SAOZ_get_ozone_vcds(dscd_S, sonde, tag,sza_range,lambda, save_fig,working_dir,code_path, filter_tag)
% This function calculates o3 VCDs using a Langley plot analysis
%ex:[dscd_S, rcd_S, avg_vcd] = get_ozone_vcds(dscd_S, sonde_fd, 'u1');
% this function need read AMF! 
%Xiaoyi's note: sonde: 1. year 2. julian day 3. hours 4. ozone data (Integrated Ozone in
% DU)
%sonde_fd: 1. faraction time(including year) 2. ozone data
% output: struct: dscd_S, rcd_S, avg_vcd
%
% modified by Kristof Bognar, March 2018: calculate sys and rand errors
% separately, according to Hendrick et al., 2011 table

if nargin == 3 % in default setting we won't save figrues
   lambda = 505; % centre wavelength for ozone DOAS fitting
   sza_range = [86,91]; % SZA range will be used to perform langley fit
   save_fig = 0;
   working_dir = pwd;
end

% Run through first time to get output scds!
n_days = 365;
if rem(dscd_S.year(1),4) == 0, n_days = 366; end
sonde_time = sonde(:,1) +(sonde(:,2)+sonde(:,3)/24)/n_days; % time when we have the sonde measurement
sonde_ozone = sonde(:,4);% the sonde ozone VCD (in DU)
measurement_time = dscd_S.year + dscd_S.fd/n_days;%our measurement time
dscd_S.o3 = interp1(sonde_time, sonde_ozone, measurement_time);

fixRCD_tag='_fixRCD';
dailyRCD_tag='_dailyRCD';

% fix RCD using SCDs as LUT input (NOT dSCDs!!)
% SAOZ RCD values for O3
year=dscd_S.year(1);
if any(year==[2005:2007])
    saoz_rcd=4e19;
elseif any(year==[2008:2010])
    saoz_rcd=5e19;
elseif year==2011
    saoz_rcd=1.6e19;
elseif any(year==[2012:2017])
    saoz_rcd=4.4e19;
end
% calculate scd for AMF LUT input (2 for SCDs, 1 for VCDs)
dscd_S.scd_for_lut=dscd_S.mol_dscd+saoz_rcd;
% [dscd_S,rcd_S]= get_all_rcds_v2016(dscd_S, 0, sza_range, 0, 2, [tag '_L2'],lambda,code_path, fixRCD_tag);
% [dscd_S,rcd_S]= get_all_rcds_v2016(dscd_S, 0, sza_range, 0, 1, [tag '_L2'],lambda,code_path, fixRCD_tag);

% Copy of GBS method with SAOZ data
% [dscd_S2,rcd_S2]= get_all_rcds_v2016(dscd_S, 0, sza_range, 0, 1, [tag '_L2'],lambda,code_path, dailyRCD_tag);
dscd_S2=dscd_S;
rcd_S2=rcd_S;

%Average RCDs and create SCDs
min_sza_range_in_langley = 2;% this is min SZA range for langley plot
min_r2=0.9;

% tmp = avg_daily_rcd_v2017(rcd_S, min_r2, 8, 70,save_fig,working_dir, fixRCD_tag, min_sza_range_in_langley);
% rcd_S.mean.day = tmp(:,1);
% rcd_S.mean.rcd = tmp(:,2);
% rcd_S.mean.diff = tmp(:,3);
% rcd_S.mean.err = tmp(:,4);
[dscd_S.scd, dscd_S.scd_err_sys, dscd_S.scd_err_rand] = SAOZ_calc_scds_o3(dscd_S);

tmp = avg_daily_rcd_v2017(rcd_S2, min_r2, 8, 70,save_fig,working_dir, dailyRCD_tag, min_sza_range_in_langley);
rcd_S2.mean.day = tmp(:,1);
rcd_S2.mean.rcd = tmp(:,2);
rcd_S2.mean.diff = tmp(:,3);
rcd_S2.mean.err = tmp(:,4);
[dscd_S2.scd, dscd_S2.scd_err_sys, dscd_S2.scd_err_rand] = calc_scds_o3_v2018(dscd_S2, rcd_S2);


% % Now calculate vcds
dscd_S.vcd = dscd_S.scd ./ dscd_S.amf;
dscd_S2.vcd = dscd_S2.scd ./ dscd_S2.amf;

% assign VCD errors
% replace 0.04 for AMF error with 0.0364 (Hendrick et al., 2011)
dscd_S.vcd_err_rand = sqrt( (dscd_S.scd_err_rand ./dscd_S.scd).^2 + 0.0364^2) .* dscd_S.vcd;
dscd_S.vcd_err_sys = (dscd_S.scd_err_sys ./dscd_S.scd) .* dscd_S.vcd;

dscd_S2.vcd_err_rand = sqrt( (dscd_S2.scd_err_rand ./dscd_S2.scd).^2 + 0.0364^2) .* dscd_S2.vcd;
dscd_S2.vcd_err_sys = (dscd_S2.scd_err_sys ./dscd_S2.scd) .* dscd_S2.vcd;

% average twilight VCDs
avg_vcd = avg_vcds_v2018(dscd_S, sza_range, 70, 8);
avg_vcd2 = avg_vcds_v2018(dscd_S2, sza_range, 70, 8);


