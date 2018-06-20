function [dscd_S, rcd_S, avg_vcd] = get_ozone_vcds_v2018(dscd_S, sonde, tag,sza_range,lambda, save_fig,working_dir,code_path, filter_tag)
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
if isnan(dscd_S.o3)
    disp(['Ozone sonde data missing for this year, please double check sonde file']);
    prompt = 'Do you want contitue with 2-step rapid delivery calculation? Y/N [Y]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'Y';
        %  o3_flag: this flag tells us whether the ozone values input to the
        %  calculation code are SCDs or
        %  VCDs.  1 for O3 VCD in DU and 2 for O3 SCD in molec/cm2.
% %         o3_flag = 2; % if we do not have ozonesonde, or just not avalible for now, we still can run the VCD code, but LUT will use input of GBS dSCDs not ozonesonde VCDs
        
        %% 2-step method (Francois)
        % use fixed VCD for first run, then use first run VCDs for second run
        o3_flag = 1;
        
        % check if first run already happened
        if exist('../../vcd_1st_run.mat','file')
            
            disp(['2nd run: loaded temporary VCD file']);
            load('../../vcd_1st_run.mat')
            
            vcd_1st_run.mean_vcd=vcd_1st_run.mean_vcd/2.687e16;
            
            vcd_1st_run.mean_vcd(isnan(vcd_1st_run.mean_vcd))=300;
            
            vcd_1st_run_time=vcd_1st_run.year + vcd_1st_run.fd/n_days;
            dscd_S.o3=interp1(vcd_1st_run_time, vcd_1st_run.mean_vcd, measurement_time);
            
            delete('../../vcd_1st_run.mat')
            
        else
            disp(['1st run: using 300 DU as AMF LUT input']);
            
            % use fix value of 300 DU for first run
            dscd_S.o3=ones(size(measurement_time))*300;
            
        end
        
        
    end
    if strcmp(str,'N')
        return;
    end
else
    o3_flag = 1; % if we have ozonesonde measurements for a year, we will use ozonesonde VCD as LUT input by default
end
%[dscd_S,rcd_S]= get_all_rcds(dscd_S, 0, [86,91], 0, 1, [tag '_L1']);% Cristen
%[dscd_S,rcd_S]= get_all_rcds(dscd_S, 0, [86,91], 0, 2, [tag '_L1'],505);% Xiaoyi: change o3 flag and lamda if need
%[dscd_S,rcd_S]= get_all_rcds_v2016(dscd_S, 0, sza_range, 0, 2, [tag '_L2'],lambda,code_path, filter_tag);% Xiaoyi: change o3 flag and lamda if need
[dscd_S,rcd_S]= get_all_rcds_v2016(dscd_S, 0, sza_range, 0, o3_flag, [tag '_L2'],lambda,code_path, filter_tag);% Xiaoyi: change o3 for LUT input(2 for SCDs, 1 for VCDs)

%Average RCDs and create SCDs
%tmp = avg_daily_rcd_v2016(rcd_S, 0.9, 8, 70,save_fig,working_dir,filter_tag);
min_sza_range_in_langley = 2;% this is min SZA range for langley plot

if dscd_S.year(1)==2003 || dscd_S.year(1)==2004
    min_r2=0.8;
else
    min_r2=0.9;
end

tmp = avg_daily_rcd_v2017(rcd_S, min_r2, 8, 70,save_fig,working_dir, filter_tag, min_sza_range_in_langley);
rcd_S.mean.day = tmp(:,1);
rcd_S.mean.rcd = tmp(:,2);
rcd_S.mean.diff = tmp(:,3);
rcd_S.mean.err = tmp(:,4);
[dscd_S.scd, dscd_S.scd_err_sys, dscd_S.scd_err_rand] = calc_scds_o3_v2018(dscd_S, rcd_S);

% % Now calculate vcds
dscd_S.vcd = dscd_S.scd ./ dscd_S.amf;

% assign VCD errors
% replace 0.04 for AMF error with 0.0364 (Hendrick et al., 2011)
dscd_S.vcd_err_rand = sqrt( (dscd_S.scd_err_rand ./dscd_S.scd).^2 + 0.0364^2) .* dscd_S.vcd;
dscd_S.vcd_err_sys = (dscd_S.scd_err_sys ./dscd_S.scd) .* dscd_S.vcd;

% average twilight VCDs
avg_vcd = avg_vcds_v2018(dscd_S, sza_range, 70, 8);


