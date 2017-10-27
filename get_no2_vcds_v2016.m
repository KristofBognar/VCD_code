function [dscd_S, rcd_S, avg_vcd] = get_no2_vcds_v2016(dscd_S, lambda, tag, sza_range, save_fig,working_dir,code_path)
% This function calculates no2 VCDs using a Langley plot analysis
%ex:[dscd_S, rcd_S, avg_vcd] = get_no2_vcds(dscd_S, 437,'u1');
%ex:[dscd_S, rcd_S, avg_vcd] = get_no2_vcds(dscd_S, 365,'p0');
% INPUT:
%      dscd_S: object described in the header of read_QDOAS.m
%      lambda: desired wavelength of AMFs (505 nm for ozone, 437 nm for
%      NO2)(in Fraser's thesis NO2 for 425nm) for NO2 in UV 365nm?
%      tag: string tag used to ID output AMF file.  e.g. 'no2_vis'
% OUTPUT:
%      dscd_S: same object as input, but with new entries for dscd_S.amf
%      rcd_s: object containing

if nargin == 3
     sza_range = [86,91];
     save_fig = 0;
     working_dir = pwd;
end
%[dscd_S,rcd_S]= get_all_rcds(dscd_S, 1, [86,91],0, 2, tag, lambda);
[dscd_S,rcd_S]= get_all_rcds_v2016(dscd_S, 1, sza_range,0, 2, tag, lambda,code_path);

%Average RCDs and create SCDs
tmp = avg_daily_rcd_v2016(rcd_S, 0, 8, 70,save_fig,working_dir);
rcd_S.mean.day = tmp(:,1);
rcd_S.mean.rcd = tmp(:,2);
rcd_S.mean.diff = tmp(:,3);
rcd_S.mean.err = tmp(:,4);

[dscd_S.scd, dscd_S.scd_err] = calc_scds_no2(dscd_S, rcd_S);

% Now average VCDs 
dscd_S.vcd = dscd_S.scd ./ dscd_S.amf;
dscd_S.vcd_err = sqrt( (dscd_S.scd_err ./dscd_S.scd).^2 + 0.05^2)...
    .* dscd_S.vcd;
%avg_vcd = avg_vcds(dscd_S, [86,90], 70, 8);
%avg_vcd = avg_vcds(dscd_S, [86,91], 70, 8);%modified by 2012 NDACC NO2 recommendation, the averaging VCD from 86-91
avg_vcd = avg_vcds(dscd_S, sza_range, 70, 8);%modified by 2012 NDACC NO2 recommendation, the averaging VCD from 86-91