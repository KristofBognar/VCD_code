function [scd,scd_err_sys,scd_err_rand] = SAOZ_calc_scds_no2(dscd_S)
% To Cristen: 
% mean data is saved in rcd_S.mean struct, not in daily_rcd! 
% [scd,scd_err] = calc_scds_no2(dscd_S, daily_rcd)
% 

%Xiaoyi's note:
%ex:
% [dscd_S.scd,dscd_S.scd_err] = calc_scds_no2(dscd_S, rcd_S);
% step5: dscd_S.vcd = dscd_S.scd ./ dscd_S.amf;dscd_S.vcd_err = sqrt((dscd_S.scd_err ./dscd_S.scd).^2 + 0.02^2).* dscd_S.vcd;


%Written by: Cristen Adams, Dec 2010
% Modified by Kristof Bognar, March 2018
%
% INPUT
%       dscd_S: structure containing dscds
%       daily_rcd: structure containing rcds
%
% OUTPUT:
%       scd: vector containing SCDs for corresponding with the indices
%       of the dscd_S structure
%       scd_err: the total error error on the scd
%
% Kristof:  redistribute errors to match table in Hendrick et al., 2011
%           SCD errors from Cristen's thesis are split: Instrument and
%           Raman error go under random, cross-section errors go under
%           systematic
%           DSCD errors go under random
%           RCD errors go under systematic
%
% modified for SAOZ retrievals, use fix RCD values

year=dscd_S.year(1);
L = length(dscd_S.day);
scd = NaN * ones(L,1);

%% SAOZ RCD values for NO2
if year==2008 || year==2009
    saoz_rcd=1.7e16;
else
    saoz_rcd=1.1e16;
end

% use 2% as RCD error (Hendrick et al., 2011) -- RCD error not calculated
% by SAOZ team
% % rcd_err = 1e15 * ones(L,1);
rcd_err = (saoz_rcd/50) * ones(L,1);


%% calculate SCD
scd=dscd_S.mol_dscd + saoz_rcd;

%% calculate SCD random error (instr effects, DSCD errors)
scd_err_rand = sqrt( (3 * dscd_S.err ./ scd).^2 + (0.051).^2).* scd;

%% calculate SCD systematic error (c-s errors and RCD error)
scd_err_sys=sqrt( (0.094*scd).^2 + rcd_err.^2 );


end