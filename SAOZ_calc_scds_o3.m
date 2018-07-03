function [scd,scd_err_sys,scd_err_rand] = SAOZ_calc_scds_o3(dscd_S)
% scd,scd_err] = calc_scds_o3(dscd_S, daily_rcd)
% Written by: Cristen Adams, Dec 2010
% modified by: Xiaoyi Zhao, Nov 2017
% modified by: Kristof Bognar, March 2018
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
% Xiaoyi --- an uncertainty calculation log will be saved
% Kristof:  redistribute errors to match table in Hendrick et al., 2011
%           SCD errors from Cristen's thesis are split: Instrument and
%           Raman error go under random, cross-section errors go under
%           systematic
%           DSCD errors go under random
%           RCD errors go under systematic


%% calculate SCD and RCD error
year=dscd_S.year(1);

L = length(dscd_S.day);
scd = NaN * ones(L,1);


%% SAOZ RCD values for O3
if any(year==[2005:2007])
    saoz_rcd=4e19;
elseif any(year==[2008:2010])
    saoz_rcd=5e19;
elseif year==2011
    saoz_rcd=1.6e19;
elseif any(year==[2012:2017])
    saoz_rcd=4.4e19;
end

% use 10% as RCD error until values are provided by SAOZ team
rcd_err = (saoz_rcd/10) * ones(L,1);

scd = dscd_S.mol_dscd + saoz_rcd;


%% calculate SCD random error (instr effects, DSCD errors, X.xs)
scd_err_rand =sqrt( (0.014).^2 + (3*dscd_S.err ./ scd).^2).*scd;

%% calculate SCD systematic error (c-s errors and RCD error)
scd_err_sys=sqrt( (0.031*scd).^2 + rcd_err.^2 );
end
