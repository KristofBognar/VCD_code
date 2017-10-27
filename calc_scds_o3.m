function [scd,scd_err] = calc_scds_o3(dscd_S, daily_rcd)
% scd,scd_err] = calc_scds_o3(dscd_S, daily_rcd)
% Written by: Cristen Adams, Dec 2010
%
% INPUT
%       dscd_S: structure containing dscds
%       daily_rcd: structure containing rcds
%
% OUTPUT:
%       scd: vector containing SCDs for corresponding with the indices
%       of the dscd_S structure
%       scd_err: the total error error on the scd

L = length(dscd_S.day);
scd = NaN * ones(L,1);
rcd_err = NaN * ones(L,1);
for i = 1:length(daily_rcd.mean.day)
    ind = find(dscd_S.day == daily_rcd.mean.day(i));
    scd(ind) = dscd_S.mol_dscd(ind) + (daily_rcd.mean.rcd(i));
    rcd_err(ind) = sqrt( 1e18 .^2 + (daily_rcd.mean.diff(i)*0.5)^2);
end
try
    %scd_sys_err =sqrt( (0.02).^2 + (3*dscd_S.err ./ scd).^2 + ((dscd_S.x *
    %1e17)./dscd_S.mol_dscd).^2 ).*scd; % this 0.02 value was first found
    %in Cristen's original code. Cristen's thesis used value of 0.034, and
    %she confirmed we should use this one. Data generated for NDACC in 2017
    %start use 0.034 (before that, data, such as the ones for ACE archive was 0.02)
    
    scd_sys_err =sqrt( (0.034).^2 + (3*dscd_S.err ./ scd).^2 + ((dscd_S.x * 1e17)./dscd_S.mol_dscd).^2 ).*scd;
catch
    %scd_sys_err =sqrt( (0.02).^2 + (3*dscd_S.err ./ scd).^2).*scd;
    scd_sys_err =sqrt( (0.034).^2 + (3*dscd_S.err ./ scd).^2).*scd;
end
scd_err = sqrt(rcd_err.^2 + scd_sys_err.^2);
%scd_err = rcd_err;
