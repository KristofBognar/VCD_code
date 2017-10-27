function daily_rcd = avg_daily_rcd_v2017(rcd_S, min_r2, min_nbr_pts, sza_min,save_fig,working_dir, filter_tag, min_SZA_range)
% Written: Cristen Adams
% December 2010
%
%ex:daily_rcd = avg_daily_rcd(rcd_S, 0, 8, 70);
% save in rcd_S
% rcd_S.mean.day = daily_rcd(:,1);rcd_S.mean.rcd = daily_rcd(:,2);rcd_S.mean.diff = daily_rcd(:,3);rcd_S.mean.err = daily_rcd(:,4);
% please do not save in daily_rcd! daily_rcd is not struct!
% dont's use:  daily_rcd.mean.day = daily_rcd(:,1);daily_rcd.mean.rcd = daily_rcd(:,2);daily_rcd.mean.diff = daily_rcd(:,3);daily_rcd.mean.err = daily_rcd(:,4);
%
% This function takes in an RCD vector and some filtering criteria for
% minimum R^2 and minimum number of points and minimum SZA and outputs a
% a daily averaged rcd
%
% INPUT
%      rcd_S: contains multiple entries
%             rcd_S.day
%             rcd_S.ampm
%             rcd_S.fd_min
%             rcd_S.fd_max
%             rcd_S.sza_min
%             rcd_S.sza_max
%             rcd_S.az_min
%             rcd_S.az_max
%             rcd_S.R2
%             rcd_S.m 
%             rcd_S.m_1sigma
%             rcd_S.y_int (-RCD)
%             rcd_S.y_int_1sigma
%             rcd_S.nbr
%    min_r2: float minimum R^2 value for Langley fit
%    min_nbr_pts: minimum number of points included in the langley plot
%    sza_min: minimum allowed SZA included in the fitting region (this is
%    used to filter out Langley plots which were produced right around
%    noon.
%
%   OUTPUT
%
%   Cristen:   daily_rcd: this vector contains the following columns
%             1. day
%             2. average rcd
%             3. difference in rcd values
%             4. combined error of rcd (from fitting error Y-int)
%  Xiaoyi:    rcd_S.mean: this vector contains the following columns
%             1. rcd_S.mean.day             day
%             2. rcd_S.mean.average         rcd
%             3. rcd_S.mean.diff            difference in rcd values
%             4. rcd_S.mean.err             combined error of rcd (from fitting error Y-int)

if nargin == 4 % do not save the RCD figure in default
   save_fig = 0;
   working_dir = pwd;
end
cd(working_dir);

daily_rcd = [];
for day = 1:366,
    SZA_range_in_langley = rcd_S.sza_max - rcd_S.sza_min;
    ind = find(rcd_S.day == day & rcd_S.R2 > min_r2 & ...
        rcd_S.nbr > min_nbr_pts & rcd_S.sza_min > sza_min & (SZA_range_in_langley >= min_SZA_range));
    if isempty(ind); continue; end
    if length(ind) == 1,
        daily_rcd = [daily_rcd;...
            day -rcd_S.y_int(ind)...
            NaN rcd_S.y_int_1sigma(ind)];
    elseif length(ind) == 2,
        daily_rcd = [daily_rcd; ...
            day -mean(rcd_S.y_int(ind)) abs(rcd_S.y_int(ind(1)) - rcd_S.y_int(ind(2))) ...
            sqrt(rcd_S.y_int_1sigma(ind(1))^2 +...
            rcd_S.y_int_1sigma(ind(1))^2 )];
    end
end

figure
hold on
errorbar(daily_rcd(:,1), daily_rcd(:,2), daily_rcd(:,3), 'r.');
errorbar(daily_rcd(:,1), daily_rcd(:,2), daily_rcd(:,4), 'b.');
legend('{\sigma}', 'combined Langley fit error');
xlabel('Day');
ylabel('RCD mol/cm^2');
fig_nm = strcat('RCDs', filter_tag);
print_setting(1/2,save_fig,fig_nm);