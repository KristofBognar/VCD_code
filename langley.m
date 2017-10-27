function langley_vec = langley(day, ampm, fd, sza, saa, dscd, rms, amf,...
    ideal_sza_range,filter_tag)
%langley_vec = langley(day, ampm, fd, sza, saa, dscd, rms, amf,ideal_sza_range)
% This function creates a langley plot and saves the output from the
% Langley plot as well as the interpolated air-mass-factors
% INPUT: 
%        day: measurement day
%        ampm: measurement twilight am = 0, pm = 1
%        fd: fractional day vector of measurements
%        sza: vector containing sza of measurements
%        saa: vector containing solar azimuth angles of measruements
%        rms: vector containing rms of measruements
%        amf: vector containing to amfs corresponding with measurements
%        (interpolated to correct SZA grid)
%        ideal_sza_range: [sza_min, sza_max]
% OUTPUT:
%        langley_vec: vector containing the following
%             1. day
%             2. ampm
%             3. fd min
%             4. fd max
%             5. sza min
%             6. sza max
%             7. az min
%             8. az max
%             9. R^2 value
%             10. slope 
%             11. slope err
%             12. y intercept (-RCD)
%             13. y intercept 1-sigma error
%             14. number of points used in Langley plot

% get string for ampm sorted out
if ampm == 0; ampm_str = 'am'; else ampm_str = 'pm'; end

% get SZA indices
[sza_range, sza_range_j]= get_SZA_indices(sza, ideal_sza_range);
min_j = min(sza_range_j);
max_j = max(sza_range_j);

% now find the indices that correspond with these values
% note that the slope and y-int variables contain two entries each:
% value and error
% calculate for alpha = 0.31731 (one-sigma)
[slope, y_int, R2, y_fit] =...
    linear_regression(amf(min_j:max_j), dscd(min_j:max_j), 0.31731, 0);

% make L-data
langley_vec = [day, ampm, fd(min_j), fd(max_j), ...
    sza_range, saa(min_j), saa(max_j), R2, slope, y_int, length(min_j:max_j)];

% make figure
fig_name = [num2str(day) ampm_str filter_tag];

subplot(3,1,1)
hold off
plot(amf(min_j:max_j), dscd(min_j:max_j), 'x')
hold on
plot(amf(min_j:max_j), y_fit,'color','r')
xlabel('amf')
ylabel('dscd')
legend('data', 'fit','location','NorthWest')
t1=['day: ' num2str(day) ' ' ampm_str ];
t2=['R^2: ' num2str(R2, 4)];
t3=['slope: '  num2str(slope(1), 2) '\pm' num2str(slope(2), 2)];
t4=['RCD: ' num2str(y_int(1), 2) '\pm' num2str(y_int(2),2)];
t5=['SZA: ' num2str(sza_range(1),3) '-' num2str(sza_range(2),3)];
title({t1; t2; t3; t4})
textbp([num2str(langley_vec(14)) ' measurements used']);

subplot(3,1,2)
hold off
plot(sza, dscd, '.')
hold on
plot([sza_range(1), sza_range(1)], [0, max(dscd)*1.3], 'm--')
plot([sza_range(2), sza_range(2)], [0, max(dscd)*1.3], 'm--')
xlabel('sza')
ylabel('dscd')

subplot(3,1,3)
hold off
plot(sza, rms, '.')
hold on
plot([sza_range(1), sza_range(1)], [0, max(rms)*1.3], 'm--')
plot([sza_range(2), sza_range(2)], [0, max(rms)*1.3], 'm--')
xlabel('sza')
ylabel('DOAS rms')


print_setting(1,1,fig_name);

end