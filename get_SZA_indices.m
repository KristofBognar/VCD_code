function [sza_range, sza_range_j]= get_SZA_indices(sza, ideal_sza_range)
% [sza_range, sza_range_j]= get_SZA_indices(sza, ideal_sza_range)
% This function finds the ideal sza range for Langley VCD retrievals.  It
% ideally chooses teh "Ideal_sza_range".  If not, it selects the closest
% interval to teh ideal sza range with the same width as the ideal sza
% range.
%  INPUT:  sza - vector of sza from morning or afternoon
%          ideal_sza_range - vector containing minimum sza and maximum sza
%  OUTPUT: sza_range - the resulting sza_range
%          sza_range_j - the indices of the resulting sza range

% determine what SZA we have available and select SZA range accordingly
sza_diff = (ideal_sza_range(2) - ideal_sza_range(1));
if min(sza) < ideal_sza_range(1)
    if max(sza) > ideal_sza_range(2)
        sza_range = ideal_sza_range;
    else
        sza_range =[max(sza)- sza_diff, max(sza)];
    end
else
    sza_range = [min(sza), min(sza) + sza_diff];
end
        
% now find the indices that correspond with these values
[tmp, min_j] = min(abs(sza - sza_range(1)));
[tmp, max_j] = min(abs(sza - sza_range(2)));
sza_range_j = [min_j, max_j];

%minimum and maximum sza
min_sza = sza(min_j);
max_sza = sza(max_j);
sza_range = [min_sza, max_sza];
end