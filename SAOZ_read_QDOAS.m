function [dscd_S, qdoas_filt, qdoas_raw, col] = SAOZ_read_QDOAS(data,trace_gas,fix_sza)
% [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS(file_nm, col, filt,2)
% version = 1 for DOAS data
% version = 2 for MAX-DOAS data, this function will be used as subfunction
% for read_QDOAS_maxdoas.m


% This function creates the daily data structure that the rest of the
% software uses to analyze VCDs.  

% for DOAS reading
% ex:[dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS('UE_2011_0660-70.asc',col_no2_p0, filt_good,1,0);
% ex: Ozone : [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS('UE_2011_600gr_fx.dat', col_o3_2, filt_good,1,0);

% INPUT:
%   file_nm:  string file-name of QDOAS output
%   col: structure containing column information listed below
    %   col.tot_nbr:  total number of columns in QDOAS output file
    %   col.fd: column for day of year
    %   col.sza: column for SZA
    %   col.saa: column for solar azimuth
    %   col.elev: column containing viewing elevation
    %   col.dscd: column for DSCD
    %   col.err: column for DSCD error
    %   col.rms: column for RMS
    %   col.x : column for the X DSCD if applicable (used for ozone error
    %   characterization!).
    %   col.shift : column containing the wavelength shift information
    %   col.stretch : column containing the wavelength shift information
    %   col.ref_sza: column containing the reference sza information
    %   col.year: column containing measurement year
%   filt: structure containing filtering information (if field is empty, do
%   not filter by this value)
    %   filt.elev: pointing elevation to be filtered by
    %   filt.rms_vec: vector containing two columns corresponding with SZA and RMS
    %   value.  This represents a threshold for RMS at the various SZA.
    %   filt.max_sza: maximum allowed SZA in the analysis (often set to 92 for
    %            AMF range)
    %   filt.shift: maximum allowed shift between reference and twilight
    %   spectrum.
%
% OUTPUT:
%    qdoas_filt: the full QDOAS data matrix corresponding to the
%               the filtered QDOAS data.
%    qdoas_raw: raw QDOAS data matrix
%    dscd_S: this is an object that contains the fields
%               which are relevant to the Langley analysis.
%               It contains columns the following columns all of which
%               should have the same length!
%                   dscd_S.day: day in local time
%                   dscd_S.fd :fractional day (UTC) -QDOAS things! NOT day
%                   of the year! 
%                   dscd_S.ampm: (am = 0, pm = 1)
%                   dscd_S.sza: SZA
%                   dscd_S.saa: solar azimuth
%                   dscd_S.rms : rms
%                   dscd_S.mol_dscd : focus molecule dscd
%                   dscd_S.err : focus molecule dscd fitting error
%                   dscd_S.x : fit of X cross-section
%                   dscd_S.shift: shift in spectrum versus reference (nm)
%                   dscd_S.stretch: strech in spectrum versus reference
%                   (nm)
%                   dscd_S.ref_sza: SZA of reference spectrum
%                   dscd_S.year: year of measurements


    save_fig = 0; % save the figures, 1 = yes, 0 = no
    working_dir = pwd; % this is the place you want save the figure

    if min(data.year)~=max(data.year), error('dSCD file contains multiple years'); end

    data_output = data ;
    N_total = size(data_output); col.tot_nbr = N_total(2);
    col.fd = find(strcmpi(data_output.Properties.VariableNames,'fd'));
    col.sza = find(strcmpi(data_output.Properties.VariableNames,'sza'));
    col.saa = find(strcmpi(data_output.Properties.VariableNames,'saa'));
    col.year = find(strcmpi(data_output.Properties.VariableNames,'year'));

    col.scans = find(strcmpi(data_output.Properties.VariableNames,'NbSomme'));
    col.tint = find(strcmpi(data_output.Properties.VariableNames,'Integration'));

    if trace_gas==1 % ozone dSCDs

        col.dscd = find(strcmpi(data_output.Properties.VariableNames,'O3'));
        col.err = find(strcmpi(data_output.Properties.VariableNames,'O3err'));

    elseif trace_gas==2 % no2 VIS 

        col.dscd = find(strcmpi(data_output.Properties.VariableNames,'NO2'));
        col.err = find(strcmpi(data_output.Properties.VariableNames,'NO2err'));

    end

    for i = 1:1:N_total(2)
        try
            qdoas_raw(:,i) = table2array(data_output(:,i));
        catch
            qdoas_raw(:,i) = repmat(-9999,N_total(1),1);
        end
    end
    

    % Now sort the data up by fractional day and ID whether there are
    % doubles of some values
    qdoas_raw = sortrows(qdoas_raw, col.fd);
    all_ind = 1:length(qdoas_raw(:,1));
    [a, unique_ind,b] = unique(qdoas_raw(:, col.fd));
    diff_ind = setdiff(all_ind, unique_ind);
    if ~isempty(diff_ind)
        disp('[WARNING]: File contains multiple entries taken at the same time.')
        disp('frac_day')
        disp('---------')
        disp(qdoas_raw(diff_ind, col.fd))
    end

    % filter by SZA only (max SZA is 92)
    if fix_sza
        qdoas_filt=qdoas_raw(qdoas_raw(:,col.sza)<91 & qdoas_raw(:,col.sza)>86,:);
        disp('Using fixed 86-91 SZA range')
    else
        qdoas_filt=qdoas_raw(qdoas_raw(:,col.sza)<92,:);
        disp('Using highest available 5 deg SZA range; max SZA=92')
    end
    
    % Make figures showing the results of the filtering
    figure
    subplot(2,1,1)
    hold on
    plot(qdoas_raw(:,col.fd), qdoas_raw(:, col.dscd), '.', 'color','r')
    plot(qdoas_filt(:,col.fd), qdoas_filt(:, col.dscd), 'b.')
    legend('Raw', 'Filt')
    xlabel('day')
    ylabel('DSCD (mol/cm^2)')
    subplot(2,1,2)
    hold on
    plot(qdoas_raw(:,col.fd), qdoas_raw(:, col.sza), '.', 'color','r')
    plot(qdoas_filt(:,col.fd), qdoas_filt(:, col.sza), 'b.')
    legend('Raw', 'Filt')
    xlabel('day')
    ylabel('SZA')
    cd(working_dir);
    print_setting(1/2,save_fig,['results_of_filtering']);

    % now sort data by day
    
    % sort data by azimuth and fractional day
    qdoas_filt = sortrows(qdoas_filt, col.saa);
    qdoas_filt = sortrows(qdoas_filt, col.fd);
    % find indices of data that should be on previous day
    prev_day_i = find(qdoas_filt(:, col.saa) > 0 ...
        & rem(qdoas_filt(:, col.fd),1) < 0.25);
    
    % find indices of data that should be pm
    pm_i = find(qdoas_filt(:, col.saa) > 0);
    
    % now create new data columns: 
    % local time day, fractional day, am/pm, sza, az, rms, dscd, err
    L = length(qdoas_filt(:,1));
    dscd_S.day = floor(qdoas_filt(:, col.fd));
    dscd_S.fd = qdoas_filt(:,col.fd); 
    dscd_S.ampm = zeros(L, 1);
    dscd_S.sza = qdoas_filt(:, col.sza);
    dscd_S.saa = qdoas_filt(:, col.saa);
    dscd_S.mol_dscd = qdoas_filt(:, col.dscd);
    dscd_S.err = qdoas_filt(:, col.err);
    dscd_S.rms = dscd_S.err*0;
%     dscd_S.tot_tint = qdoas_filt(:, col.scans).*qdoas_filt(:, col.tint);
    
    try
        dscd_S.year = qdoas_filt(:, col.year);
    catch
        disp('No column for year supplied')
    end
    dscd_S.ampm(pm_i) = 1;
    dscd_S.day(prev_day_i) = dscd_S.day(prev_day_i) - 1;
    
    % make figure
    figure
    subplot(2,1,1)
    hold on
    i1 = find(dscd_S.ampm == 0);
    plot(dscd_S.fd(i1), dscd_S.mol_dscd(i1),'.','color','b') 
    i2 = find(dscd_S.ampm == 1);
    plot(dscd_S.fd(i2), dscd_S.mol_dscd(i2),'.','color','r') 
    xlabel('Fractional day')
    ylabel('DSCD (mol/cm^2)')
    legend('am','pm')
    subplot(2,1,2)
    hold on
    plot(dscd_S.fd(i1), dscd_S.sza(i1),'.','color','b') 
    plot(dscd_S.fd(i2), dscd_S.sza(i2),'.','color','r')
    xlabel('Fractional day')
    ylabel('SZA')
    legend('am','pm')
    print_setting(1/2,save_fig,['filtered_DSCDs']);
    len_fin = length(dscd_S.sza);
    
end