function [rcd_vec, amf] = ... 
    rcd_ndacc_no2_v2016(year, day, ampm, dscd_S, sza_range, tag, lambda,code_path,NO2_AMF_version, filter_tag)

%  This function calculates the RCDs for ozone for a given twilight, SZA range, 
%  and dscd_vecset.
%  *** NOTE directory is hard-programed.  Change the amf_dir variable if
%  would like to set up in another directory structure.
%
% Xiaoyi's note: before using this M function, the working folder should
% have /AMF sub folder, normally runing this function with Current Directory as QDOAS/Output folder
%
% INPUT:
%       year: integer year
%       day: integer day that you would like to get RCD for
%       ampm: 0 = am, 1 = pm
%       dscd_S: structure described in header of read_QDOAS.m
%       sza_range: [sza_min, sza_max]
%       tag: string tag used to ID output AMF file (e.g. "no2vis")
%       lambda: wavelength of desired AMFs
%     
%
% OUTPUT:
%       rcd_vec: vector containing the following
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
%        amf: corresponding interpolated air-mass factor

if nargin < 7,
    lambda = 412;
end


% change this is running on a new computer
% / works for paths on Windows as well
if NO2_AMF_version == 1
    amf_dir = [code_path '/AMF_LUT/no2_amf_lut_v0_1/'];
elseif NO2_AMF_version ==2
    amf_dir = [code_path '/AMF_LUT/no2_amf_lut_v1_0/'];
end

rcd_vec = [];
amf = [];
% set am/pm string-name and sort things by SZA
if ampm == 0; ampm_str = 'am'; else ampm_str = 'pm'; end

% find indices for desired twilight and time of day
ind = find(dscd_S.day == day & dscd_S.ampm == ampm);
    
% do some early error catching
if isempty(ind)
    disp(['No DSCDs for day ' num2str(day) ampm_str])
    return
end

disp(['Calculating RCD for day: ' num2str(day) ' ' ampm_str])

% and make the amf output file name and the directory name for the amfs
out_file = ['AMF/amf_no2_v' num2str(NO2_AMF_version) '_' tag '_'...
    num2str(year) '_' num2str(day) '_' ampm_str '.dat'];% the working folder should have /AMF sub folder, Xiaoyi



% print the ozone data to the sza input file
fid = fopen([amf_dir 'DAY_SZA.dat'], 'w');
fprintf(fid, '%u\t%2.4f\t%2.4f\n',...
    [year*ones(length(ind),1) dscd_S.fd(ind) dscd_S.sza(ind)]');
fclose(fid);

% now print up to input no2 file
fid = fopen([amf_dir 'input_file_no2_amf.dat'], 'w');
fprintf(fid, '%s\n', '*Input file for NO2 AMF interpolation program');
fprintf(fid, '%s\n', '*');
fprintf(fid, '%s\n', '*Wavelength (350-550 nm) ?');
fprintf(fid, '%s\n', num2str(lambda));
fprintf(fid, '%s\n', '*Latitude (-90 (SH) to +90 (NH)) ?');
fprintf(fid, '%s\n', '80.05');
fprintf(fid, '%s\n', '*Longitude (-180 (- for W) to +180 (+ for E)) ?');
fprintf(fid, '%s\n', '-86.42');
fprintf(fid, '%s\n', '*Ground albedo flag: 1 for Koelemeijer dscd_vecbase and 2 for albedo value defined by the user');
fprintf(fid, '%s\n', '1');
fprintf(fid, '%s\n', '*Ground albedo value (if albedo flag = 1, put -99)');
fprintf(fid, '%s\n', '-99');
fprintf(fid, '%s\n', '*Name of the file with SZA values for interpolation (less than 30 characters) ?');
fprintf(fid, '%s\n', 'DAY_SZA.dat');
fprintf(fid, '%s\n', '*Interpolation results appearing on the screen: 1 -> yes, 0 -> no');
fprintf(fid, '%s\n', '1');
fclose(fid);

cur_dir = pwd;
cd(amf_dir);

% check operating system and run LUT executable accordingly
if ispc
    if NO2_AMF_version == 1
        executable = [amf_dir 'no2_amf_interpolation_dos.exe'];% for the old version of AMF LUT
    elseif NO2_AMF_version == 2
        executable = [amf_dir 'no2_amf_interpolation_v1_0.exe'];
    end
elseif isunix
    if NO2_AMF_version == 1
        executable = ['wine ' amf_dir 'no2_amf_interpolation_dos.exe'];% for the old version of AMF LUT
    elseif NO2_AMF_version == 2
        executable = ['wine ' amf_dir 'no2_amf_interpolation_v1_0.exe'];
    end
end
    
[status, result] = dos(executable, '-echo');

% copy output to desired file-name.  Read output sza, o3, and AMF
cd(cur_dir);
%fid=fopen(out_file,'w');
%fclose(fid);

copyfile([amf_dir 'no2_amf_output.dat'], out_file);
%rename(f,'out_file.dat',out_file);


[tmp_year, tmp_day, tmp_sza, amf] = textread(out_file, '%f%f%f%f', 'headerlines',5);


% make plots and get langley dscd_vec for morning and afternoon
close all;
figure;
print_setting(1,0,['temp']);
rcd_vec = langley(day, ampm, dscd_S.fd(ind), dscd_S.sza(ind),...
    dscd_S.saa(ind), dscd_S.mol_dscd(ind), dscd_S.rms(ind), amf, sza_range, filter_tag,dscd_S.err);

end