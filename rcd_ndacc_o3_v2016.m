function [rcd_vec, amf] = ... 
    rcd_ndacc_o3_v2016(year, day, ampm, dscd_S,sza_range, o3_flag, tag,code_path,O3_AMF_version, filter_tag, lambda)

%  This function calculates the RCDs for ozone for a given day, SZA range, 
%  and dscd_vecset.
%  *** NOTE directory is hard-programed.  Change the amf_dir variable if
%  would like to set up in another directory structure.
%  *** NOTE have also hard-programmed Eureka.
%
% INPUT:
%    year: year
%    day: day number that would like to analyze
%    ampm: analyze am or pm twilight?  (am = 0, pm = 1)
%    dscd_S: this is an object that contains the fields
%               which are relevant to the Langley analysis.
%               It contains columns the following columns all of which
%               should have the same length!
%                   dscd_S.day: day in local time
%                   dscd_S.fd :fractional day (UTC)
%                   dscd_S.ampm: (am = 0, pm = 1)
%                   dscd_S.sza: SZA
%                   dscd_S.saa: solar azimuth
%                   dscd_S.rms : rms
%                   dscd_S.mol_dscd : focus molecule dscd
%                   dscd_S.err : focus molecule dscd fitting error
%                   dscd_S.x : fit of X cross-section
%                   dscd_S.o3 : ozone for input to radiative transfer
%       sza_range: this is a vector with two indices [sza_min, sza_max]
%                  representing the ideal langley plot range
%       o3_flag: this flag tells us whether the ozone values are SCDs or
%                VCDs.  1 for O3 VCD in DU and 2 for O3 SCD in molec/cm2.
%       tag: string tag is appended to the AMF file-name to identify it.
%            Typical tags would be "L1" or "L2"
%
% OUTPUT:
%       rcd_vec: vector containing the following entries
%             1 : day in local time
%             2 : am =0, pm = 1
%             3 : fd_min : minimum frac day of fits
%             4 : rcd_S.fd_max : maximum frac day of the fits
%             5 : sza_min : minimum SZA of fits
%             6 : sza_max : maximum SZA of the fits
%             7 : saa_min : minimum sol az of fits
%             8 : saa_max : maximum sol az of the fits
%             9 : r^2 value
%             10 : slope : slope of langley plot 
%             11 : slope_err : error in langley plot
%             12 :rcd : rcd
%             13 :rcd_err : rcd 1-sigma error in y-int
%             14 : number of points included in langley plot
%       amf : AMF corresponding with each index

% put LUT folder pareller to VCD code folder!
if ispc
    working_dir = pwd; % record current working dir
    mkdir('AMF');% make folder to save all AMF outputs from LUT
    cd(code_path);
    cd ..
    code_path = [pwd];% this should be the path for LUT folder!
    cd(working_dir);% go back to original working dir
end

% / works for paths on Windows as well
if O3_AMF_version == 1
    amf_dir = [code_path '/AMF_LUT/o3_amf_lut_v1_0/'];
elseif O3_AMF_version == 2
    amf_dir = [code_path '/AMF_LUT/o3_amf_lut_v2_0/'];
end

rcd_vec = [];
amf = [];
% set am/pm string-name and sort things by SZA
if ampm == 0; ampm_str = 'am'; else ampm_str = 'pm'; end

% find indices for desired twilight and time of day
ind = find(dscd_S.day == day & dscd_S.ampm == ampm);% Cirsten's code. Use local day of the year
%ind = find(floor(dscd_S.fd) == day & dscd_S.ampm == ampm); % modified by Xiaoyi, 2017-10-13, test to use UTC day of the year
    
% do some early error catching
if isempty(ind)
    disp(['No DSCDs for day ' num2str(day) ampm_str])
    return
end

if o3_flag==1
    condition=isnan(sum(dscd_S.o3(ind)));
elseif o3_flag==2
    condition=isnan(sum(dscd_S.mol_dscd(ind)));
end
    
if condition
    amf = ones(length(ind),1) * NaN;
    disp(['No input ozone VCDs for day ' num2str(day) ampm_str])
    return
end
disp(['Calculating RCD for day: ' num2str(day) ' ' ampm_str])

% Make amf output file name
output_file_nm = ['AMF/amf_o3_v' num2str(O3_AMF_version) '_' tag '_'...
    num2str(year) '_' num2str(day) '_' ampm_str '_' num2str(lambda) 'nm_.dat'];

% print the ozone data to the sza input file
L = length(dscd_S.day(ind));
fid = fopen([amf_dir 'sza_file_amf.dat'], 'w');

% Cristen: use total columns from sonde data as input
if o3_flag==1
    fprintf(fid, '%2.2f\t%2.2f\n', [dscd_S.sza(ind) dscd_S.o3(ind)]');
% Xiaoyi: use dSCDs as input 
elseif o3_flag==2, 
    fprintf(fid, '%2.2f\t%2.2f\n', [dscd_S.sza(ind) dscd_S.mol_dscd(ind)]'); 
end

fclose(fid);

% now print up to input o3 file
fid = fopen([amf_dir 'input_file_o3_amf.dat'], 'w');

fprintf(fid, '%s\n', '*Input file for O3 AMF interpolation program');
fprintf(fid, '%s\n', '*');
fprintf(fid, '%s\n', '*Wavelength (440-580 nm) ?');
fprintf(fid, '%s\n', num2str(lambda));
fprintf(fid, '%s\n', '*Day number (1-365 or 366 for leap year) ?');
fprintf(fid, '%d\n', day);
fprintf(fid, '%s\n', '*Latitude (-90 (SH) to +90 (NH)) ?');
fprintf(fid, '%s\n', '80.05');
fprintf(fid, '%s\n', '*Longitude (-180 (- for W) to +180 (+ for E)) ?');
fprintf(fid, '%s\n', '-86.42');
fprintf(fid, '%s\n', '*Ground albedo flag: 1 for Koelemeijer dscd_vecbase and 2 for albedo value defined by the user');
fprintf(fid, '%s\n', '1');
fprintf(fid, '%s\n', '*Ground albedo value (if albedo flag = 1, put -99)');
fprintf(fid, '%s\n', '-99');
fprintf(fid, '%s\n', '*Altitude of the station (between 0 and 4 km)');
fprintf(fid, '%s\n', '0.6 ');
fprintf(fid, '%s\n', '*Name of the file with SZA values for interpolation (less than 30 characters) ?');
fprintf(fid, '%s\n', 'sza_file_amf.dat');
fprintf(fid, '%s\n', '*O3 column interpolation: 1 for O3 VCD in DU and 2 for O3 SCD in molec/cm2');
fprintf(fid, '%d\n', o3_flag);
fclose(fid);

% change to amf directory to excecute interpolation scheme
cur_dir = pwd;
cd(amf_dir);

% check operating system and run LUT executable accordingly
if ispc
    if O3_AMF_version == 1
        executable = [amf_dir 'o3_amf_interpolation_dos.exe'];
    elseif O3_AMF_version == 2
        executable = [amf_dir 'o3_amf_interpolation_v2_0_dos.exe'];
    end
elseif isunix
    if O3_AMF_version == 1
        executable = ['wine ' amf_dir 'o3_amf_interpolation_dos.exe'];
    elseif O3_AMF_version == 2
        executable = ['wine ' amf_dir 'o3_amf_interpolation_v2_0_dos.exe'];
    end
end

try
    [status, result] = dos(executable, '-echo');
catch
    if isunix
        error('Please install Wine or another windows emulator');
    elseif ispc
        error('Could not run LUT executable')
    end
end
% copy output to desired file-name.  Read output sza, o3, and AMF
cd(cur_dir);
%copyfile([amf_dir 'o3_amf_output.dat'], output_file_nm);
copyfile([amf_dir 'o3_amf_output.dat'], output_file_nm);

%[tmp1_sza, tmp2_o3, amf] = textread(output_file_nm, '%f%f%f', 'headerlines',8);
[~, ~, amf] = textread(output_file_nm, '%f%f%f', 'headerlines',8);

% make plots and get langley dscd_vec for morning and afternoon
close all;
figure;
print_setting(1,0,['temp']);
rcd_vec = langley(day, ampm, dscd_S.fd(ind), dscd_S.sza(ind),...
    dscd_S.saa(ind), dscd_S.mol_dscd(ind), dscd_S.rms(ind), amf, sza_range, filter_tag,dscd_S.err);

end