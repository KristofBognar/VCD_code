function [VCD_table, dscd_S_tmp, rcd_S, avg_vcd, qdoas_filt, VCD_table2, dscd_S2, rcd_S2, avg_vcd2] = DSCD_to_VCD_modify_dscd_S(year,code_path,plot_path,save_fig,QDOAS_data_dir,QDOAS_data_file,sonde)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 0. input options %%%%%%%%%%%%%%%%%%%%%%%%

% check whether code is called by CF package (all inputs given), or as a
% standalone function (no input given)
if nargin==0
    % if called as standalone function, set input parameters here
    
    % only run VCD code once
    CF_run=false;
    
    % Change the current folder to the folder of this m-file.
    if(~isdeployed)
      cd(fileparts(which(mfilename)));
    end    
    
    % read input file
    input_table = read_input_file_VCD();
    
    % assign year
    year=str2double(input_table.year);
    
    % assign trace gas number
    trace_gas=input_table.tg;
    
    % assign instrument number
    if strcmp(input_table.instrument,'UT-GBS')
        instrument=1;
    elseif strcmp(input_table.instrument,'PEARL-GBS')
        instrument=2;
    end
    
    code_path=input_table.VCD_code_path;
    plot_path=input_table.plot_path;
    QDOAS_data_file=input_table.data_input;
    load(input_table.data_sonde);
    save_fig=input_table.save_fig;
    
    % check if input file is a .mat file (table)
    if isempty(strfind(QDOAS_data_file,'.mat'))
        error('QDOAS data must be in matlab table format');
    end

    % select filter
    try
        if strcmp(input_table.filter,'bad')
            goodfilt=false;
        else
            goodfilt=true;
        end
    catch
        goodfilt=true;
    end
    
elseif nargin~=7
    error('Need 7 input variables');
else
    CF_run=true;
    goodfilt=true;
    
    % trace gas type
    trace_gas = 1; % 1 = ozone, 2 = NO2

    % instrument
    instrument = 1; % 1 = UT-GBS, 2 = PEARL-GBS
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 1. input variables %%%%%%%%%%%%%%%%%%%%%%%%
% the place where you have the VCD code package
if CF_run 
    working_dir = strcat(plot_path, 'VCD');
else
    working_dir = plot_path;
end

if ~exist(working_dir, 'dir'), mkdir(working_dir); end

% Solar Zenith Angle Range for VCD calculation (please reference the guide line from NDACC manual)
sza_range_ozone = [86,91]; % for example, SZA from 86 to 90 degree
sza_range_no2 = [86,91]; % for example, SZA from 86 to 91 degree

% wavelength will be used in the NDACC AMF-LUT
lambda_ozone = 505; % the centre wavlenth used in QDOAS fitting for ozone
lambda_no2 = 437;% for UT-GBS in Vis
lambda_no2_UV = 365;% for PEARL-GBS in UV

% save output figures or not
%save_fig = 1; % 1 = yes, 0 = no
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 2. extra inputs %%%%%%%%%%%%%%%%%%%%%%%%%%
% the following are some extra improtant inputs which can not be set up in
% the input file above, please visit the sub functions and change to what
% you need

% 2.1 NDACC AMF-LUT version
% get_all_rcds_v2016.m ---> need select which NDACC AMF LUT will be used
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% 3. read in QDOAS data %%%%%%%
cd(code_path);
load('input_sample_2017_no_sonde.mat'); % use the default structure of QDOAS output
if CF_run, cd(QDOAS_data_dir); end
% if trace_gas == 1 % we ONLY retrieve ozone from UT-GBS Vis data! (follow NDACC recommendation)
    
% choose filter
if goodfilt
    filt=filt_good;
    disp('Data read using good RMS filter')
    
    % set Different RMS filter for NO2 UV (from Cristen's thesis) 
    if trace_gas==3, filt.rms_vec=[[87,0.002];[92,0.003]]; end
else
    filt=filt_bad;
    disp('Data read using bad RMS filter')
    
    if trace_gas==3, filt.rms_vec=[[87,0.004];[92,0.006]]; end
end
    
% trace gas is selected in read_QDOAS_v2017 when reading from table input
[dscd_S_tmp, ~, ~, ~] = read_QDOAS_v2017(QDOAS_data_file, col_o3_3, filt,3,save_fig,working_dir, trace_gas, CF_run);
    

if ~CF_run
    cd('../../')

    if trace_gas == 1
        trace_gas_nm ='_O3_';
    elseif trace_gas==2
        trace_gas_nm = '_NO2_';
    elseif trace_gas==3
        trace_gas_nm = '_NO2_UV_';
    end
    savename=[input_table.instrument trace_gas_nm 'VCD_' input_table.year '.mat'];

    % modify saved file to add total integration timeinfo to dscd_S
    % DON'T overwrite dscd_S in saved file!
    
    tmp=matfile(savename,'Writable',true);
    
    dscd_S_orig=tmp.dscd_S;
    dscd_S_orig.tot_tint=dscd_S_tmp.tot_tint;
    tmp.dscd_S=dscd_S_orig;
    
    close all;
    
end

end

