function [VCD_table, dscd_S, rcd_S, avg_vcd, qdoas_filt, VCD_table2, dscd_S2, rcd_S2, avg_vcd2] = DSCD_to_VCD(year,code_path,plot_path,save_fig,QDOAS_data_dir,QDOAS_data_file,sonde)
% this is a function package to convert QDOAS output (DSCD) to NDACC required VCD
% writen by Xiaoyi Zhao 18/04/2016

% eg : 
% [VCD,VCD_CF, dscd_S, dscd_SCF, rcd_S, rcd_SCF, avg_vcd, avg_vcdCF, qdoas_filt] = DSCD_to_VCD(year,VCD_code_path,plot_path,save_fig,QDOAS_data_dir,QDOAS_data_file,sonde);

% Notes: 
% 1. this function needs few sub-functions, all of them should be placed in the directory of code_path
% 2. this function needs NDACC airmass factor look up table (AMF_LUT)
% 3. this function needs read in QDOAS output file in a defualt format (otherwise you need change the column structure in the input_sampe_2016.mat)
% 4. if this is the first time you use this package, I suggest you to go through each sub-functions (they are already properly commented) to make sure you understand the structure of this package (to avoid misusage)

% Output format: the main data output of this code is a VCD table, which
% contents:
%        1. year
%        2. day: julian day in local time
%        3. ampm: 0=am, 1=pm
%        4. fd_min: minimum fractional day UTC
%        5. fd_max: maximum fractional day UTC
%        6. fd: weighted mean fractional day UTC
%        7. sza_min: minimum solar zenith angle
%        8. sza_max: maximum solar zenith angle
%        9. sza: weighted mean solar zenith angle
%        10. saa_min: minimum solar azimuth angle
%        11. saa_max: maximum solar azimuth angle
%        12. saa: weighted mean solar azimuth angle
%        13. mean_vcd: twilight's VCD
%        14. sigma_mean_vcd: mean systematic error in VCDs
%        15. std_vcd: standard deviation of VCDs + other random error components
%        16. sigma_w_vcd: weighted mean error of VCDs
%        17. langley_vcd: slop of the langley fit
%        18. langley_vcd_err: err in slop of the langley fit
%
% Kristof Bognar: modified to return error budget based on Hendrick et al.,
% 2011. For old error budget, uncomment 2016 versions of functions in step 4

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
    
    % add tag to saved file if submission is for RD data
    batch_tag='';
    if any(strcmpi(input_table.Properties.VariableNames,'batch'))
        batch_tag=['_' input_table.batch];        
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
%VCD_code_path = 'F:\Work\VCD\VCD_code\VCD_code_2017';
%addpath(code_path);

% the place you want save the output
%working_dir = 'F:\Work\GBS\UT_GBS\Eureka2015\O3_reprocess';
%working_dir = 'F:\Work\GBS\UT_GBS\Eureka2015\NO2_reprocess';
%working_dir = 'H:\work\Eureka\GBS\CI\2011\UTGBS\VCD3';
if CF_run 
    working_dir = strcat(plot_path, 'VCD');
else
    working_dir = plot_path;
end

if ~exist(working_dir, 'dir'), mkdir(working_dir); end

% the place you have the QDOAS output data
%QDOAS_data_dir = 'F:\Work\QDOAS\Output\good\2015';
%QDOAS_data_dir = 'H:\work\Eureka\GBS\CI\2011\UTGBS\plots\CI_callibration';

% the name of the QDOAS output file
%QDOAS_data_file = 'UTGBS2015_zenith_43_128.dat';
%QDOAS_data_file = 'u1_2011_cal_CI_allflags_tests.mat';
%QDOAS_data_file = 'u1_2011_cal_CI_allflags.mat'; % CI&smooth flag SZA<85 
%QDOAS_data_file = 'u1_2011_cal_CIflag85_smoothflag92.mat'; % CI flag SZA<85, smooth flag <92 


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
[dscd_S, qdoas_filt, qdoas_raw, col] = read_QDOAS_v2017(QDOAS_data_file, col_o3_3, filt,3,save_fig,working_dir, trace_gas, CF_run);
    

%[dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS_v2016(QDOAS_data_file, col_o3_3, filt_good,1,save_fig,working_dir);

% elseif trace_gas == 2
%     if instrument == 1 % we retrieve NO2 from UT-GBS Vis data (follow NDACC recommendation)
%         [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS_v2016(QDOAS_data_file,col_no2_3, filt_good,1,save_fig,working_dir);
%     elseif instrument == 2 % we retrieve NO2 from PEARL-GBS UV data
%         [dscd_S, qdoas_filt, qdoas_raw] = read_QDOAS_v2016(QDOAS_data_file,col_no2_p0, filt_good,1,save_fig,working_dir);
%     end
% end

dscd_S = time_from_SZA_dscd_S(instrument, dscd_S); % this one will add time for the ref spec. 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 4. convert Ozone DSCDs to VCDs %%%%%
cd(working_dir);
if ~exist('AMF', 'dir'), mkdir('AMF'); end

if instrument==1 
    if trace_gas==1 || trace_gas==2
        instr_str='u1';
    elseif trace_gas==3
        instr_str='u0';
    end
elseif instrument==2
    if trace_gas==1 || trace_gas==2
        instr_str='p1';
    elseif trace_gas==3
        instr_str='p0';
    end
end    
    
if trace_gas == 1
    filter_tag = '';
% %     [dscd_S, rcd_S, avg_vcd] = get_ozone_vcds_v2016(dscd_S, sonde, instr_str, sza_range_ozone,lambda_ozone, save_fig,working_dir,code_path,filter_tag);
    [dscd_S, rcd_S, avg_vcd] = get_ozone_vcds_v2018(dscd_S, sonde, instr_str, sza_range_ozone,lambda_ozone, save_fig,working_dir,code_path,filter_tag);
    [avg_vcd,rcd_S] = assign_refspec_time_to_avgvcd(dscd_S,avg_vcd,rcd_S);% this will assgin the ref spec time and sza to rcd_S and avg_vcd
   
    if CF_run
        disp('Traditional Langley fits finished ... wait for 120 s ... ');pause(120);
        
        filter_tag = 'nocloud';
        error('modify get_ozone_vcds_v2017.m to include new error budget')
        [dscd_S2, rcd_S2, avg_vcd2] = get_ozone_vcds_v2017(dscd_S, sonde,instr_str,sza_range_ozone,lambda_ozone, save_fig,working_dir,code_path,filter_tag);
        [avg_vcd2,rcd_S2] = assign_refspec_time_to_avgvcd(dscd_S2,avg_vcd2,rcd_S2);% this will assgin the ref spec time and sza to rcd_S and avg_vcd
        disp('Cloud-screened Langley fits finished ... wait for 120 s ... ');pause(120);
    end
    
elseif trace_gas == 2
    filter_tag = '';
% %     [dscd_S, rcd_S, avg_vcd] = get_no2_vcds_v2016(dscd_S, lambda_no2, instr_str, sza_range_no2, save_fig,working_dir,code_path,filter_tag);
    [dscd_S, rcd_S, avg_vcd] = get_no2_vcds_v2018(dscd_S, lambda_no2, instr_str, sza_range_no2, save_fig,working_dir,code_path,filter_tag);
elseif trace_gas == 3
    filter_tag = '';
% %     [dscd_S, rcd_S, avg_vcd] = get_no2_vcds_v2016(dscd_S, lambda_no2_UV, instr_str, sza_range_no2, save_fig,working_dir,code_path,filter_tag);
    [dscd_S, rcd_S, avg_vcd] = get_no2_vcds_v2018(dscd_S, lambda_no2_UV, instr_str, sza_range_no2, save_fig,working_dir,code_path,filter_tag);

else
    disp('Warning: wrong set up in input file!');
    disp('Please double check input settings for this code!');
end

[avg_vcd] = add_langley_vcd(rcd_S, avg_vcd); % we include the slop in langley fit as VCD
if CF_run, [avg_vcd2] = add_langley_vcd(rcd_S2, avg_vcd2); end % we include the slop in langley fit as VCD
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 5. dispaly VCDs and print data %%%%%%%
if trace_gas == 1
    trace_gas_nm ='ozone';
elseif trace_gas==2
    trace_gas_nm = 'no2';
elseif trace_gas==3
    trace_gas_nm = 'no2uv';
end

VCD_table = struct2table(avg_vcd); % convert VCD output (structure) to table format
N = size(VCD_table);
%VCD_table.year = repmat(year,[N(1),1]);
year_table = table(repmat(year,[N(1),1]),'VariableNames',{'year'});
%VCD_table = [VCD_table(:,16),VCD_table(:,1:15)];% this is just move year to the 1st column
%VCD_table = [VCD_table(:,18),VCD_table(:,1:17)];% this is just move year to the 1st column
VCD_table = [year_table,VCD_table];

if CF_run 
    VCD_table2 = struct2table(avg_vcd2); % convert VCD output (structure) to table format
    N = size(VCD_table2);
    %VCD_table2.year = repmat(year,[N(1),1]);
    year_table2 = table(repmat(year,[N(1),1]),'VariableNames',{'year'});
    %VCD_table2 = [VCD_table2(:,16),VCD_table2(:,1:15)];% this is just move year to the 1st column
    %VCD_table = [VCD_table(:,18),VCD_table(:,1:17)];% this is just move year to the 1st column
    VCD_table2 = [year_table2,VCD_table2];
end

%% plot timeserise
figure;hold all; % plot VCDs as a function of fractional day of the year
gscatter(VCD_table.fd,VCD_table.mean_vcd,VCD_table.ampm,'br','..');
if CF_run
    gscatter(VCD_table2.fd,VCD_table2.mean_vcd,VCD_table2.ampm,'br','..');
    legend('V1 a.m.','V1 p.m.','V2 a.m.','V2 p.m.');
end
ylabel([trace_gas_nm ' molec/cm^2']);
xlabel('Day of the year (UTC)');
print_setting(1/2,save_fig,[trace_gas_nm 'VCDs']);


%% print data to tables
fid = fopen([trace_gas_nm 'VCDs.csv'],'w+');% print the VCD data to txt file
writetable(VCD_table,[trace_gas_nm 'VCDs.csv'],'Delimiter',',');
if CF_run
    fid2 = fopen([trace_gas_nm 'VCDs2.csv'],'w+');% print the VCD data to txt file
    writetable(VCD_table2,[trace_gas_nm 'VCDs2.csv'],'Delimiter',',');
end

%% save final variables
if ~CF_run
    cd('../../')

    if trace_gas == 1
        trace_gas_nm ='_O3_';
    elseif trace_gas==2
        trace_gas_nm = '_NO2_';
    elseif trace_gas==3
        trace_gas_nm = '_NO2_UV_';
    end
    savename=[input_table.instrument trace_gas_nm 'VCD_' input_table.year batch_tag '.mat'];

    save(savename,'avg_vcd','dscd_S','qdoas_filt','rcd_S','VCD_table');
end
