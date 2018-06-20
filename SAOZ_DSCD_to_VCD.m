function [VCD_table, dscd_S, rcd_S, avg_vcd, qdoas_filt, VCD_table2, dscd_S2, rcd_S2, avg_vcd2] = SAOZ_DSCD_to_VCD(year,code_path,plot_path,save_fig,QDOAS_data_dir,QDOAS_data_file,sonde)
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

for twostep=[1,2]

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
    
    code_path=input_table.VCD_code_path;
    plot_path=input_table.plot_path;
    QDOAS_data_file=input_table.data_input;
    load(input_table.data_sonde);
    save_fig=input_table.save_fig;
    
    % check if input file is a .mat file (table)
    if isempty(strfind(QDOAS_data_file,'.mat'))
        error('QDOAS data must be in matlab table format');
    end

    % add tag to saved file if submission is for RD data
    batch_tag='';
    rd_run=0;
    if any(strcmpi(input_table.Properties.VariableNames,'batch'))
        batch_tag=['_' input_table.batch];
        rd_run=1;
    end
    
    
elseif nargin~=7
    error('Need 7 input variables');
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
working_dir = plot_path;

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
% lambda_no2 = 437;% for UT-GBS in Vis, 425-450 window
lambda_no2 = 457;% for UT-GBS in Vis, 425-490 window

if trace_gas==2
    disp(['Wavelength for AMF LUT set to ' num2str(lambda_no2)])
end

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
% if trace_gas == 1 % we ONLY retrieve ozone from UT-GBS Vis data! (follow NDACC recommendation)
    
    
load(QDOAS_data_file);
QDOAS_data=dscd(dscd.year==year,:);

% trace gas is selected in read_QDOAS_v2017 when reading from table input
[dscd_S, qdoas_filt, qdoas_raw, col] = SAOZ_read_QDOAS(QDOAS_data,trace_gas);
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 4. convert DSCDs to VCDs %%%%%
cd(working_dir);
if ~exist('AMF', 'dir'), mkdir('AMF'); end

instr_str='saoz';
filter_tag = '';

if trace_gas == 1
    error('Change code to use fix RCD');
    [dscd_S, rcd_S, avg_vcd] = get_ozone_vcds_v2018(dscd_S, sonde, instr_str, sza_range_ozone,lambda_ozone, save_fig,working_dir,code_path,filter_tag);
   
elseif trace_gas == 2

    % get VCDs with fix RCD, and also with daily RCD (dscd_S2 and avg_vcd2)
    [dscd_S, rcd_S, avg_vcd, dscd_S2, avg_vcd2] = SAOZ_get_no2_vcds(dscd_S, lambda_no2, instr_str, sza_range_no2, save_fig,working_dir,code_path,filter_tag);

else
    disp('Warning: wrong set up in input file!');
    disp('Please double check input settings for this code!');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 5. dispaly VCDs and print data %%%%%%%
if trace_gas == 1
    trace_gas_nm ='ozone';
elseif trace_gas==2
    trace_gas_nm = 'no2';
end

VCD_table = struct2table(avg_vcd); % convert VCD output (structure) to table format
N = size(VCD_table);
%VCD_table.year = repmat(year,[N(1),1]);
year_table = table(repmat(year,[N(1),1]),'VariableNames',{'year'});
%VCD_table = [VCD_table(:,16),VCD_table(:,1:15)];% this is just move year to the 1st column
%VCD_table = [VCD_table(:,18),VCD_table(:,1:17)];% this is just move year to the 1st column
VCD_table = [year_table,VCD_table];

VCD_table2 = struct2table(avg_vcd2); % convert VCD output (structure) to table format
N2 = size(VCD_table2);
%VCD_table.year = repmat(year,[N(1),1]);
year_table2 = table(repmat(year,[N2(1),1]),'VariableNames',{'year'});
%VCD_table = [VCD_table(:,16),VCD_table(:,1:15)];% this is just move year to the 1st column
%VCD_table = [VCD_table(:,18),VCD_table(:,1:17)];% this is just move year to the 1st column
VCD_table2 = [year_table2,VCD_table2];

%% plot timeserise
figure;hold all; % plot VCDs as a function of fractional day of the year
gscatter(VCD_table.fd,VCD_table.mean_vcd,VCD_table.ampm,'br','..');
gscatter(VCD_table2.fd,VCD_table2.mean_vcd,VCD_table2.ampm,'br','xx');
ylabel([trace_gas_nm ' molec/cm^2']);
xlabel('Day of the year (UTC)');
print_setting(1/2,save_fig,[trace_gas_nm 'VCDs']);


%% print data to tables
writetable(VCD_table,[trace_gas_nm 'VCDs.csv'],'Delimiter',',');

writetable(VCD_table2,[trace_gas_nm 'VCDs_dailyRCD.csv'],'Delimiter',',');

%% save final variables
cd('../')

if trace_gas == 1
    trace_gas_nm ='_O3_';
elseif trace_gas==2
    trace_gas_nm = '_NO2_';
end

if rd_run==1 && twostep==1 && trace_gas==1
    savename='vcd_1st_run.mat';
    vcd_1st_run=VCD_table;
    save(savename,'vcd_1st_run');
else
    savename=[input_table.instrument trace_gas_nm 'VCD_' input_table.year batch_tag '.mat'];
    save(savename,'avg_vcd','dscd_S','qdoas_filt','rcd_S','VCD_table',...
                  'avg_vcd2','dscd_S2','VCD_table2');

    return
end
    
    
end
    