function [dscd_S, rcd_S] = get_all_rcds_v2016(dscd_S, species,...
    sza_range, pause, o3_flag, tag, lambda,code_path, filter_tag)
%  [dscd_S, rcd_S] = get_all_rcds(dscd_S, species, sza_range, pause, o3_flag, tag)

%ex:NO2 : [dscd_S,rcd_S]= get_all_rcds(dscd_S, 1, [86,91],0, 2, 'L2',437);
%  Ozone:  [dscd_S,rcd_S]= get_all_rcds(dscd_S, 0, [86,91],0, 2, 'L2',505);

%  This function runs langley plots of all data in a dscd vector
%  INPUT: 
%      dscd_S: standard DSCD structure (see read_QDOAS.m for details)
%      species: (0 = ozone, 1 = no2)
%      sza_range: [sza_min, sza_max]
%      pause: 0 = don't pause, 1 = do pause
%      o3_flag: this flag tells us whether the ozone values input to the
%            calculation code are SCDs or
%            VCDs.  1 for O3 VCD in DU and 2 for O3 SCD in molec/cm2.
%      tag: string tag is appended to the AMF file-name to identify it.
%            Typical tags would be "L1" or "L2"
%      lambda: optional (for variou NO2 fit regions!)
%  OUTPUT:
%      dscd_S: updated with dscd_S.amf entry
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
%             rcd_S.nbr -- number of data points used in calculation

%%%%%%%%% select Ozone/NO2 NDACC AMF version %%%%%%%%%
    O3_AMF_version = 2; % 1 = v1.0 , 2 = v2.0
    NO2_AMF_version = 2; % 1 = v0_1, 2 = v1_0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if species == 1
    if NO2_AMF_version == 1
        version_str = 'v0_1';
    elseif NO2_AMF_version == 2
        version_str ='v1_0';
    end
    disp(['Notice: You have selected NDACC NO_2 AMF version ' version_str]);
    disp('If want use different AMF, please check get_all_rcds_v2016.m');
    str = input(['Continue?: Y/N [Y]: '],'s');
    if strcmp(str,'N')
        return;
    end
elseif species == 0
    if O3_AMF_version == 1
        version_str = 'v1_0';
    elseif O3_AMF_version == 2
        version_str ='v2_0';
    end
    disp(['Notice: You have selected NDACC O_3 AMF version ' version_str]);
    disp('If want use different AMF, please check get_all_rcds_v2016.m');

%    str = input(['Continue?: Y/N [Y]: '],'s');
%    if strcmp(str,'N')
%        return;
%    end
end

if nargin < 7,
    lambda = 412;
end

rcd_vec = [];
dscd_S.amf = [];
for day = 1:366,
    for ampm = 0:1
        if species == 0,
            [rcd_vec_ln, amf] = ...
                rcd_ndacc_o3_v2016(dscd_S.year(1), day, ampm, dscd_S,sza_range, o3_flag, tag,code_path,O3_AMF_version, filter_tag, lambda);
        elseif species == 1,
            [rcd_vec_ln, amf] = ... 
                rcd_ndacc_no2_v2016(dscd_S.year(1), day, ampm, dscd_S, sza_range, tag, lambda,code_path,NO2_AMF_version);
        end
        rcd_vec = [rcd_vec; rcd_vec_ln];
        dscd_S.amf = [dscd_S.amf; amf];
    end
    if pause == 1; user_entry = input('Press enter to continue'); end
end
if ~isempty(rcd_vec)
    rcd_S.day = rcd_vec(:,1);
    rcd_S.ampm = rcd_vec(:,2);
    rcd_S.fd_min = rcd_vec(:,3);
    rcd_S.fd_max = rcd_vec(:,4);
    rcd_S.sza_min = rcd_vec(:,5);
    rcd_S.sza_max = rcd_vec(:,6);
    rcd_S.az_min = rcd_vec(:,7);
    rcd_S.az_max = rcd_vec(:,8);
    rcd_S.R2 = rcd_vec(:,9);
    rcd_S.m  = rcd_vec(:,10);
    rcd_S.m_1sigma = rcd_vec(:,11);
    rcd_S.y_int = rcd_vec(:,12);
    rcd_S.y_int_1sigma = rcd_vec(:,13);
    rcd_S.nbr  = rcd_vec(:,14);
else
    rcd_S = {};
end