function [avg_vcd,rcd_S] = assign_refspec_time_to_avgvcd(dscd_S,avg_vcd,rcd_S)
%% assign ref spec time to avg_vcd
N = size(avg_vcd.day);
disp('Assign reference spec time to VCD structure ...');
for i = 1:1:N(1)
    TF = dscd_S.day == avg_vcd.day(i);
    try 
        if sum(TF) > 0 
            ref_sza = dscd_S.ref_sza(TF,:);
            ref_sza_utc1 = dscd_S.ref_sza_utc1(TF,:);
            ref_sza_utc2 = dscd_S.ref_sza_utc2(TF,:);
        else
            disp('Warning: a day number in avg_vcd is not found in dscd_S! Pls check these two structure!');
        end

        avg_vcd.ref_sza(i,:) = ref_sza(1);
        avg_vcd.ref_sza_utc1(i,:) = ref_sza_utc1(1);
        avg_vcd.ref_sza_utc2(i,:) = ref_sza_utc2(1);
    catch
        disp('Warning: Failled to assign reference spec time to VCD structure!');
        disp('Warning: No reference spec time will be stored to the final VCD files!');
    end
        
end
%% assign ref spec time to rcd_S
N = size(rcd_S.day);
disp('Assign reference spec time to RCD structure ...');
for i = 1:1:N(1)
    TF = dscd_S.day == rcd_S.day(i);
    try 
        if sum(TF) > 0 
            ref_sza = dscd_S.ref_sza(TF,:);
            ref_sza_utc1 = dscd_S.ref_sza_utc1(TF,:);
            ref_sza_utc2 = dscd_S.ref_sza_utc2(TF,:);
        else
            disp('Warning: a day number in rcd_S is not found in dscd_S! Pls check these two structure!');
        end

        rcd_S.ref_sza(i,:) = ref_sza(1);
        rcd_S.ref_sza_utc1(i,:) = ref_sza_utc1(1);
        rcd_S.ref_sza_utc2(i,:) = ref_sza_utc2(1);
    catch
        disp('Warning: Failled to assign reference spec time to RCD structure!');
    end
        
end