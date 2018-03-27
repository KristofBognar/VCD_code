function input_table = read_input_file_VCD()

input_table = table;

cd ..
cd('CF_package_local');

if ~exist('VCD_input_file.txt', 'file')
    error('Input file not found')
end

fid = fopen('VCD_input_file.txt','r');
%tline = fgetl(fid);

while ~feof(fid)
    tline = fgetl(fid);
    try
        eval(tline);
    catch
    end
end

if ~exist(input_table.plot_path,'dir'), mkdir(input_table.plot_path); end

copyfile('VCD_input_file.txt', [input_table.plot_path 'VCD_input_file.txt'], 'f');

% try
%     mkdir(input_table.plot_path);
%     status = copyfile('CF_input_file.txt', [input_table.plot_path 'CF_input_file_archive.txt'], 'f');
%     if status == 1
%         disp('input file has been archived');
%     else
%         disp('Warning: input file not archived, pls check "read_input_file.m" !');
%     end
% catch
%     disp('Warning: input file not archived, pls check "read_input_file.m" !');
% end