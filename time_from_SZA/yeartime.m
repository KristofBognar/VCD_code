function [ year_time ] = yeartime( year )
%yeartime Calculate matlab time on Jan. 1 00:00 of given year

% matlab datenum runs from jan 1, 0000, 00:00:00 which is equal to 1!!

formstr='dd-mm-yyyy HH:MM:SS';

if max(size(year))==1
    datestr=['01-01-',num2str(year),' 00:00:00'];
    year_time=datenum(datestr,formstr);
else
    year_time=NaN(size(year));
    for i=1:max(size(year))
        datestr=['01-01-',num2str(year(i)),' 00:00:00'];
        year_time(i)=datenum(datestr,formstr);
    end
end
    
end

