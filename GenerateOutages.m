function outages = GenerateOutages(start_time, end_time, len, interval, num, file_name)
% This script is used to generate outage file
% start_time: start time of first outage
% end_time: end time of last outage
% len: length of single outage
% num: number of outages
% file_name: file name to save outage file

outages = [];
outage_start = start_time;
for i = 1:num
    outage_end = outage_start+len;
    if outage_end > end_time
        break;
    end
    outage = [i, outage_start, len];
    outages = [outages;outage];
    
    outage_start = outage_end + interval;
end

fid = fopen(file_name,'w');
for i = 1:size(outages, 1)
    fprintf(fid, '%d,%.1f,%d\n', outages(i,:));
end
