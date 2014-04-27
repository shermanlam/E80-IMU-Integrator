% sec - the number of sections at the beginning of the data over which to
% get the average
% fs - sampling frequency
% data - the array of data points
function ave = average(sec,fs,data)
    index = sec*fs;
    subset = data(1:index);
    ave = mean(subset);
end