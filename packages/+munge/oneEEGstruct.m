function tmp = oneEEGstruct(avgeeg)
% tmp = oneEEGstruct(avgeeg)
%
% one long eeg struct from a cell array of eeg structs

if iscell(avgeeg)
    avgeeg = ndb.toNd(avgeeg);
end
tmp = nd.cat(avgeeg, 1, 2);
tmp = tmp(end);
tmp.delta = nd.cat(tmp.delta, 1, 1);
tmp.theta = nd.cat(tmp.theta, 1, 1);
tmp.ripple = nd.cat(tmp.ripple, 1, 1);
tmp.starttime = sort(tmp.starttime);
tmp.starttime = tmp.starttime(1);
tmp.endtime = sort(tmp.endtime);
tmp.endtime = tmp.endtime(end);
tmp.samprate = median(tmp.samprate);
tmp.time = geteegtimes(tmp);
tmp.data = tmp.data(:);


if size(tmp.data, 1) == size(tmp.theta.data, 1)
    tmp.theta  = tmp.theta.data(:,1);
    tmp.delta  = tmp.delta.data(:,1);
    tmp.ripple = tmp.ripple.data(:,1);
else
    % Load filters
    disp("Sizes unequal, loading filters");
    load('thetafilter.mat');
    load('deltafilter.mat'); % assuming you have a deltafilter.mat similar to thetafilter.mat
    load('ripplefilter.mat'); % assuming you have a ripplefilter.mat similar to thetafilter.mat
    disp("Filtering");
    % Apply the theta, delta, and ripple filters to tmp.data using filtfilt
    tmp.theta = filtfilt(thetafilter.tf.num, thetafilter.tf.den, tmp.data);
    tmp.delta = filtfilt(deltafilter.tf.num, deltafilter.tf.den, tmp.data);
    tmp.ripple = filtfilt(ripplefilter.kernel, 1, tmp.data);
end
%
% tmp.theta.starttime = tmp.starttime(1);
% tmp.theta.endtime = tmp.endtime(end);
% tmp.theta.samprate = tmp.samprate;
% tmp.theta.time = tmp.time;
%
% tmp.delta.starttime = tmp.starttime(1);
% tmp.delta.endtime = tmp.endtime(end);
% tmp.delta.samprate = tmp.samprate;
% tmp.delta.time = tmp.time;
%
% tmp.ripple.starttime = tmp.starttime(1);
% tmp.ripple.endtime = tmp.endtime(end);
% tmp.ripple.samprate = tmp.samprate;
% tmp.ripple.time = tmp.time;
