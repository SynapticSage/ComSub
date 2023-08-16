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
tmp.starttime = tmp.starttime(1);
tmp.endtime = tmp.endtime(end);
tmp.samprate = median(tmp.samprate);
tmp.time = geteegtimes(tmp);

tmp.theta.starttime = tmp.starttime(1);
tmp.theta.endtime = tmp.endtime(end);
tmp.theta.samprate = tmp.samprate;
tmp.theta.time = tmp.time;

tmp.delta.starttime = tmp.starttime(1);
tmp.delta.endtime = tmp.endtime(end);
tmp.delta.samprate = tmp.samprate;
tmp.delta.time = tmp.time;

tmp.ripple.starttime = tmp.starttime(1);
tmp.ripple.endtime = tmp.endtime(end);
tmp.ripple.samprate = tmp.samprate;
tmp.ripple.time = tmp.time;
