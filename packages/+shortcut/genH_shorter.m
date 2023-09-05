function genH = genH_shorter(genH)
% genH_shorter(genH) - shorten the genH file to make it easier to read
%

genH = replace(genH, "fromCoherence  fromRipTimes", "coh");
genH = replace(genH, "fromSpectra  fromRipTimes", "spec");
genH = replace(genH, "fromWpli  fromRipTimes", "wpli");
genH = replace(genH, "coherence", "coh");
genH = replace(genH, "spectra", "spec");

