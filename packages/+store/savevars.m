function savevars(Option, Events, Spk, saveVars)

if nargin < 4
	saveVars = struct();
end
saveVars.Option = Option;
hash = store.gethash(Option);

thisFile = fullfile(hashdefine(), hash + ".mat");
disp("Saving ...");
tic; save(thisFile, '-struct', 'saveVars','-v7.3', '-nocompression');
disp("... " + toc + " seconds");
% link most recent state
pushd(hashdefine());
recencyName = Option.animal + "_" + replace(Option.generateH," ", "") + ...
	    "_zscore_" +  Option.preProcess_zscore + ...
	    "_midpattern_" + Option.midpattern + ...
	    "_mostRecentState.mat";
system("ln -sf " + hash + ".mat " + recencyName);
popd()
% save raw?
if Option.saveRaw
	disp("Saving raw...");
	tic; save(thisFile, "Events", "Spk",'-v7.3', '-append', '-nocompression');
	disp("... " + toc + " seconds");
else
	disp("Saving raw scrubbed")
	Events = nd.flexrmfield(Events, ["H","Hvals","H","times"]);
	tic; save(thisFile, "Events", '-v7.3', '-append', '-nocompression');
	disp("... " + toc + " seconds");
end

datetime_=datetime();
save(thisFile, "datetime_",'-v7.3', '-append', '-nocompression');
