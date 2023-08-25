fommsubspaceToPath
% Grab fields from cgramcnew and place into existing efizz structure
% (rather than from scratch)

% grabfields = ["Cavg", "Ctoppair"];
% grabfields = ["wpli_avg", "wpli_toppair"];
% grabfields = ["phi"];
grabfields = ["wpli_avg", "phi"];
const = option.constants();
animal_list = ["ER1"]
animal = animal_list{1};
for i = progress(1:length(animal_list), 'Title', 'FieldsToEfizz')
    animal = animal_list{i};
    m = matfile(animal + "spectralBehavior.mat", 'Writable', true);
    efizz = m.efizz;
    assert(ndbFile.exist(animal, 'cgramcnew'), "cgramcnew does not exist for " + animal);
    % cgramcnew = ndb.load(animal, 'cgramcnew', 'inds', [1,2]);
    cgramcnew = ndb.load(animal, 'cgramcnew', 'asNd', true);
    cgramcnew = cgramcnew(1, 2:2:end);
    cgramcnew = nd.cat(cgramcnew, 1, [], 'removeEmpty', true);
    
    if ~isempty(cgramcnew(1).t)
        tc=numel(cgramcnew.t);
    else
        tc=size(cgramcnew.S1,1);
    end
    te=numel(efizz.t);
    assert(all(tc==te), "time vectors do not match for " + animal);

    for j = 1:length(grabfields)
        field = grabfields(j);
        efizz.(field) = cgramcnew.(field);
    end

    m.efizz = efizz;
    
end
!pushover-cli "FieldsToEfizz done"
% RunAll_wpli
