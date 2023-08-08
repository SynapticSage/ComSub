function tab = hackyconversions(tab)
% just performs a few hacky conversions on tables for datatype compatability purposes

if ismember('numWindowsCut',fieldnames(tab)) && numel(tab.numWindowsCut(1,:)) > 1
    tmp = tab.numWindowsCut
    tmp = num2cell(tab.numWindowsCut,2);
    tmp = cellfun(@string, tmp, 'UniformOutput', false);
    tmp = cellfun(@(x)strjoin(x,'-'), tmp, 'UniformOutput', false);
    tmp = string(tmp);
    tab.numWindowsCut = tmp;
end
