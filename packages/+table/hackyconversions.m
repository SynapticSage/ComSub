function tab = hackyconversions(tab)
% just performs a few hacky conversions on tables for datatype compatability purposes

tab = stringify(tab, 'numWindowsCut');
tab = stringify(tab, 'cutoffs');
tab = stringify(tab, 'lowcutoffs');
tab = stringify(tab, 'patternNames');
tab = stringify(tab, 'patternNamesFull');

function tab = stringify(tab, fieldname)
    if ismember(fieldname,fieldnames(tab)) && numel(tab.(fieldname)(1,:)) > 1
        disp("stringifying " + fieldname)
        tmp = tab.(fieldname);
        tmp = num2cell(tab.(fieldname),2);
        tmp = cellfun(@string, tmp, 'UniformOutput', false);
        tmp = cellfun(@(x)strjoin(x,'__'), tmp, 'UniformOutput', false);
        tmp = string(tmp);
        tab.(fieldname) = tmp;
    end
