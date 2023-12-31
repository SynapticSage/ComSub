function [obj, fields] = unnest(obj, nestfield, fields, unnestAddressChar)
% Unnests all fields in nestfield, or optionally, certain fields
%
% obj = unnest(obj, nestfield, fields, unnestAddressChar)
%
% Unnests all fields in nestfield, or optionally, certain fields. If
% fields is specified, then only those fields will be unnested. If
% unnestAddressChar is specified, then the unnested fields will be
% addressed by joining the nestfield and field names with the
% unnestAddressChar. If unnestAddressChar is not specified, then the
% unnested fields will be addressed by the field name alone.
%
% INPUTS
%   obj                 NdArray or struct
%   nestfield           string
%   fields              string or string array
%   unnestAddressChar   string (optional)
%
% OUTPUTS
%   obj                 NdArray or struct
%   fields              string array

if nargin < 3
    fields = [];
end
if nargin < 4
    unnestAddressChar = [];
end

if numel(obj) > 1
    %% NON-SCALAR : RECURSE
    %% --------------------
    
    obj = arrayfun(@(x) nd.unnest(x, nestfield, fields), obj, 'UniformOutput', false);
    obj = ndb.toNd(obj);
else
    %% WHAT TO DO WITH A SCALAR STRUCT
    %% -------------------------------
    
    if ~isfield(obj, nestfield) || ~isstruct(obj.(nestfield))
        obj = [];
        fields = [];
        return;
    end

    nFields = numel(fieldnames(obj.(nestfield)));
    if nargin <= 2 || isempty(fields)
        fields = string(fieldnames(obj.(nestfield)));
    else
        fields = string(fields);
    end

    if numel(fields) == nFields
        deleteNestField  = true;
    else
        deleteNestField = false;
    end

    for field = fields(:)'
        if isfield(obj.(nestfield), field)
            if isempty(unnestAddressChar)
                obj.(field) = obj.(nestfield).(field);
            else
                obj.(join([nestfield, field], unnestAddressChar)) = ...
                    obj.(nestfield).(field);
            end
            obj.(nestfield) = rmfield(obj.(nestfield),field);
        end
    end

    if deleteNestField
        obj = rmfield(obj, nestfield);
    end
end
