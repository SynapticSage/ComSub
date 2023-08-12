function [t_uv, t_spec] = trigspec(s, time_labels, frequency_labels, scalar_info)
    % TRIGSPEC Convert a struct array of trigspec structs to two tables
    if numel(s) > 1
        t_uv = cell(1, numel(s));
        t_spec = cell(1, numel(s));
        for i = 1:numel(s)
            if nd.isEmpty(s(i))
                continue
            end
            [t_uv{i}, t_spec{i}] = table.analyses.trigspec(s(i), time_labels, frequency_labels, scalar_info);
        end
        t_uv   = t_uv(~cellfun(@isempty, t_uv));
        t_spec = t_spec(~cellfun(@isempty, t_spec));
        t_uv   = vertcat(t_uv{:});
        t_spec = vertcat(t_spec{:});
    else

        s = nd.unnest(s, 'spec_avg');
        s = nd.broadcast(s);
        s = rmfield(s, ["threshold_crossed_times", "spec_stderr"]);
        if nargin < 4 || isempty(scalar_info)
            scalar_info = struct();
        end
        if nargin < 3 || isempty(frequency_labels)
            frequency_labels = 1:size(s(1).S1, 2);
        end
        if nargin < 2 || isempty(time_labels)
            time_labels = 1:size(s(1).S1, 1);
        end

        % Define the UV and spectral component field names
        % NOTE: ALL FIELDS HAVE TO BE SPECIFIED HERE, otherwise the setdiff below causes problems
        uv_fields = {'u_average', 'v_average', 'u_stderr', 'v_stderr', 'u_threshold', 'v_threshold'};
        common_fields = {'name', 'comp', 'direction', 'time_avg', 'epoch'};
        spec_fields = setdiff(fieldnames(s), [uv_fields, common_fields]);

        % Initialize the cell arrays for the UV and spectral component columns
        uv_columns = cell(1, numel(uv_fields) + numel(common_fields) + 2);
        uv_names = [uv_fields(:)', common_fields(:)', {'time', 'components'}];
        
        spec_columns = cell(1, numel(spec_fields) + numel(common_fields) + 2);
        spec_names = [spec_fields(:)', common_fields(:)', {'time', 'frequency'}];

        % Add the UV component columns
        min_2nd_dimension = structfun(@(x) size(x, 2), s);
        min_2nd_dimension = min(min_2nd_dimension);
        for i = 1:numel(uv_fields)
            field = uv_fields{i};
            array_uv = s.(field);
            if size(array_uv, 2) > min_2nd_dimension
                disp("Reshaping array_uv")
                array_uv = array_uv(:, 1:min_2nd_dimension);
            end
            assert(size(array_uv, 2) == min_2nd_dimension, 'The size of array_uv is not equal to the min_2nd_dimension.')
            array = reshape(array_uv, [], 1);
            uv_columns{i} = array;
        end

        % Add the spectral component columns
        shapes = cell(numel(spec_fields), 1);
        for i = 1:numel(spec_fields)
            field = spec_fields{i};
            array = reshape(s.(field), [], 1);
            spec_columns{i} = array;
            shapes{i} = size(s.(field));
        end

        % Add the common fields to both tables
        max_2nd_dimension = structfun(@(x) size(x, 2), s);
        max_2nd_dimension = max(max_2nd_dimension);
        for i = 1:numel(common_fields)
            field = common_fields{i};
            array_uv = reshape(s.(field)(:, 1:size(s.u_average, 2)), [], 1);
            array_spec = s.(field);
            if size(array_spec, 2) < max_2nd_dimension
                disp("Reshaping array_spec")
                array_spec = repmat(array_spec(:, 1), [1, max_2nd_dimension]);
            end
            disp("Size of array_spec: " + size(array_spec))
            assert(size(array_spec, 2) == max_2nd_dimension, 'The size of array_spec is not equal to the max_2nd_dimension.')
            array_spec = reshape(array_spec, [], 1);
            uv_columns{numel(uv_fields) + i} = array_uv;
            spec_columns{numel(spec_fields) + i} = array_spec;
        end
        % nums = cellfun(@numel, spec_columns);
        % nums = num2cell(nums);
        % assert(isequal(nums{:}), 'The number of rows in the spectral table is not equal to the number of rows in the UV table.')
        % q={spec_names{:}; spec_columns{:}}
        % q=struct(q{:})


        % Create the time and frequency labels
        uv_columns{end-1} = repmat(time_labels(:), [size(s.u_average, 2), 1]);
        uv_columns{end}   = reshape(repmat((1:size(s.u_average, 2)), [numel(time_labels), 1]), [], 1);

        spec_columns{end-1} = repmat(time_labels(:), [size(s.(spec_fields{1}), 2), 1]);
        spec_columns{end} = reshape(repmat(frequency_labels(:)', [numel(time_labels), 1]), [], 1);

        % Add scalar info to each table
        scalar_fields = fieldnames(scalar_info);
        for i = 1:numel(scalar_fields)
            field = scalar_fields{i};
            scalar_value = scalar_info.(field);
            uv_columns{end+1} = repmat(scalar_value, [numel(uv_columns{1}), 1]);
            uv_names{end+1} = field;
            spec_columns{end+1} = repmat(scalar_value, [numel(spec_columns{1}), 1]);
            spec_names{end+1} = field;
        end
        
        % Check that all columns have the same number of rows
        uv_lengths   = cellfun(@numel, uv_columns);
        spec_lengths = cellfun(@numel, spec_columns);
        
        if range(uv_lengths) ~= 0
            error('The UV columns do not all have the same number of rows.')
        end
        
        if range(spec_lengths) ~= 0
            error('The spectral columns do not all have the same number of rows.')
        end

        % Construct the tables
        t_uv   = table(uv_columns{:}, 'VariableNames', uv_names);
        t_spec = table(spec_columns{:}, 'VariableNames', spec_names);
    end
end

