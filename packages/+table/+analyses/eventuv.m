function t_event = eventuv(s, pattern_labels, event_labels, uv_labels, scalar_info)
    % EVENTSPEC Convert a struct array of event structs to a table

    first_nonempty = find(arrayfun(@(x) ~nd.isEmpty(x), s), 1);
    if nargin < 5 || isempty(scalar_info)
        scalar_info = struct();
    end
    if nargin < 4 || isempty(uv_labels)
        uv_labels = 1:size(s(first_nonempty).event_u_values, 3);
        uv_labels = uv_labels(:);
        % uv_labels = permute(uv_labels, [3, 2, 1]);
        % uv_labels = repmat(uv_labels, [size(s(first_nonempty).event_u_values, 1), size(s(first_nonempty).event_u_values, 2), 1]);
        % uv_labels = reshape(uv_labels, [], 1);
    end
    if nargin < 2 || isempty(pattern_labels)
        pattern_labels = 1:size(s(first_nonempty).event_u_values, 1);
        pattern_labels = pattern_labels(:);
        % pattern_labels = repmat(pattern_labels, [1, size(s(first_nonempty).event_u_values, 2), size(s(first_nonempty).event_u_values, 3)]);
        % pattern_labels = reshape(pattern_labels, [], 1);
    end

    if numel(s) == 1
        % Define the event field names
        event_fields = {'event_u_values', 'event_v_values', 'patterns', 'events', 'uv_components', 'event_time', 'epoch', 'lindist', 'vel', 'trajbound', 'correct'};
        if nargin < 3 || isempty(event_labels)
            % event_labels = 1:size(s(first_nonempty).event_u_values, 2);
            event_labels = 1:size(s.event_u_values, 2);
            event_labels = event_labels(:);
            % event_labels = repmat(event_labels, [size(s(first_nonempty).event_u_values, 1), 1, size(s(first_nonempty).event_u_values, 3)]);
            % event_labels = reshape(event_labels, [], 1);
        end

        s.patterns = pattern_labels(:);
        s.events = event_labels(:)';
        s.uv_components = permute(uv_labels(:), [3, 2, 1]);
        s = rmfield(s, setdiff(fieldnames(s), event_fields));
        q = nd.broadcast(s);
        sz = structfun(@size, q, 'UniformOutput', false);
        sz = struct2cell(sz);
        if ~isequal(sz{:})
            error('The event fields do not all have the same size.')
        end
        s = q;

        % Initialize the cell arrays for the event columns
        event_columns = cell(1, numel(event_fields));  
        event_names = event_fields;

        % Add the event columns
        for i = 1:numel(event_fields)
            field = event_fields{i};
            array = reshape(s.(field), [], 1);
            event_columns{i} = array;
        end

        % Add scalar info to the table
        scalar_fields = fieldnames(scalar_info);
        for i = 1:numel(scalar_fields)
            field = scalar_fields{i};
            scalar_value = scalar_info.(field);
            event_columns{end+1} = repmat(scalar_value, [numel(event_columns{1}), 1]);
            event_names{end+1} = field;
        end
        
        % Check that all columns have the same number of rows
        event_lengths = cellfun(@numel, event_columns);
        
        if range(event_lengths) ~= 0
            error('The event columns do not all have the same number of rows.')
        end

        % Construct the table
        t_event = table(event_columns{:}, 'VariableNames', event_names);
    else
        t_event = cell(1, numel(s));
        for i = 1:numel(s)
            if nd.isEmpty(s(i))
                continue
            end
            [i1,i2] = ind2sub(size(s), i);
            scalar_info.pattern_cca = i;
            scalar_info.pattern_cca1 = i1;
            scalar_info.pattern_cca2 = i2;
            t_event{i} = table.analyses.eventuv(s(i), pattern_labels, event_labels, uv_labels, scalar_info);
        end
        t_event = vertcat(t_event{:});
    end
end

