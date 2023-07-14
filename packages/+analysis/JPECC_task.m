function Jpecc_Results = JPECC_task(Source, Target, Option)
% Computes lagged JPECC between patterns in two brain araes

%Patterns(i).X_source = [Neurons, Times, Trials]

% Number of cross validation folds.
cvNumFolds = Option.jpecc.cvnum;

for i = 1:size(Source,1)
for j = 1:size(Source,2)

    disp("----------")
    disp([i,j])
    disp("----------")

    S = Source{i,j};

for p = 1:size(Target,1)

    disp("----------")
    disp(p)
    disp("----------")

for q = 1:size(Target,2)

    T = Target{p,q};

    tens_curr_source = S.val;
    tens_curr_target = T.val;
    source_time = S.time;
    target_time = T.time;

    % sampling
    num_samples_source = size(tens_curr_source, 1);
    num_samples_target = size(tens_curr_target, 1);

    if num_samples_source > num_samples_target

        % Sort the time vectors
        [source_time, source_idx] = sort(source_time);
        [target_time, ~] = sort(target_time);

        % Initialize the selected source indices
        selected_source_idx = zeros(1, num_samples_target);

        % For each of the first n target times
        for d = 1:num_samples_target
            % Calculate the time differences between the current target time and all source times
            time_diffs = abs(source_time - target_time(d));
    
            % Find the source time that is closest to the current target time
            [~, closest_source_idx] = min(time_diffs);
    
            % Add the index of the closest source time to the selected source indices
            selected_source_idx(d) = source_idx(closest_source_idx);
    
            % Remove the closest source time from consideration in future iterations
            source_time(closest_source_idx) = [];
            source_idx(closest_source_idx) = [];
        end

        tens_curr_source = tens_curr_source(selected_source_idx, :);

    elseif num_samples_source < num_samples_target

        % Sort the time vectors
        [source_time, ~] = sort(source_time);
        [target_time, target_idx] = sort(target_time);

        % Initialize the selected source indices
        selected_target_idx = zeros(1, num_samples_source);

        % For each of the first n target times
        for d = 1:num_samples_source
            % Calculate the time differences between the current target time and all source times
            time_diffs = abs(target_time - source_time(d));
    
            % Find the source time that is closest to the current target time
            [~, closest_target_idx] = min(time_diffs);
    
            % Add the index of the closest source time to the selected source indices
            selected_target_idx(d) = target_idx(closest_target_idx);
    
            % Remove the closest source time from consideration in future iterations
            target_time(closest_target_idx) = [];
            target_idx(closest_target_idx) = [];
        end

        tens_curr_target = tens_curr_target(selected_target_idx, :);

    end

    cs = tens_curr_source;
    ct = tens_curr_target;

    nan_rows = any(isnan(cs), 2) | ...
               any(isnan(ct), 2);

    cs = cs(~nan_rows,:);
    ct = ct(~nan_rows,:);

    lessThan3Samples = sum(~nan_rows) <= 3;
    if lessThan3Samples
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val1    = nan;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p1    = nan;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val2    = nan;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p2    = nan;
        continue
    end

    N = size(cs,1); %number of trials
    cvp = cvpartition(N, 'KFold', cvNumFolds);
        
    % use cross-validation
    U1 = zeros(size(cs,1),1);
    V1 = zeros(size(cs,1),1);
    U2 = zeros(size(cs,1),1);
    V2 = zeros(size(cs,1),1);

    for k = 1:cvNumFolds
            
    % generate canonical dimension of each matrix on the training set  

        [A,B] = canoncorr(...
                cs(cvp.training(k),:),...
                ct(cvp.training(k),:));


        % project the test set onto the canonical dimension
        U1(cvp.test(k)) = cs(cvp.test(k),:)*A(:,1);
        V1(cvp.test(k)) = ct(cvp.test(k),:)*B(:,1);

        U2(cvp.test(k)) = cs(cvp.test(k),:)*A(:,2);
        V2(cvp.test(k)) = ct(cvp.test(k),:)*B(:,2);
  
   end

        % correlate the projections. since each test set's projections will
        % be zero mean already, we can just combine them all here
        [rval1, pval1] = corr(U1,V1);
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val1 = rval1;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p1 = pval1;
        [rval2, pval2] = corr(U2,V2);
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val2 = rval2;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p2 = pval2;

end
end

end
end
