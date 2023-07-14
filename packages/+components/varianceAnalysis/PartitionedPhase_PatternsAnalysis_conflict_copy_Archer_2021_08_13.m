%% variance by phase and patterns
%sv: strength_variance

% Obtain all of the keys
if ~ispc
    folder = "/Volumes/Colliculus/commsubspace/hash_forTesting/";
    keys = dir('/Volumes/Colliculus/commsubspace/hash_forTesting/*.mat');
else
    folder = 'C:\Users\BrainMaker\commsubspace\hash_forTesting\';
    keys = dir('C:\Users\BrainMaker\commsubspace\hash_forTesting\*.mat');
end

keys = string({ keys.name });
keys = keys(~contains(keys, "Poly"));
[Patterns, Behaviors, Options] = query.combinePatterns(folder + keys);
%%
clear beh
for a = 1:size(Behaviors,1)
    animal = Options(a).animal;
    for phase = 1:4
        beh{a, phase} = table.behavior.lookup(animal, Behaviors(a,1,1).compAnalysisByBeh(phase).component_time);
    end
end
%%
%Patterns = reshape(Patterns, size(Patterns,1) * size(Patterns,2), size(Patterns, 3), size(Patterns,4));

plot_high_patterns_only = true;
normalize               = false;
color_by = ["pattern", "component"];
% animal  = Option.animal;
if ispc
   savedir = "C:\Users\BrainMaker\commsubspace\figures\componentOverEpochTperf";
else
savedir = "/Volumes/Colliculus/commsubspace/figures/componentOverEpochTperf_allAnim";
end
mkdir(savedir);
patternNames = ["theta", "delta", "ripple"];
patternSym = shortcut.patternSymbols(patternNames);
TIME = 2;
COMP = 1;
K = 5;

%%
combined_hpc = cell(size(Behaviors,1), 4, size(Behaviors,2)); 
combined_pfc = cell(size(Behaviors,1), 4, size(Behaviors,2));
usegramm = true;
for a = 1:size(Behaviors,1)
    for phase = 1:4
        groups = findgroups(beh{a,phase}.epoch);
        num_epochs = numel(unique(groups));
        for p = 1:size(Behaviors,2)
            if isempty(Behaviors(a,p,1)) || all(structfun(@isempty, Behaviors(a,p,1)))
                continue
            end
            hpccomp = Behaviors(a,p,1).compAnalysisByBeh(phase).component_strength';
            pfccomp = Behaviors(a,p,2).compAnalysisByBeh(phase).component_strength';
            pfc.med = median(pfccomp, TIME);
            hpc.med = median(hpccomp, TIME);
            uGroups = unique(groups);
            for u = [1:numel(uGroups)]
                C_hpc = hpccomp(:,groups == u);
                C_pfc = pfccomp(:,groups == u);
                pat = @(x)  ceil(x ./ K);
                com = @(x)  mod(x-1, K) +  1;
                
                
                % For epoch analysis
                b.tperf(a,phase,u) = nanmean(beh{a,phase}.tperf_timewise(groups == u));
                for k = 1:size(hpccomp,1)
                    pfc.sv(a,p, phase, k, u,:) = var(bsxfun(@minus, C_pfc(k,:), pfc.med(k)),[],TIME);
                    hpc.sv(a,p, phase, k, u,:) = var(bsxfun(@minus, C_hpc(k,:), hpc.med(k)),[],TIME);
                    pfc.m(a,p,phase, k, u,:)   = mean(C_pfc(k,:), TIME);
                    hpc.m(a,p,phase, k, u,:)   = mean(C_hpc(k,:), TIME);
                    pfc.am(a,p,phase, k, u,:)  = mean(abs(C_pfc(k,:)), TIME);
                    hpc.am(a,p,phase, k, u,:)  = mean(abs(C_hpc(k,:)), TIME);
                end
                if usegramm
                    % For gramm
                    [KK,TT] = size(C_hpc);
                    N = numel(C_hpc);
                    tperf = beh{a, phase}.tperf_timewise(groups==u);
                    D = C_hpc'; D = D(:); 
                    combined_hpc{a,phase,p,u} = [repmat([a, u, phase], N ,1), repelem(pat(1:KK)',TT), repelem(com(1:KK)',TT), repmat(tperf, KK, 1), D];
                    D = C_pfc'; D = D(:);
                    combined_pfc{a,phase,p,u} = [repmat([a, u, phase], N ,1), repelem(pat(1:KK)',TT), repelem(com(1:KK)',TT), repmat(tperf, KK, 1), D];
                    if u == numel(uGroups)
                        disp(join(["finished", a, phase, p], " "))
                    end
                end
            end
        end
    end
end
if usegramm
    combined_hpc = cat(1, combined_hpc{:});
    combined_pfc = cat(1, combined_pfc{:});
    combined_hpc = num2cell(combined_hpc, 1);
    combined_pfc = num2cell(combined_pfc, 1);
    combined_hpc = table(combined_hpc{:}, 'VariableNames', {'animal','epoch','phase','pattern','component','tperf','val'});
    combined_pfc = table(combined_pfc{:}, 'VariableNames', {'animal','epoch','phase','pattern','component','tperf','val'});
    if downsample ~= 0
        filt = rand(height(combined_hpc), 1) <= downsample;
        combined_hpc = combined_hpc(filt,:);
        combined_pfc = combined_pfc(filt,:);
    end
end

%%
for a = 1:size(Behaviors, 1)
    for phase = 1:4
        for field  = string(fieldnames(pfc))'
            if ndims(pfc.(field)) ~= 5 || a > size(pfc.(field),1)
                continue
            end
            if all(pfc.(field)(a,:,phase,:,:,:) == 0, 'all')
                pfc.(field)(a,:,phase,:,:,:) = nan;
                hpc.(field)(a,:,phase,:,:,:) = nan;
            end
        end
    end
end

%% Just get a sense of where at all the partitions are
%%
pfc.org.sv = reshape(pfc.sv, [size(pfc.sv, 1:3), size(pfc.sv, 4)/K, K, size(pfc.sv, 5)]);
hpc.org.sv = reshape(hpc.sv, [size(hpc.sv, 1:3), size(hpc.sv, 4)/K, K, size(hpc.sv, 5)]);
pfc.org.m  = reshape(pfc.m,  [size(pfc.m,  1:3), size(pfc.m,  4)/K, K, size(pfc.m,  5)]);
hpc.org.m  = reshape(hpc.m,  [size(hpc.m,  1:3), size(hpc.m,  4)/K, K, size(hpc.m,  5)]);
pfc.org.am = reshape(pfc.am,  [size(pfc.am,  1:3), size(pfc.am,  4)/K, K, size(pfc.am,  5)]);
hpc.org.am = reshape(hpc.am,  [size(hpc.am,  1:3), size(hpc.am,  4)/K, K, size(hpc.am,  5)]);
%%
% Epoch-wise
close all
kws = {'savedir', savedir, 'patternSym', patternSym, 'patternNames', patternNames};
components.plotStatOverX(b, hpc, pfc, "sv", "epoch", kws{:});
components.plotStatOverX(b, hpc, pfc, "am", "epoch", kws{:});
components.plotStatOverX(b, hpc, pfc, "m", "epoch",  kws{:});

if animals == 1
    % State space performance wise
    components.plotStatOverX(b, hpc, pfc, "sv", "performance", kws{:});
    components.plotStatOverX(b, hpc, pfc, "am", "performance", kws{:});
end

g = findgroups();

% Tesing 1 2 3 4
%%
for a = 1:size(Behaviors, 1)
    for phase = 1:4
        for u = 1:numel(unique(groups))
        
            T = beh{a,phase}(groups == u, :).tperf_timewise;

            var_tperf{anim,phase}(u)         = var(bsxfun(@minus, T, median(T,1)),[],1); % var( fluctuation of component above its median )
            mu_tperf{anim,phase}(u)          = mean(T); % var( fluctuation of component above its median )
            abs_diff_tperf{anim,phase}(u)    = mean(abs(diff(T))); % var( fluctuation of component above its median )
            mu_diff_tperf{anim,phase}(u)     = abs(mean(diff(T))); % var( fluctuation of component above its median )
            % RY : I fixed this above : runs
            
            % RYAN : still need ci for this variance thing above.
            bootfun = @(x) mean(var(bsxfun(@minus, x, median(x,COMP)),[],TIME));
            var_ci{anim,phase}(u,:) = bootci(1000, bootfun, comp(:,groups == u));
            
        end
    end
end
%%
animalName_byBeh = repmat(animal, 4*2*Option.numPartition,1);
inxName_byBeh    = [repmat("hpc-hpc", 4*Option.numPartition,1);repmat("hpc-pfc", 4*Option.numPartition,1)];
behName          = repmat([repmat("reward", Option.numPartition,1);repmat("error", Option.numPartition,1);...
    repmat("inBoundChoice", Option.numPartition,1);repmat("outBoundChoice", Option.numPartition,1)],2,1);


temp_hpc = permute(slope_hpc_sv,[3,1,2]);
temp_hpc = reshape(temp_hpc, [200,12]);
temp_pfc = permute(slope_pfc_sv,[3,1,2]);
temp_pfc = reshape(temp_pfc, [200,12]);
SlopeTable_byBeh = array2table([animalName_byBeh, inxName_byBeh, behName, [temp_hpc; temp_pfc]]);
SlopeTable_byBeh.Properties.VariableNames(1:15) = ...
    {'animal','inx','beh','all','theta','delta','ripple','low-theta','low-delta','low-ripple','dim1','dim2','dim3','dim4','dim5'};
load('PolyFitTable_ByBeh.mat')
PolyFitTable_ByBeh = [PolyFitTable_ByBeh; SlopeTable_byBeh];
save('PolyFitTable_ByBeh','PolyFitTable_ByBeh','-v7.3')
%%
slope_hpc_sv_mean = mean(slope_hpc_sv, 3);
slope_pfc_sv_mean = mean(slope_pfc_sv, 3);
slope_hpc_sv_std  = std (slope_hpc_sv,0, 3);
slope_pfc_sv_std  = std (slope_pfc_sv,0, 3);
%%
fig('variance change fit across epochs, by behavior and pattern, hpc'); clf
for i = 1:4
    ax1 = nexttile;
    bar(ax1, slope_hpc_sv_mean(i, 1:4));
    hold on
end

fig('variance change fit across epochs, by behavior and pattern, pfc'); clf
for i = 1:4
    ax1 = nexttile;
    errorbar(ax1, slope_pfc_sv_mean(i, 1:4),slope_pfc_sv_std(i, 1:4));
    hold on
end

fig('variance change fit across epochs, by behavior and dim, hpc'); clf
for i = 1:3
    ax1 = nexttile;
    errorbar(ax1, slope_hpc_sv_mean(i, 8:10),slope_hpc_sv_std(i, 8:10));
    hold on
end

fig('variance change fit across epochs, by behavior and dim, pfc'); clf
for i = 1:3
    ax1 = nexttile;
    errorbar(ax1, slope_pfc_sv_mean(i, 8:10),slope_pfc_sv_std(i, 8:10));
    hold on
end


%% Plot by behavior and pattern
fig('variance of hpc comp strength across epochs, split by behavior and pattern '+animal);
clf; tiledlayout(2,2);
titles_byBeh = ["reward", "error", "inBoundChoice","outBoundChoice"];
titles_bypattern = ["all","theta","delta","ripple"];
for i = 1:4 % behavior
    
    nexttile
    for k = 1:4 % pattern
        for j = 1:numel(mean_hpc_sv_org{i, k})
            hold on
            stem(j:j,mean_hpc_sv_org{i, k, j});
            
        end
    end
    legend("all","theta","delta","ripple")
    title(titles_byBeh(i))
    hold on
    
    ylabel("variance")
    xlabel("epoch")
end

fig('variance of pfc comp strength across epochs, split by behavior and pattern '+animal);
clf; tiledlayout(2,2);
titles_byBeh = ["reward", "error", "inBoundChoice","outBoundChoice"];
titles_bypattern = ["all","theta","delta","ripple"];
for i = 1:4 % behavior
    
    nexttile
    for k = 1:4 % pattern
        for j = 1:numel(mean_pfc_sv_org{i, k})
            hold on
            stem(j:j,mean_pfc_sv_org{i, k, j});
            
        end
    end
    legend("all","theta","delta","ripple")
    title(titles_byBeh(i))
    hold on
    
    ylabel("variance")
    xlabel("epoch")
end
%% Plot by behavior and dimension
fig('variance of hpc comp strength across epochs, split by behavior and dimension '+animal);
subplot(2,1,1)
tiledlayout(2,2);
titles_byBeh = ["reward", "error", "inBoundChoice","outBoundChoice"];
titles_bypattern = ["all","theta","delta","ripple"];
for i = 1:4 % behavior
    
    nexttile
    for k = 8:10 % pattern
        for j = 1:numel(mean_hpc_sv_org{i, k})
            hold on
            stem(j:j,mean_hpc_sv_org{i, k, j});
            
        end
    end
    legend("dim 1","dim 2","dim 3")
    title(titles_byBeh(i))
    hold on
    
    ylabel("variance")
    xlabel("epoch")
end

fig('variance of pfc comp strength across epochs, split by behavior and dimension '+animal);
tiledlayout(2,2);
titles_byBeh = ["reward", "error", "inBoundChoice","outBoundChoice"];
titles_bypattern = ["all","theta","delta","ripple"];
for i = 1:4 % behavior
    
    nexttile
    for k = 8:10 % pattern
        for j = 1:numel(mean_pfc_sv_org{i, k})
            hold on
            stem(j:j,mean_pfc_sv_org{i, k, j}, 'color', "black", 'MarkerSize',(k/4)^2 );
            
        end
    end
    legend("dim 1","dim 2","dim 3")
    title(titles_byBeh(i))
    hold on
    
    ylabel("variance")
    xlabel("epoch")
end

%% look at comp strength mean

