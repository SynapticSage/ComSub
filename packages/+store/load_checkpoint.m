function varargout = load_checkpoint(hash, varargin)

    tic
    if ~endsWith(hash, ".mat")
	hash = hash + ".mat";
    end
    m = matfile(hash);
    % m = matfile("bef0923.mat", "Writable", true);
    % Events             = util.matfile.getdefault(m, 'Events', []);
    % Spk                = util.matfile.getdefault(m, 'Spk', []);
    % Patterns           = util.matfile.getdefault(m, 'Patterns', []);
    % Patterns_overall   = util.matfile.getdefault(m, 'Patterns_overall', []);
    % Components         = util.matfile.getdefault(m, 'Components', []);
    % Components_overall = util.matfile.getdefault(m, 'Components_overall', []);
    % Option             = util.matfile.getdefault(m, 'Option', []);
    varargout = cell(numel(varargin), 1);
    disp("Loaded variables: ")
    cnt = 0;
    for v = varargin(:)'
	cnt = cnt+1;
	disp(v{1})
	varargout{cnt} = util.matfile.getdefault(m, v{1}, []);
    end
    disp("...done in " + num2str(toc) + " seconds")

end
