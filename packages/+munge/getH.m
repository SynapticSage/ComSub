function [Data, Option] = getH(Data, Option, genh)
    genH = nd.fieldGet(Option, 'generateH');
    genH = replace(genH, "fromCoherence  fromRipTimes", "coh");
    genH = replace(genH, "fromSpectra  fromRipTimes", "spec");
    genH = replace(genH, "fromWpli  fromRipTimes", "wpli");
    genH = replace(genH, "coherence", "coh");
    genH = replace(genH, "spectra", "spec");
    indz = bsxfun(@times, genH==genh, ones(size(Data)));
    % Data = Data(genH == genh, :, :, :, :, :, :);
    origsz = size(Data);
    Data = Data(logical(indz));
    if ~isempty(Data)
        Data = reshape(Data, [1 origsz(2:end)]);
        Option = Option(genH == genh);
    else
        Data = [];
        Option = [];
    end
end
