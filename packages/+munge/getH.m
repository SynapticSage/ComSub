function [Data, Option] = getH(Data, Option, genh)
    genH = nd.fieldGet(Option, 'generateH');
    genH = replace(genH, "fromCoherence  fromRipTimes", "coh");
    genH = replace(genH, "fromSpectra  fromRipTimes", "spec");
    genH = replace(genH, "fromWpli  fromRipTimes", "wpli");
    indz = bsxfun(@times, genH==genh, ones(size(Data)));
    % Data = Data(genH == genh, :, :, :, :, :, :);
    origsz = size(Data);
    Data = Data(logical(indz));
    Data = reshape(Data, [1 origsz(2:end)]);
    Option = Option(genH == genh);
end
