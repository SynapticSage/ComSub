function out = getLogical(obj, index)
% getLogical  Get the logical value of a struct object and hold it's shape
%
%  take an arbitrary N-dim object and imagine we have an index that grabs
%  all or none of each of the N-dimensions ... let's return it holding
%  it's shape as much in tact as possible (rather than linearing/raveling)
%
%  say we have obj = [1 2 3; 4 5 6; 7 8 9] and we have a logical 
%  index of [1 0 1; 1 0 1; 1 0 1] ... we want to return [1 3; 4 6; 7 9]
%
% except obj can be any number of dimensions

% find which dimensions that index is truly grabbing all of
% and which dimensions it is grabbing none of
all_ = false(1, ndims(obj));
dim_selected = false(1, ndims(obj));
for i = 1:ndims(obj)
    all_(i) = mean(index, i);
    dim_selected(i) = all_  == 0 | all_ == 1;
end

getind = {}
for i = 1:ndims(obj)
end

