function T = gethash(T, hash)
% acquire the row(s) matching the hash

T = T(T.hash == hash,:);

end
