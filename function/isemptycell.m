function idx=isemptycell(A)
    idx=find(cellfun(@isempty,A));
end