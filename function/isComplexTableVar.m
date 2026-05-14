function type=isComplexTableVar(T,vars)
if nargin<2 || isempty(vars)
    %     vars=T.Properties.Variablenames;
else
    T=get_subtable(T,[],vars);
end
type=varfun(@iscomplex, T, 'OutputFormat', 'uniform');



end