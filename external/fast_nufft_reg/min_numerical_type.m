function [type,byte]=min_numerical_type(minval,maxval,redundancy,isunsign,isfloating)
% [type,byte]=min_numerical_type(0,1000);
if nargin<1 || isempty(minval)
    minval = 0;
end
if nargin<2 || isempty(maxval)
    error('maxval is required')
end
if nargin<3 || isempty(redundancy)
    redundancy = 1;
end
if nargin<4 || isempty(isunsign)
    isunsign = minval>=0;
end
if nargin<5 || isempty(isfloating)
    isfloating = ~all(isintegernum([minval,maxval]));
end
unint_types={'uint8'; 'uint16'; 'uint32'; 'uint64'};
int_types={'int8'; 'int16'; 'int32'; 'int64';};
float_types={'single'; 'double';};
if isunsign && ~isfloating
    types=[unint_types;int_types;float_types];
    types_max_val=cellfun(@(x) double(intmax(x)),[unint_types;int_types]);
    types_max_val=[types_max_val;cellfun(@(x) double(realmax(x)),[float_types])];
elseif ~isunsign && ~isfloating
    types=[int_types;float_types];
    types_max_val=cellfun(@(x) double(intmax(x)),[int_types]);
    types_max_val=[types_max_val;cellfun(@(x) double(realmax(x)),[float_types])];
elseif isfloating
    types=[float_types];
    types_max_val=[cellfun(@(x) double(realmax(x)),[float_types])];
end
idx_avalible=(types_max_val*redundancy-maxval)>0;
if isempty(idx_avalible)
    error('overflow for your maxval')
end
[~,idx_min]=min(types_max_val(idx_avalible));
type=types(idx_avalible);
type=type{idx_min};
byte=bytesize(1,'B',type);
byte=byte.B;
end