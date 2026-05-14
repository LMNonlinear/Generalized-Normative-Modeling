function varargout=bit_convert(accuracy,varargin)
% accuracy lessthan 7 will turn to single
% large than 8 will turn to double
if isnumeric(accuracy)
    if accuracy<23
        varargout=cellfun(@(x) single(x), varargin,'UniformOutput',false);
    else %23<=acc<52
        varargout=cellfun(@(x) double(x), varargin,'UniformOutput',false);
    end
elseif ischar(accuracy)
    switch accuracy
        case 'single'
            varargout=cellfun(@(x) single(x), varargin,'UniformOutput',false);
        case 'double'
            varargout=cellfun(@(x) double(x), varargin,'UniformOutput',false);
    end
end

%%
% d = log2(eps('single'))
% d = log2(eps('double'))