function [bsize,scale] = bytesize( b, unit, type )
% 
%[bsize,scale] = bytesize( 1, 'b','int8')
% b: bytesize or number of element or data to measure
% unit: 'B' or 'b', out put use the same
% type: 'double', 'sinlge',et al. numerical type, only useful when b is scalar


% [bsize,scale] = osl_util.byte_size( b, unit='B' )
% Format input bytesize as a struct with conversions to kB, MB, GB, TB and PB
% If it is a BITsize (not BYTEsize), you can set the unit symbol to 'b' instead.
%
% If the input is not a scalar number, then we compute the bytesize of the argument
% itself using whos.
%
% NOTE: Due to limitations of Matlab's whos function, bytesize estimation for
% instances of handle classes is not accurate.
%
%
% The second output is the recommended unit to represent the input size.
% For example, the recommended unit for bytesize(123456) is 'kB'.
%
% If no output is collected, the size is printed to the console with the
% recommended unit.
%
% JH
if nargin < 3 || isempty(type), type=['']; end
if nargin < 2|| isempty(unit), unit='B'; end

if ~isnumeric(b) || ~isscalar(b)
    w = whos('b');
    b = w.bytes;
else % isnumeric(b) && isscalar(b)
    if ischar(type)
        switch type
            case {'uint64','int64','double'}
                b = 8*b;
            case {'uint32','int32','single'}
                b = 4*b;
            case {'uint16','int16'}
                b = 2*b;
            case {'uint8','int8'}
                
            otherwise
                error('not include')                
        end
    elseif isnumeric(type)
        w = whos('type');
        b = w.bytes*b;
    end
end

if strcmp(unit,'b')
    b=b*8;
end


assert( isnumeric(b) && isscalar(b), 'Bad input bytesize.' );

units = mapfun( @(x) [x unit], {'','k','M','G','T','P'}, false );
scale = units{min( numel(units), 1+floor( log(b)/log(1024) ) )};

for i = 1:numel(units)
    u = units{i};
    bsize.(u) = b / 1024^(i-1);
end

if nargout == 0
    fprintf('%.2f %s\n', bsize.(scale), scale );
end

end
