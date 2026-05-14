function s = set_defaults(s, field, value)
% set_defaults  Set a default value for a struct field if it does not exist.
%
%   s = set_defaults(s, field, value)
%
%   If s.(field) does not exist, set it to value. Otherwise keep existing.

if ~isfield(s, field)
    s.(field) = value;
end
end
