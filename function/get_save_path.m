function path=get_save_path(path_input,addon,tag,ext)
% addon is just make sure path will be difference with path_input, 
% if tag exist, addon will not use.
[filepath,filename,extRaw] = fileparts(path_input);
if isempty(tag)
    filename=[filename,addon];
end
filename=add_tag(filename,tag);
if nargin<4 || isempty(ext)
    ext=extRaw;
end
path=[filepath,filesep,filename,ext];

end