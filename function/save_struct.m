function path_struct=save_struct(s,path_struct,addon,tag)
path_struct=get_save_path(path_struct,addon,{tag},'.mat');
save(path_struct,'-struct','s','-v7.3')
end
