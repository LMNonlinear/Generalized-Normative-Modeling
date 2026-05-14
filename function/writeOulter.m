function     t=writeOulter(name,OutIdx,tag)
con_name=unique(name);%batch{1}
t=table(name);
t.(tag)=false(length(name),1);
% T=initializeTableVar(T,title,{'logical'});
idx_title=find_char(t.Properties.VariableNames,tag);
idx_outer=find_char(name,con_name(OutIdx));
t(idx_outer,idx_title)=array2table(true(length(idx_outer),1));
disp(con_name(OutIdx));
