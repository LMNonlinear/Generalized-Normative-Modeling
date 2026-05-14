function g=set_gscatter_para(T,group,low_d)
if nargin<3|| isempty(low_d)
    low_d=2;
end
low_d=2;
t= unique(get_subtable(T,[],{'name',group}));
groupset=unique(t.(group));
% load('.\result\batch_id.mat');
% ind=find_char(batch_id,groupset);
% color=plt.batch14;
color=parula(length(groupset));%turbo

% gcol=zeros(size(g,1),3);
% g.col=repmat(zero(length(g,1)))
shape=repmat('.',length(groupset),1);



if low_d==2
   g.color= color;%color(ind,:);
   g.shape=shape;    
else
    t=initializeTableVar(t,'color','cell');
    t=initializeTableVar(t,'shape','cell');
    for k=1:length(groupset)
        if isnumeric(groupset)
            idx_group=find(t.(group),groupset(k));
        else
            idx_group=find_char(t.(group),groupset{k});
        end
        
        g.color(idx_group)={color(k,:)};
        g.shape(idx_group)={shape(k)};
    end
end

g.group=t.(group);


end