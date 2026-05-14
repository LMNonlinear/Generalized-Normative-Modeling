
function idx=gen_loopgroup(num_y,num_y_loop,col_in_y,isflip)
if nargin<4 || isempty(isflip)
    isflip=false;
end
if nargin<3 || isempty(col_in_y)
    col_in_y=1;
end
assert(num_y_loop>=1,'each loop at least has one item')
num_y_loop=floor(num_y_loop);
% num_loop=ceil(num_loop);
% num_y_loop=floor(num_y/num_loop);
num_loop=ceil(num_y/num_y_loop);
num_y_loop=ceil(num_y/num_loop);

up=num_y_loop:num_y_loop:num_y;
down=1:num_y_loop:num_y;

if isempty(up)
    up(1)=num_y;
end 
if up(end)<num_y
    up(end+1)=num_y;
end

idx=cell(num_loop,1);
for i=1:num_loop
    idx{i}=down(i):up(i);
end
if isflip
    idx=flip(idx,1);
    idx=cellfun(@(x) flip(x),idx,'UniformOutput',false);
end

% idx=cellfun(@(x) kron(1:col_in_y,x),idx,'UniformOutput',false);
idx=cellfun(@(x) col_in_y*(min(x)-1)+1:col_in_y*max(x),idx,'UniformOutput',false);

end

