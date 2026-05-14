function [T,tregs]=cv_tnureg_zmap(T,y,tregs,opt)
if nargin<4 || isempty(opt)
    opt=[];
end

if ~isfield(opt,'maxibf')
    opt.maxibf=1;
end
if ~isfield(opt,'global_h')
    opt.global_h=false;
end
if isstruct(tregs)
    tregs={tregs};
end
dy=numel(y);
n=size(T,1);

%% standardize
if opt.standardize
    % muy=mean(reshape(table2array(get_subtable(T,[],y)),n,[],dy),1);
    % sigmay=var(reshape(table2array(get_subtable(T,[],y)),n,[],dy),1);
    func=@(x) [mean(x,1);std(real(x),1)+1i*std(imag(x))];
    mu=func(reshape(table2array(get_subtable(T,[],y)),n,[],dy));
    stdy=mu(2,:,:);%/2warning('enlarge variance')
    mu=mu(1,:,:);
    for i=1:dy
        T.(y{i})=T.(y{i})-mu(:,:,i);
        T.(y{i})=real(T.(y{i}))./real(stdy(:,:,i))+1i*imag(T.(y{i}))./imag(stdy(:,:,i));
    end
end
%% initialize resp to regress
num_level=size(tregs,1);
for i=1:num_level
    tregs{i,1}=tnureg_create(tregs{i,1});
end
%% regression

for ibf=1:opt.maxibf
    for i=1:num_level
        disp(['runing level: ',num2str(i)])
        %% mu (and musigma)
        if i==1 && ibf==1 % first layer and first time, assign response=y
            T=clone_table_var(T,tregs{i,1}.resp,y);
        elseif i==1 && ibf>1 % first layer and >1 time(ibf is not necessary, remove in the future)
            if tregs{num_level,1}.opt.calc_z
                T=clone_table_var(T,tregs{i,1}.resp,tregs{num_level,1}.opt.resp_z);
            else
                T=clone_table_var(T,tregs{i,1}.resp,tregs{num_level,1}.opt.resp_eps);
            end
        else % >1 layer, assign response=z
            if tregs{i-1,1}.opt.calc_z
                T=clone_table_var(T,tregs{i,1}.resp,tregs{i-1,1}.opt.resp_z);
            else
                T=clone_table_var(T,tregs{i,1}.resp,tregs{i-1,1}.opt.resp_eps);
            end
            %             tregs{i,1}.opt.inherit_bandwidth=cell2mat(reshape(arrayfun(@(x) find(ismember(tregs{1}.hList{1}(:,find_char(tregs{i-1}.factor,tregs{i}.factor)),tregs{2}.hList{1}(x,:),'row')),1:size(tregs{2}.hList{1},1),'UniformOutput',false),1,[])).';
            if isempty(tregs{i}.factor_mu)
                tregs{i,1}.opt.inherit_bandwidth=[];
            elseif ~opt.y_keep_all
                % When y_keep_all=false, the previous level's output is
                % already GCV-selected (single column). No bandwidth
                % inheritance needed — each row maps to column 1.
                dh_child = size(tregs{i}.hList{1}, 1);
                tregs{i,1}.opt.inherit_bandwidth = ones(dh_child, 1);
            else
                tregs{i,1}.opt.inherit_bandwidth=cell2mat(reshape(arrayfun(@(x) find(ismember(tregs{i-1}.hList{1}(:,find_char(tregs{i-1}.factor_mu,tregs{i}.factor_mu)),tregs{i}.hList{1}(x,:),'row')),1:size(tregs{i}.hList{1},1),'UniformOutput',false),1,[])).';
            end
            %             warning('test 2nd level')
            %             dh=size(tregs{i,1}.hList{1, 1},1);
            %             for iresp=1:num_resp
            %                 for ih=1:dh
            %                     T.(tregs{i,1}.resp{iresp})(:,ih)=T.(tregs{i,1}.resp{iresp})(:,1);
            %                 end
            %             end
        end
        [T,tregs{i,1}]=cv_tnureg(T,tregs{i,1});
    end
end
%% check tregs
check_tregs(T,tregs);

%% compact to save memory

if opt.compact.y
    T=asnarray2table(T,[],y,[]);
end
if opt.compact.sigma
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_sigma,[]);
    end
end
%% get ystar
% y=mg+sg*ms+sg*ss*mb+sg*ss*sb*z
% ystar=mg+sg*ms+0+sg*ss*1*z

% y=mg+sg(ms+ss(mb+sb*z))
% ystar=mg+sg(ms+ss(0+1*z))

ystar=append('ystar',y);% T=asnarray2table(T,[],ystar,1);
tregs{1,1}.opt.ystar=ystar;

if tregs{num_level}.opt.calc_musigma
    T=clone_table_var(T,ystar,tregs{num_level,1}.opt.resp_z);
else
    T=clone_table_var(T,ystar,tregs{num_level,1}.opt.resp_eps);
end

for i=num_level:-1:1
    if tregs{i}.opt.calc_musigma
        if ~opt.isremove(i)
            func=@(mu,sigma,ystar) (real(mu)+real(sigma).*real(ystar))+...
                1i*(imag(mu)+imag(sigma).*imag(ystar));
            T = tablefun(func,T,ystar,[],...
                tregs{i}.opt.resp_mu,...
                tregs{i}.opt.resp_musigma,...
                ystar);
        end
    else
        if ~opt.isremove(i)
            func=@(mu,ystar) real(mu)+real(ystar)+...
                1i*(imag(mu)+imag(ystar));
            T = tablefun(func,T,ystar,[],...
                tregs{i}.opt.resp_mu,...
                ystar);
        end
    end
end

%% unstandardize for ystar
% for i=1:dy
%     T.(ystar{i})=T.(ystar{i})./sqrt(sigmay(:,:,i));
% end
if opt.unstandardize
    for i=1:dy
        T.(ystar{i})=real(T.(ystar{i})).*real(stdy(:,:,i))+1i*imag(T.(ystar{i})).*imag(stdy(:,:,i));
        T.(ystar{i})=T.(ystar{i})+mu(:,:,i);
    end
end
%% compact to save memory
if opt.compact.mu
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_mu,[]);
    end
end

if opt.compact.eps
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_eps,[]);
    end
end
if opt.compact.musigma
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_musigma,[]);
    end
end

if opt.compact.z
    for i=1:num_level
        T=asnarray2table(T,[],tregs{i,1}.opt.resp_z,[]);
    end
end
if opt.compact.ystar
    for i=1:num_level
        T=asnarray2table(T,[],ystar,[]);
    end
end


end

%%
% a=T.(tregs{i,1}.resp{3});
% x(:,1)=T.(tregs{1, 1}.factor{1});
% x(:,2)=T.(tregs{1, 1}.factor{2});
% scatter3(x(:,1),x(:,2),real(a(:,1)))
%











