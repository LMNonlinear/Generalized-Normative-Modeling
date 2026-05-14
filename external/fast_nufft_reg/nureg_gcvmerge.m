function [gcvout]=nureg_gcvmerge(gcvin)
if isfield(gcvin,'yhat')
    dimy=ndims(gcvin(1).yhat);
elseif isfield(gcvin,'yhat_1se')
    dimy=ndims(gcvin(1).yhat_1se);
end
fields_all=fieldnames(gcvin);
fields_fpp={'fpp'};
fields_uni={'metric'};% later may consider the metric for gcv also different for each y
fields_array=setdiff(fields_all,[fields_fpp(:);fields_uni(:)]);

%% save unique fields
for ifield=1:length(fields_uni)
    if isfield(gcvin,fields_uni{ifield})
        gcvout.(fields_uni{ifield})=gcvin.(fields_uni{ifield});
    end
end
%% concatenate array type data
for ifield=1:length(fields_array)
    gcvout.(fields_array{ifield})=cat(dimy,gcvin.(fields_array{ifield}));
end
%% concatenate fpp
num_loop=length(gcvin);
gcvout.(fields_fpp{1})=gcvin(1).(fields_fpp{1});
dx=length(gcvin(1).(fields_fpp{1}).GridVectors);
dy=0;
num_y_loop=zeros(num_loop,1);
for iloop=1:num_loop
    dyi=size(gcvin(iloop).(fields_fpp{1}).Values);
    num_y_loop(iloop)=prod(dyi(dx+1:end));
    dy=dy+num_y_loop(iloop);
end

N=length(gcvin(1).(fields_fpp{1}).GridVectors{1});
ndcolon(1:dx) = {':'};
gcvout.(fields_fpp{1}).Values=zeros([N*ones(1,dx),dy]);
ky=0;
for iloop=1:num_loop
    iy=ky+[1:num_y_loop(iloop)];
    gcvout.(fields_fpp{1}).Values(ndcolon{:},iy)=gcvin(iloop).(fields_fpp{1}).Values;
    ky=num_y_loop(iloop)+ky;
end
end