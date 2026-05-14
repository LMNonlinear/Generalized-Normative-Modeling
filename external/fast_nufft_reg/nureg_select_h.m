function [regs]=nureg_select_h(regs)
% if ~regs.opt.y_corr_bandwith
if ischar(regs.opt.select_h)
    if strcmp(regs.opt.select_h,'id1se')
        %% in hinput
        regs.h=regs.h(regs.gcv_yhat.id1se,:);
        
        %% in the h(ihgood)
        ihselect=false(regs.dh_input,1);
        ihselect(regs.gcv_yhat.id1se)=true;
        ihgood=(~regs.ihbad);
        idhselect=find(ihselect(ihgood));
        regs.kdf=regs.kdf(regs.ndcolon{:},idhselect);
        regs.s=regs.s(regs.ndcolon{:},idhselect);
        
        regs.pdof_inv_m=regs.pdof_inv_m(idhselect);
        regs.pdof_inv_sd=regs.pdof_inv_sd(idhselect);
        regs.pdof_m=regs.pdof_m(idhselect);
        regs.pdof_sd=regs.pdof_sd(idhselect);
        
        regs.gcv_yhat.gcv_lo=regs.gcv_yhat.gcv_lo(idhselect);
        regs.gcv_yhat.gcv_m=regs.gcv_yhat.gcv_m(idhselect);
        regs.gcv_yhat.gcv_sd=regs.gcv_yhat.gcv_sd(idhselect);
        regs.gcv_yhat.gcv_up=regs.gcv_yhat.gcv_up(idhselect);
        
        regs.ihbad=regs.ihbad(regs.gcv_yhat.id1se);
        regs.gcv_yhat.id1se=1;
        regs.gcv_yhat.idmin=1;
        regs.dh=1;
        regs.dh_input=1;
    elseif strcmp(regs.opt.select_h,'idmin')
        %% in hinput
        regs.h=regs.h(regs.gcv_yhat.idmin,:);
        
        %% in the h(ihgood)
        ihselect=false(regs.dh_input,1);
        ihselect(regs.gcv_yhat.idmin)=true;
        ihgood=(~regs.ihbad);
        idhselect=find(ihselect(ihgood));
        regs.kdf=regs.kdf(regs.ndcolon{:},idhselect);
        regs.s=regs.s(regs.ndcolon{:},idhselect);
        
        regs.pdof_inv_m=regs.pdof_inv_m(idhselect);
        regs.pdof_inv_sd=regs.pdof_inv_sd(idhselect);
        regs.pdof_m=regs.pdof_m(idhselect);
        regs.pdof_sd=regs.pdof_sd(idhselect);
        if ~isempty(regs.gcv_yhat)
            regs.gcv_yhat.gcv_lo=regs.gcv_yhat.gcv_lo(idhselect);
            regs.gcv_yhat.gcv_m=regs.gcv_yhat.gcv_m(idhselect);
            regs.gcv_yhat.gcv_sd=regs.gcv_yhat.gcv_sd(idhselect);
            regs.gcv_yhat.gcv_up=regs.gcv_yhat.gcv_up(idhselect);
        end
        regs.ihbad=regs.ihbad(regs.gcv_yhat.idmin);
        regs.gcv_yhat.id1se=1;
        regs.gcv_yhat.idmin=1;
        regs.dh=1;
        regs.dh_input=1;
    end
elseif isnumeric(regs.opt.select_h) && ~isempty(regs.opt.select_h) %isvector(regs.opt.select_h) ||~isempty(regs.opt.select_h)
    %% in hinput
    idx=ismember(regs.h,regs.opt.select_h,'row');
    
    id=find(idx);
    regs.h=regs.h(id,:);
    
    %% in the h(ihgood)
    ihselect=false(regs.dh_input,1);
    ihselect(id)=true;
    ihgood=(~regs.ihbad);
    idhselect=find(ihselect(ihgood));%only in good bandwidth
    regs.kdf=regs.kdf(regs.ndcolon{:},idhselect);
    regs.s=regs.s(regs.ndcolon{:},idhselect);
    
    regs.pdof_inv_m=regs.pdof_inv_m(idhselect);
    regs.pdof_inv_sd=regs.pdof_inv_sd(idhselect);
    regs.pdof_m=regs.pdof_m(idhselect);
    regs.pdof_sd=regs.pdof_sd(idhselect);
    regs.df_m=regs.df_m(idhselect);
    regs.df_sd=regs.df_sd(idhselect);
    if ~isempty(regs.gcv_yhat)
        regs.gcv_yhat.gcv_lo=regs.gcv_yhat.gcv_lo(idhselect);
        regs.gcv_yhat.gcv_m=regs.gcv_yhat.gcv_m(idhselect);
        regs.gcv_yhat.gcv_sd=regs.gcv_yhat.gcv_sd(idhselect);
        regs.gcv_yhat.gcv_up=regs.gcv_yhat.gcv_up(idhselect);
    end
    regs.ihbad=regs.ihbad(id);
    idxoutfale=false(length(id),1);
    idx1se=idxoutfale;
    idxmin=idxoutfale;
    if ~isempty(regs.gcv_yhat)
        if regs.gcv_yhat.id1se<length(id)
            idx1se(regs.gcv_yhat.id1se)=true;
        end
        
        if regs.gcv_yhat.idmin<length(id)
            idxmin(regs.gcv_yhat.idmin)=true;
        end
        regs.gcv_yhat.id1se=find(idx1se);%may empty
        regs.gcv_yhat.idmin=find(idxmin);
    end
    regs.dh=length(id);
    regs.dh_input=length(id);
    
end

%     regs.h1se= regs.h(regs.gcv_yhat.id1se,:);
%     regs.hmin= regs.h(regs.gcv_yhat.idmin,:);
% end
end