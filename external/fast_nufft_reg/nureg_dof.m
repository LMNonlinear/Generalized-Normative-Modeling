function [regs]=nureg_dof(regs)
%(1/n)w'*(w-what)== (1/n)tr(I-A)
if regs.opt.calc_dof
    p=randn(regs.Ty,regs.opt.num_dof_sample);
    
    [mq]=nureg_reg(regs,p);%N*ones(1,regs.dx) regs.dh regs.opt.num_dof_sample
    
    mq=reshape(mq,[regs.N*ones(1,regs.dx),regs.dh*regs.opt.num_dof_sample]);
    method ='griddedInterpolant';
    opt.Method='spline';
    opt.ExtrapolationMethod ='linear';
    [~,phat]=nureg_gridinterp(regs,mq,[],method,opt);
    phat=reshape(phat,[regs.Ty,regs.dh,regs.opt.num_dof_sample]);
    phat=permute(phat,[1,3,2]);
    
    pdof=(sum(p.*(p-phat))./sum(p.*p));
    pdof=permute(pdof,[2,3,1]);
    
    df=regs.Ty-regs.Ty*pdof;
    if sum(df<0,[1,2])>0
        df(df<=0)=0;
        warning('df  is negative');
    end
    
    regs.pdof_inv_sd=std(1./pdof.^2,[],1).';
    regs.pdof_inv_m=mean(1./pdof.^2,1).';
    regs.df_m=mean(df,1).';
    regs.df_sd=std(df,[],1).';
    regs.pdof_m=mean(pdof,1).';
    regs.pdof_sd=std(pdof,[],1).';
else
    regs.pdof_inv_sd=[];
    regs.pdof_inv_m=[];
    regs.df_m=[];
    regs.df_sd=[];
    regs.pdof_m=[];
    regs.pdof_sd=[];
end

end
