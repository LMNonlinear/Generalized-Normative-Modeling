function [regs]=nureg_y(regs)

% prepareing y
% if tensor
if regs.opt.y_corr_bandwith
    dhy=size(regs.y,2);
    if dhy==1
        y=repmat(regs.y,[1,regs.dh,1]);
    else
        y=regs.y;
    end
    y=reshape(y,[regs.Ty,regs.dh*regs.dy]);
    dy=regs.dh*regs.dy;
else
    y=regs.y;
    dy=regs.dy;
end

% if variance
if strcmp(regs.opt.y_type_out,'variance')
    yraw=y;
    y=log(y+1/regs.Ty);
end

%Todo!:
%in case memory decreasing with loops, the number of regression in each
%loop better be decreasing.
if isunix
    redundancy=4;
else
    redundancy=6;
end
dtype_nufft=min_numerical_type(0,2*max(regs.L));
if regs.opt.nufft_loop
    byte_nufft_xindx=bytesize(((2*regs.Msp+1)^regs.dx)*regs.Tx*(regs.dx),'B',dtype_nufft);
else
    byte_nufft_xindx=bytesize(((2*regs.Msp+1)^regs.dx)*regs.Tx*dy*(regs.dx+1),'B',dtype_nufft);
end
byte_nufft_data=bytesize(max(((2*regs.Msp+1)^regs.dx)*regs.Tx*dy,2*2*((2*max(regs.L))^regs.dx)*dy),'B',regs.dtype);%ysp and Ftau
byte_nufft=byte_nufft_xindx.B+byte_nufft_data.B;
byte_reg=bytesize((regs.L)^regs.dx*regs.dy*regs.dh,'B',dtype_nufft);
byte_reg=byte_reg.B;
bsize = max(byte_nufft,byte_reg)*redundancy;

if ~regs.opt.isgpu
    [ispara,freemem]=test_memory(bsize,true); %#ok<ASGLU>
else
    [ispara,freemem]=test_gpu(bsize,1); %#ok<ASGLU>
    regs.s=gpuArray(regs.s);
    regs.kdf=gpuArray(regs.kdf);
end
if ~regs.opt.y_corr_bandwith
    idx=gen_loopgroup(dy,freemem/(bsize/dy));
else
    idx=gen_loopgroup(regs.dy,freemem/(bsize/regs.dy),regs.dh);
end
num_loop=numel(idx);


%% loop for y regression
if ~regs.opt.y_corr_bandwith && regs.opt.calc_dof
    for iloop=num_loop:-1:1
        iy=idx{iloop};
        %% y regression (important part)
        % fft y
        y_ft=nureg_nufft(regs,y(:,iy));
        % regression
        [mq]=nureg_reg(regs,y_ft);
        if regs.opt.y_keep_all
            % grid interpolation
            [gcv(1).fpp,gcv(1).yhat]=nureg_gridinterp(regs,mq,regs.xlist,regs.opt.y_grid_method,regs.opt.y_grid_opt);
            gcv(2)=nureg_gcvmerge(gcv);
            if iloop==1
                gcv=gcv(2);
            end
        else
            % grid interpolation
            [fpp,yhat]=nureg_gridinterp(regs,mq,regs.xlist,regs.opt.y_grid_method,regs.opt.y_grid_opt);
            % gcv
            [gcv(iloop)]=nureg_gcv(regs,y(:,iy),yhat,fpp,regs.opt.y_gcv_metric,regs.opt.y_gcv_dstd);            
        end
    end
    if regs.opt.y_keep_all
        regs.gcv_yhat=[];
        regs.fpp_yhat=gcv.fpp;
        regs.yhat=gcv.yhat;
    else
        regs.gcv_yhat=nureg_gcvmerge(gcv);
        regs.fpp_yhat=regs.gcv_yhat.fpp;
        regs.gcv_yhat.fpp=[];
        regs.yhat=regs.gcv_yhat.yhat_1se;
        regs.gcv_yhat.yhat_1se=[];
    end
else
    for iloop=num_loop:-1:1
        iy=idx{iloop};
        % fft y
        if ~regs.opt.isgpu
            y_ft=nureg_nufft(regs,y(:,iy));
        else
            y_ft=gather(nureg_nufft(regs,gpuArray(y(:,iy))));
        end
        % regression
        [mq]=nureg_reg(regs,y_ft);
        % grid interpolation
        [gcv(1).fpp,gcv(1).yhat]=nureg_gridinterp(regs,mq,regs.xlist,regs.opt.y_grid_method,regs.opt.y_grid_opt);
        gcv(2)=nureg_gcvmerge(gcv);
        if iloop==1
            gcv=gcv(2);
        end
    end
    regs.gcv_yhat=[];
    regs.fpp_yhat=gcv.fpp;
    regs.yhat=gcv.yhat;
end
if regs.opt.y_corr_bandwith
    regs.yhat=reshape(regs.yhat,[regs.Ty,regs.dh,regs.dy]);
end
%%
regs.d=[];
if strcmp(regs.opt.y_type_out,'variance')
    % regs.yhat shape depends on the code path taken above:
    %   - gcv+calc_dof path (line 61-93): [Tx, dy] — already band-selected
    %   - loop path (y_corr_bandwith or !calc_dof, line 94-115 + reshape 117):
    %                                     [Tx, dh, dy] after reshape
    % Thus we normalize to [Tx, dh_eff, dy_eff] using regs.yhat's actual shape.
    yhat_sz  = size(regs.yhat);
    if numel(yhat_sz) < 3
        yhat_sz(end+1:3) = 1;  % pad trailing singletons for indexing
    end
    dh_eff = yhat_sz(2);
    dy_eff = yhat_sz(3);

    % Align yraw to [Ty, dh_eff, dy_eff]:
    if regs.opt.y_corr_bandwith
        yraw = reshape(yraw, [regs.Ty, regs.dh, regs.dy]);
    else
        % yraw is [Ty, dy]; broadcast over dh_eff (usually 1 here anyway)
        yraw = reshape(yraw, [regs.Ty, 1, regs.dy]);
    end

    d = 1./((1/regs.Ty).*sum(yraw.*exp(-regs.yhat),1));  % [1, dh_eff, dy_eff]
    regs.yhat = exp(regs.yhat) ./ d;

    % fpp_yhat.Values layout mirrors regs.yhat's trailing dims
    n_out = dh_eff * dy_eff;
    d_flat = reshape(d, [1, n_out]);
    n_grid = regs.N^regs.dx;
    n_out_fpp = numel(regs.fpp_yhat.Values) / n_grid;
    if n_out_fpp ~= n_out
        % In the gcv path fpp may still contain all bands (dh*dy original
        % cv bands, but dy_eff has been collapsed to dy). Use fpp's own n_out.
        n_out = n_out_fpp;
        % expand d across bands if needed (cyclic/broadcast)
        if mod(n_out, dh_eff*dy_eff) == 0
            rep = n_out / (dh_eff * dy_eff);
            d_flat = repmat(d_flat, 1, rep);
        end
    end
    fpp_V = reshape(regs.fpp_yhat.Values, [regs.N * ones(1, regs.dx), n_out]);
    fpp_V = exp(fpp_V) ./ reshape(d_flat, [ones(1, regs.dx), n_out]);
    regs.fpp_yhat.Values = reshape(fpp_V, size(regs.fpp_yhat.Values));
    regs.d = d;
end



if any(isinf(regs.yhat(:)))|| any(isnan(regs.yhat(:)))
    error('yhat have bad value, try enlarge opt.Nratio')
end

%% Clamp fpp_yhat grid values to prevent spline oscillation artifacts
%  fpp grid values can overflow between sparse data points (Runge's phenomenon).
%  Clamp to the range of actual data-point values with a safety margin.
%  Layout notes:
%    fpp_yhat.Values : [N,...,N, n_out_fpp] where first dx axes are grid;
%                      n_out_fpp depends on code path:
%                        - gcv calc_dof path: n_out_fpp = dh * dy (all bands)
%                        - loop path        : n_out_fpp = dy
%                      Trailing singleton may be squeezed by MATLAB.
%    regs.yhat       : [Tx, n_out_y] with n_out_y = dy in the gcv path
%                      (already band-selected) or varies in the loop path.
%  To be robust we derive n_out_fpp from numel(Values)/N^dx, and loop over
%  yhat columns cyclically (for the gcv path, where yhat has fewer columns
%  than fpp_vals — each yhat column corresponds to dh fpp-value columns).
n_grid   = regs.N^regs.dx;
fpp_sz   = size(regs.fpp_yhat.Values);
n_out_fpp = numel(regs.fpp_yhat.Values) / n_grid;
if mod(n_out_fpp,1) ~= 0
    error('fpp_yhat.Values size %s inconsistent with grid %d', mat2str(fpp_sz), n_grid);
end
fpp_vals = reshape(regs.fpp_yhat.Values, [regs.N * ones(1, regs.dx), n_out_fpp]);
yhat_ref = reshape(regs.yhat, regs.Tx, []);
n_out_y  = size(yhat_ref, 2);
margin   = 3;  % allow fpp to extend 3x beyond data range
ndcolon  = repmat({':'}, 1, regs.dx);
for ir = 1:n_out_fpp
    v   = fpp_vals(ndcolon{:}, ir);
    y_r = yhat_ref(:, mod(ir-1, n_out_y) + 1);
    y_range = max(y_r(:)) - min(y_r(:));
    y_range = max(y_range, eps);
    lo = min(y_r(:)) - margin * y_range;
    hi = max(y_r(:)) + margin * y_range;
    v_clamped = max(min(v, hi), lo);
    % Replace Inf/NaN with median
    bad = isinf(v_clamped) | isnan(v_clamped);
    if any(bad(:))
        v_clamped(bad) = median(v_clamped(~bad), 'all');
    end
    fpp_vals(ndcolon{:}, ir) = v_clamped;
end
% Restore the original .Values shape (MATLAB may squeeze trailing singleton)
fpp_vals_out = reshape(fpp_vals, fpp_sz);
if any(fpp_vals_out(:) ~= regs.fpp_yhat.Values(:))
    n_fixed = sum(fpp_vals_out(:) ~= regs.fpp_yhat.Values(:));
    warning('nureg_y:fpp_clamp', ...
        'Clamped %d/%d fpp_yhat grid values (spline oscillation).', ...
        n_fixed, numel(fpp_vals_out));
    regs.fpp_yhat.Values = fpp_vals_out;
end
