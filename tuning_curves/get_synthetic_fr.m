function [fr] = get_synthetic_fr(cellid, fr,ft,model_mean,norm_type);
% 1. Flat firing rate
% 2. firing rate = model mean
% 3. step function encoding of model
%
normflag = strcmp(norm_type, 'div');
fr_in = fr;
if cellid == -1
    fr = ones(size(fr_in))+randn(size(fr_in));
    fr(isnan(fr_in)) = NaN;
elseif cellid == -11
    fr = 10.*ones(size(fr_in))+randn(size(fr_in));
    fr(isnan(fr_in)) = NaN;
elseif cellid < -1
    dex = ft > 0 & ft < length(model_mean)*1e-3;
    fr = ones(size(fr_in))+randn(size(fr_in));
    dt = min(diff(ft));
    T = length(model_mean)*1e-3;
    model_dex = 0:dt:T;
    model_dex = ceil(model_dex/1e-3) +1;
    model_dex(model_dex > length(model_mean)) = [];
    if length(model_dex) < sum(dex);
        model_dex = [model_dex length(model_mean)];
    end
    model_samples = model_mean(model_dex);
    
    if cellid == -2
        model_samples = model_samples + 10;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        if normflag;
        fr = fr ./10;
        end
    elseif cellid == -20
        model_samples = model_samples + 10;
        model_samples(model_samples < 0) = 0;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -21
        model_samples = model_samples + 5;
        model_samples(model_samples < 0) = 0;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -22
        model_samples = model_samples + 20;
        model_samples(model_samples < 0) = 0;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -23
        model_samples = model_samples + 10;
        model_samples(model_samples < 0) = 0;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        fr = fr./5;
    elseif cellid == -24
        model_samples = model_samples;
        model_samples(model_samples < 0) = 0;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        fr = fr.*5;
    elseif cellid == -3;
        model_samples(model_samples > 0 ) = 1.5;
        model_samples(model_samples < 0 ) = 0.5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -30;
        model_samples(model_samples > 0 ) = 2.5;
        model_samples(model_samples < 0 ) = 0.5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -31;
        model_samples(model_samples > 0 ) = 4.5;
        model_samples(model_samples < 0 ) = 2.5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -32;
        model_samples(model_samples > 0 ) = 10.5;
        model_samples(model_samples < 0 ) = 5.5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
    elseif cellid == -4;
        model_samples(model_samples > 0 ) = 10;
        model_samples(model_samples < 0 ) = 5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        if normflag;
        fr = fr./7.5;
        end
    elseif cellid == -5;
        fr = ones(size(fr_in))+randn(size(fr_in));
        fr = fr + (1:length(ft))/50;
        fr(isnan(fr_in)) = NaN;
        if normflag;
        fr = fr ./ mean((1:length(ft))/50);
        end
    elseif cellid == -6
        model_samples(model_samples > 0 ) = 10;
        model_samples(model_samples < 0 ) = 5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        if normflag;
        fr = fr./7.5;
        end
    elseif cellid == -7
        model_samples(model_samples > 0 ) = 10;
        model_samples(model_samples < 0 ) = 5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        if normflag;
        fr = fr./7.5;
        end
    elseif cellid == -8
        model_samples(model_samples > 0 ) = 10;
        model_samples(model_samples < 0 ) = 5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
        if normflag;
        fr = fr./7.5;
        end
    elseif cellid == -9
        model_samples(model_samples > 0 ) = 10;
        model_samples(model_samples < 0 ) = 5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
         if normflag;       
        fr = fr./7.5;
        end
    elseif cellid == -10
        model_samples(model_samples > 0 ) = 10;
        model_samples(model_samples < 0 ) = 5;
        fr(dex) = model_samples;
        fr(isnan(fr_in)) = NaN;
                if normflag;
        fr = fr./7.5;
        end
    end
else
    error('Synthetic neuron ID not found')
end




