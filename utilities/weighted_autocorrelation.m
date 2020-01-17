%Create Weighted Average Autocorrelation
function [wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors, wt_ddd, ddd_boot_errors] = ...
    weighted_autocorrelation(traces, lags, bootstrap,n_boots,trace_weights)
    %traces: array of traces with zeros preceeding and succeeding period of
    %        activity. Oriented column-wise
    %lags: num lags to use
    %bootstrap: binary var. If 1, returns array of bootstrap errors
    auto_array = zeros(lags+1,size(traces,2));
    time_steps = zeros(lags+1,size(traces,2));
    %Convert NaNs to zeros
    traces(isnan(traces)) = 0;
    t_filter = sum(traces>0)>1;
    traces = traces(:,t_filter);
    trace_weights = trace_weights(t_filter);
    if ~bootstrap
        n_boots = 1;
    end
    samples = NaN(lags+1,n_boots);    
    dd_samples = NaN(lags-1,n_boots);
    ddd_samples = NaN(lags-2,n_boots);
    for b = 1:n_boots
        %If bootstrap errors are desired, take n_boots samples with
        %replacement, each the size of original array
        if bootstrap
            s_vec = 1:size(traces,2);
            s = randsample(s_vec,length(s_vec),true,trace_weights);
            trace_sample = traces(:,s);
        else
            trace_sample = traces;
        end        
        for col = 1:size(trace_sample,2)            
            trace = trace_sample(:,col);
            if sum(isnan(trace)) > 0 
                error('Problem with NaN filtering');
            end            
            trace_active = trace(find(trace,1):find(trace,1,'last'));
            if length(trace_active) < lags + 1
%                 warning('Length of input trace insufficient for specified number of lags')
                t_lags = length(trace_active) - 1;
            else
                t_lags = lags;
            end       
            xcv = xcov(trace_active,t_lags,'Normalized');
            auto_array(1:t_lags+1,col) = xcv(t_lags+1:end);
            time_steps(1:t_lags+1,col) = fliplr((length(trace_active)-t_lags):length(trace_active));
        end
        %Take weighted mean. Traces with length < lags should just not be
        %factored in for lags beyond their length
        numerator = sum(auto_array.*time_steps,2);
        denominator = sum(time_steps,2);
        samples(:,b) = numerator ./ denominator;
        dd_samples(:,b) = diff(diff(samples(:,b)));
        ddd_samples(:,b) = diff(samples(:,b),3);
    end
    wt_autocorr = mean(samples,2);
    wt_dd = mean(dd_samples,2);
    wt_ddd = mean(ddd_samples,2);
    a_boot_errors = std(samples')';
    dd_boot_errors = std(dd_samples')';
    ddd_boot_errors = std(ddd_samples')';