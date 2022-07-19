dt = 60;
spot_struct = main01_compile_traces('hbBAC-MS2-17_5C-Approved', 'dt', dt);

times_cell = {};
good_traces = [];
max_times = [];
for i = 1:length(spot_struct)
    if ~all(isnan(spot_struct(i).timeInterp))
        good_traces(end+1) = i;
        times_cell{end+1} = spot_struct(i).timeInterp;
        max_times(end+1) = max(times_cell{end});
    end
end
traces_array = zeros(max(max_times)/dt+1, length(good_traces));
time_vec = 0:dt:max(max_times);
for trace_index = 1:length(good_traces)
   traces_array(ismember(time_vec, times_cell{trace_index}), trace_index) = spot_struct(good_traces(trace_index)).fluoInterp;
end

[wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors, wt_ddd, ddd_boot_errors]  = weighted_autocorrelation(traces_array, 20,true, 100, ones(1,size(traces_array, 2)));

figure(1)
scatter(1:19, wt_dd)

figure(2)
scatter(1:18, wt_ddd)