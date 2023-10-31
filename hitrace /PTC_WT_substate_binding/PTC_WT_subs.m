filenames = { ...
'raw_data', ...
};
ylimit = [800, 6000];
reorder = [];
[d_align, d_ref_align, ylimit, labels] = quick_look(filenames, ylimit, reorder);
%%
sequence = 'GGCCGTCAGCGAGTAGCTGACAACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACAAAAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGATGAAACAAACGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC';
sequence = strrep(sequence, 'T', 'U');
structure = '';
data_types = [ ...
repmat({'nomod'}, 1, 8),... 
repmat({'SHAPE'}, 1, 12), ...
repmat({'ddATP'}, 1, 1), ...
repmat({'ddTTP'}, 1, 1), ...
repmat({'ddCTP'}, 1, 1), ...
repmat({'ddGTP'}, 1, 1),...
];
offset = -23;
first_RT_nucleotide = length(sequence) - 20 + offset;
xsel = []; clf;
[xsel, seqpos, area_pred] = annotate_sequence(d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure);
savefig('TT022823_PTC_WT_bound_7SequenceAssign.fig')
%%
[area_peak, darea_peak] = fit_to_gaussians(d_align, xsel);
savefig('TT022823_PTC_WT_bound_8Annotate.fig')
%%
ref_segment = 'GAGUA';
ref_peak = get_ref_peak(sequence, ref_segment, offset);
sd_cutoff = 1.5;
saturated_idx = [9,10,11,12,13,14,15,16,17,18,19,20];
saturated_array_nomod = [mean(area_peak(:, 1:2), 2), mean(area_peak(:, 3:4), 2), mean(area_peak(:, 5:6), 2), mean(area_peak(:, 7:8), 2) ];
saturated_array = [saturated_array_nomod, area_peak(:, saturated_idx)];
saturated_error = [ ...
    mean(darea_peak(:, 1:2), 2), ...
    mean(darea_peak(:, 3:4), 2), ... 
    mean(darea_peak(:, 5:6), 2), ...
    mean(darea_peak(:, 7:8), 2), ...
    darea_peak(:, saturated_idx), ...
];
diluted_idx = saturated_idx;
diluted_array_nomod = saturated_array_nomod;
diluted_array = saturated_array;
diluted_error = saturated_error;
bkg_col = [1,2,3,4,1,1,1,2,2,2,3,3,3,4,4,4];
[normalized_reactivity, normalized_error, seqpos_out] = get_reactivities( ...
    saturated_array, diluted_array, saturated_error, diluted_error, ...
    bkg_col, ref_peak, seqpos, [], data_types([1,2,3,4, saturated_idx]), sequence, offset, sd_cutoff);
savefig('TT100923_PTC_WT_bound_9Normalized.fig')
%%
[d_SHAPE_0, da_SHAPE_0, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, [6,7]), normalized_error(:, [6,7]), [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_WT_bound_0uM_10SHAPE.fig')
%%
[d_SHAPE_50, da_SHAPE_50, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, [9,10]), normalized_error(:, [9,10]), [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_WT_bound_50uM_10SHAPE.fig')
%%
[d_SHAPE_100, da_SHAPE_100, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, [11,12]), normalized_error(:, [11,12]), [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_WT_bound_100uM_10SHAPE.fig')
%%
[d_SHAPE_200, da_SHAPE_200, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, [15,16]), normalized_error(:, [15,16]) , [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_WT_bound_200uM_10SHAPE_2reps.fig')
%%
plot(seqpos_out, d_SHAPE_0, 'k', 'linewidth', 2); hold on;
%plot(seqpos_out, d_SHAPE_50, 'r', 'linewidth', 2); hold on;
plot(seqpos_out, d_SHAPE_100, 'g', 'linewidth', 2); hold on;
%plot(seqpos_out, d_SHAPE_200, 'b', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_SHAPE_0, da_SHAPE_0, 'k'); hold on;
%errorbar(seqpos_out, d_SHAPE_50, da_SHAPE_50, 'r'); hold on;
errorbar(seqpos_out, d_SHAPE_100, da_SHAPE_50, 'g'); hold on;
%errorbar(seqpos_out, d_SHAPE_200, da_SHAPE_50, 'b'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 284:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-30 330 -0.5 10]);
legend('PTC WT');
title('1D SHAPE PTC WT', 'fontweight', 'bold', 'fontsize', 20);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out, first_RT_nucleotide + 1:20], 'xticklabel', sequence', 'fontsize', 5);
xtickangle(0)
%make_lines([0.5:50:380] - 0.5, 'k', 0.5, 0, 1);
%% no error bars
plot(seqpos_out(24:307), d_SHAPE_0(24:307), 'k', 'linewidth', 2); hold on;
plot(seqpos_out(24:307), d_SHAPE_50(24:307), 'r', 'linewidth', 2); hold on;
plot(seqpos_out(24:307), d_SHAPE_100(24:307), 'g', 'linewidth', 2); hold on;
plot(seqpos_out(24:307), d_SHAPE_200(24:307), 'b', 'linewidth', 2); hold on;
errorbar(seqpos_out(24:307), d_SHAPE_0(24:307), da_SHAPE_0(24:307), 'k'); hold on;
errorbar(seqpos_out(24:307), d_SHAPE_50(24:307), da_SHAPE_50(24:307), 'r'); hold on;
errorbar(seqpos_out(24:307), d_SHAPE_100(24:307), da_SHAPE_100(24:307), 'g'); hold on;
errorbar(seqpos_out(24:307), d_SHAPE_200(24:307), da_SHAPE_200(24:307), 'b'); hold on;
make_lines_horizontal(-0.5, 'k');
%make_lines([(seqpos_out(24:307))] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-5 290 -0.5 8]);
legend('PTC WT with 0 uM CCA-pcb and 0 uM C-pmn', 'PTC WT with 50 uM CCA-pcb and 50 uM C-pmn', 'PTC WT with 100 uM CCA-pcb and 100 uM C-pmn','PTC WT with 200 uM CCA-pcb and 200 uM C-pmn');
title('1D SHAPE PTC WT with CCA-pcb and C-pmn', 'fontweight', 'bold', 'fontsize', 50);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out(24:307)], 'xticklabel', sequence(24:307)', 'fontsize', 5);
xtickangle(0)
savefig('TT100923_PTC_WT_bound_summary.fig')
save('TT100923_PTC_WT_bound_workspace.mat')
