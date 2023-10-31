filenames = { ...
'combined_0216_0227', ...
};
ylimit = [800, 6000];
reorder = [];
[d_align, d_ref_align, ylimit, labels] = quick_look(filenames, ylimit, reorder);
%%
sequence = 'GGCCGTCAGCGAGTAGCTGACAATCCGCGGCAAGACGGAAAGACCCCGTGCACCTTTACTATAACTTGACACAAAAGTGTCTGGTAGGTAGTTTGACTGGGATGGTCACGGATAAAAGGTACGCCGGGGATAACAGGCTAATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATTTCGGCTCATCGCATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGCGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGAATGAAACAAACGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC';
sequence = strrep(sequence, 'T', 'U');
structure = '';
data_types = [ ...
repmat({'nomod'}, 1, 12),... 
repmat({'SHAPE'}, 1, 20), ...
repmat({'ddATP'}, 1, 1), ...
repmat({'ddTTP'}, 1, 1), ...
repmat({'ddCTP'}, 1, 1), ...
repmat({'ddGTP'}, 1, 1),...
];
offset = -23;
first_RT_nucleotide = length(sequence) - 20 + offset;
%xsel = []; clf;
[xsel, seqpos, area_pred] = annotate_sequence(d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure);
savefig('TT100923_PTC_1-1_bound_7SequenceAssign.fig')
%%
[area_peak, darea_peak] = fit_to_gaussians(d_align, xsel);
savefig('TT100923_PTC_1-1_bound_8Annotate.fig')
%%
ref_segment = 'GAGUA';
ref_peak = get_ref_peak(sequence, ref_segment, offset);
sd_cutoff = 1.5;
saturated_idx = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32];
saturated_array_nomod = [mean(area_peak(:, 1:3), 2), mean(area_peak(:, 4:6), 2), mean(area_peak(:, 7:9), 2), mean(area_peak(:, 10:12), 2) ];
saturated_array = [saturated_array_nomod, area_peak(:, saturated_idx)];
saturated_error = [ ...
    mean(darea_peak(:, 1:3), 2), ...
    mean(darea_peak(:, 4:6), 2), ... 
    mean(darea_peak(:, 7:9), 2), ...
    mean(darea_peak(:, 10:12), 2), ...
    darea_peak(:, saturated_idx), ...
];
diluted_idx = saturated_idx;
diluted_array_nomod = saturated_array_nomod;
diluted_array = saturated_array;
diluted_error = saturated_error;
bkg_col = [1, 2, 3, 4, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4];
[normalized_reactivity, normalized_error, seqpos_out] = get_reactivities( ...
    saturated_array, diluted_array, saturated_error, diluted_error, ...
    bkg_col, ref_peak, seqpos, [], data_types([1,2,3,4, saturated_idx]), sequence, offset, sd_cutoff);
savefig('TT100923_PTC_1-1_bound_9Normalized.fig')
%%
[d_SHAPE_0, da_SHAPE_0, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 5:8), normalized_error(:, 5:8), [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_1-1_bound_0uM_10SHAPE.fig')
%%
[d_SHAPE_50, da_SHAPE_50, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, [10,11,12,14]), normalized_error(:, [10,11,12,14]), [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_1-1_bound_50uM_10SHAPE.fig')
%%
[d_SHAPE_100, da_SHAPE_100, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 16:19), normalized_error(:, 16:19), [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_1-1_bound_100uM_10SHAPE.fig')
%%
[d_SHAPE_200, da_SHAPE_200, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 20:23), normalized_error(:, 20:23) , [], seqpos_out, sequence, offset); 
savefig('TT100923_PTC_1-1_bound_200uM_10SHAPE_2reps.fig')
%%
plot(seqpos_out, d_SHAPE_0, 'k', 'linewidth', 2); hold on;
%plot(seqpos_out, d_SHAPE_50, 'r', 'linewidth', 2); hold on;
%plot(seqpos_out, d_SHAPE_100, 'g', 'linewidth', 2); hold on;
plot(seqpos_out, d_SHAPE_200, 'b', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_SHAPE_0, da_SHAPE_0, 'k'); hold on;
%errorbar(seqpos_out, d_SHAPE_100, da_SHAPE_50, 'g'); hold on;
errorbar(seqpos_out, d_SHAPE_200, da_SHAPE_50, 'b'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 284:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-30 330 -0.5 10]);
legend('PTC 1.1');
title('1D SHAPE PTC 1.1', 'fontweight', 'bold', 'fontsize', 20);
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
axis([-10 290 -10 10]);
legend('PTC 1.1 with 0 uM CCA-pcb and 0 uM C-pmn', 'PTC 1.1 with 50 uM CCA-pcb and 50 uM C-pmn', 'PTC 1.1 with 100 uM CCA-pcb and 100 uM C-pmn','PTC 1.1 with 200 uM CCA-pcb and 200 uM C-pmn');
title('1D SHAPE PTC 1.1 with CCA-pcb and C-pmn', 'fontweight', 'bold', 'fontsize', 50);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out(24:307)], 'xticklabel', sequence(24:307)', 'fontsize', 5);
xtickangle(0)
savefig('TT100923_PTC1.1_bound_summary.fig')
save('TT100923_PTC_1.1_bound_workspace.mat')
