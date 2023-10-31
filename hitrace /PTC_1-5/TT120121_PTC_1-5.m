filenames = { ...
'TT_111521_PTC_1-5', ...
};
ylimit = [800, 6000];
reorder = [];
[d_align, d_ref_align, ylimit, labels] = quick_look(filenames, ylimit, reorder);
%%
sequence = 'GGCCGTCAGCGAGTAGCTGACAACCCGCGGCAAGACGGAAAGACCCCGTGCACCTTTACTATAACTTGATACAAAAGTGTCTGGTAGGTAGTTTGACTGGGGCGGTCACGGATAAAAGGTACGCCGGGGATAACAGGCTAATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATTTCGGCTCATCGCATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGCGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGATGAAACAAACGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC';
sequence = strrep(sequence, 'T', 'U');
structure = '';
data_types = [ ...
  repmat({'nomod'}, 1, 2), repmat({'SHAPE'}, 1, 3), repmat({'DMS'}, 1, 3), repmat({'ddATP'}, 1, 1), repmat({'ddTTP'}, 1, 1), ...
  repmat({'ddCTP'}, 1, 1), repmat({'ddGTP'}, 1, 1)
];
offset = -23;
first_RT_nucleotide = length(sequence) - 20 + offset;
xsel = []; clf;
[xsel, seqpos, area_pred] = annotate_sequence(d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure);
savefig('TT120121_PTC_1-5_7SequenceAssign.fig')
%%
[area_peak, darea_peak] = fit_to_gaussians(d_align, xsel);
savefig('TT120121_PTC_1-5_8Annotate.fig')
%%
ref_segment = 'GAGUA';
ref_peak = get_ref_peak(sequence, ref_segment, offset);
sd_cutoff = 1.5;
saturated_idx = [3, 4, 5, 6, 7, 8];
saturated_array_nomod = [mean(area_peak(:, 1:2), 2)];
saturated_array = [saturated_array_nomod, area_peak(:, saturated_idx)];
saturated_error = [ ...
    mean(darea_peak(:, 1:2), 2), ...
    darea_peak(:, saturated_idx), ...
];
diluted_idx = saturated_idx;
diluted_array_nomod = saturated_array_nomod;
diluted_array = saturated_array;
diluted_error = saturated_error;
bkg_col = [1, 1, 1, 1, 1, 1, 1];
[normalized_reactivity, normalized_error, seqpos_out] = get_reactivities( ...
    saturated_array, diluted_array, saturated_error, diluted_error, ...
    bkg_col, ref_peak, seqpos, [], data_types([1, saturated_idx]), sequence, offset, sd_cutoff);
savefig('TT120121_PTC_1-5_9Normalized.fig')
%%
[d_SHAPE, da_SHAPE, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 2:4), normalized_error(:, 2:4), [], seqpos_out, sequence, offset); 
savefig('TT120121_PTC_1-5_10SHAPE.fig')
%%
[d_DMS, da_DMS, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:,5:7), normalized_error(:, 5:7), [], seqpos_out, sequence, offset);
savefig('TT120121_PTC_1-5_11DMS.fig')
%%
d_DMS_minus_seq = d_DMS_minus(24:307);
d_SHAPE_minus_seq = d_SHAPE_minus(24:307);
%%
plot(seqpos_out, d_SHAPE, 'b', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_SHAPE, da_SHAPE, 'b'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 284:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-30 330 -0.5 8]);
legend('PTC 1.5');
title('1D SHAPE PTC 1.5', 'fontweight', 'bold', 'fontsize', 20);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out, first_RT_nucleotide + 1:20], 'xticklabel', sequence', 'fontsize', 5);
xtickangle(0)
%make_lines([0.5:50:380] - 0.5, 'k', 0.5, 0, 1);
savefig('TT120121_PTC_1-5_SHAPE_combined.fig')
%%
plot(seqpos_out, d_DMS, 'r', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_DMS, da_SHAPE, 'r'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 284:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-30 330 -2 15]);
legend('PTC 1.5');
title('1D DMS PTC 1.5', 'fontweight', 'bold', 'fontsize', 20);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out, first_RT_nucleotide + 1:20], 'xticklabel', sequence', 'fontsize', 5);
xtickangle(0)
savefig('TT120121_PTC_1-5_DMS_combined.fig')
