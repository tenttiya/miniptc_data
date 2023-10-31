filenames = { ...
'TT_111521_PTC_1-3', ...
};
ylimit = [800, 6000];
reorder = [];
[d_align, d_ref_align, ylimit, labels] = quick_look(filenames, ylimit, reorder);
%%
sequence = 'GGCCGTCAGCGAGTAGCTGACAACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACAAAAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGTTTTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGATCTCTATCTGCCGTGGGATGAAACAAACGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC';
sequence = strrep(sequence, 'T', 'U');
structure = '';
data_types = [ ...
  repmat({'nomod'}, 1, 2), repmat({'SHAPE'}, 1, 3), repmat({'DMS'}, 1, 3), repmat({'ddATP'}, 1, 1), repmat({'ddTTP'}, 1, 1), ...
  repmat({'ddCTP'}, 1, 1), repmat({'ddGTP'}, 1, 1)
];
offset = 20;
first_RT_nucleotide = length(sequence) - 20 + offset;
%xsel = []; clf;
[xsel, seqpos, area_pred] = annotate_sequence(d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure);
savefig('TT111521_PTC_1-3_7SequenceAssign.fig')
%%
[area_peak, darea_peak] = fit_to_gaussians(d_align, xsel);
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
%%
[d_SHAPE, da_SHAPE, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 2:4), normalized_error(:, 2:4), [], seqpos_out, sequence, offset); 
%%
[d_DMS, da_DMS, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 5:7), normalized_error(:, 5:7), [], seqpos_out, sequence, offset); 