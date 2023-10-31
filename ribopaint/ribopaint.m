%% Load .mat from workspace and residue locations to set up
%change this to other miniPTC TIF file
mat_name = 'TT102023_1.1.mat';
imagex = imread('miniPTC_1.1_v2.tif');
square_width = 80;
%change this to other miniPTC sequence
%sequence = 'CCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACAAAAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGG';
seqpos = 1:284;
%residue_locations = [];
offset = 0;
whichres = seqpos - offset;
%residue_locations = pick_points(imagex_2, offset, residue_locations, square_width);
%%
color_scheme = 1;
whattoplot = d_SHAPE(24:307);
color_profile = color_palette(whattoplot, 1, -1, color_scheme);

%%
is_circle = 1;
%color_residues(imagex, residue_locations, whichres, whattoplot, color_profile, square_width, is_circle);
color_circles(imagex, residue_locations, whichres, whattoplot, color_profile, square_width, 'miniPTC-1.1_diff');
save(mat_name)
