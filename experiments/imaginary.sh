cd ..
make
sage trees_generation.sage -3 -7 -11 --gm-add=100 > d8_im_trees.log 2>&1
sage testrelations.sage > experiments/d8_im_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_im_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_im_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_im_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_im_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_im_dlog_5.log 2>&1

sage trees_generation.sage -3 -7 -11 -19 --gm-add=100 > d16_im_trees.log 2>&1
sage testrelations.sage > experiments/d16_im_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_im_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_im_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_im_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_im_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_im_dlog_5.log 2>&1

sage trees_generation.sage -3 -7 -11 -19 -23 --gm-add=100 > d32_im_trees.log 2>&1
sage testrelations.sage > experiments/d32_im_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_im_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_im_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_im_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_im_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_im_dlog_5.log 2>&1

sage trees_generation.sage -3 -7 -11 -19 -23 -31 --gm-add=100 > d64_im_trees.log 2>&1
sage testrelations.sage > experiments/d64_im_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_im_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_im_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_im_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_im_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_im_dlog_5.log 2>&1