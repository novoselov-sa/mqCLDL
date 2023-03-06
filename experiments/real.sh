cd ..
make
sage trees_generation.sage 5 13 17 --gm-add=100 > experiments/d8_re_trees.log 2>&1
sage testrelations.sage > experiments/d8_re_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_re_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_re_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_re_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_re_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d8_re_dlog_5.log 2>&1

sage trees_generation.sage 5 13 17 29 --gm-add=100 > experiments/d16_re_trees.log 2>&1
sage testrelations.sage > experiments/d16_re_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_re_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_re_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_re_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_re_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d16_re_dlog_5.log 2>&1

sage trees_generation.sage 5 13 17 29 37 --gm-add=100 > experiments/d32_re_trees.log 2>&1
sage testrelations.sage > experiments/d32_re_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_re_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_re_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_re_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_re_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d32_re_dlog_5.log 2>&1

sage trees_generation.sage 5 13 17 29 37 41 --gm-add=100 > experiments/d64_re_trees.log 2>&1
sage testrelations.sage > experiments/d64_re_rels.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_re_dlog_1.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_re_dlog_2.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_re_dlog_3.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_re_dlog_4.log 2>&1
sage dlog_cyc_PoC.sage > experiments/d64_re_dlog_5.log 2>&1