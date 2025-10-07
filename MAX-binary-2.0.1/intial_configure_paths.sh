#!/bin/bash -l
##### Date:04 June 2020
##### configure.sh: set the current path for loading R functions
istr="/cfs/klemming/projects/snic/snic2020-6-4/Meghana/MAX-binary-2.0.1"
ostr="/path/to"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/AEM_update_X_beta.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/Create_count_matrix.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/buildCRP.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/mergeMutSingleSample.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/estimateBeta_MAX2.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/splitFasta.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' new_constructMutXmatrix.sh"
eval "sed -i -e 's#"$istr"#"$ostr"#g' runMAX_step1_Xmat.sh"
eval "sed -i -e 's#"$istr"#"$ostr"#g' runMAX_step2_alignment.sh"
eval "sed -i -e 's#"$istr"#"$ostr"#g' runMAX_step3_createYcount.sh"
eval "sed -i -e 's#"$istr"#"$ostr"#g' runMAX_step4_quant.sh"

echo "Done!"
