### GEM_01_Data_Loading.ipynd 
Script to load the relevant data into jupyter notebook and creates one single GEM input file ('gem_input_data.pkl') to be used in the model. 
'gem_input_data.pkl' contains: 
- metadata (18×16)
- nmr (18×38)
- microbiome (18×44)
- pathways (18×250)
- mean_expr (23017×11) (single-cell)
- mean_expr_high (23017×11) (single-cell high SA >50)
- mean_expr_low (23017×11) (single-cell low SA <20)
- metabolite_map (dict) (GEM model IDs for metabolites)
- gem_metabolites (list of 10)
- common_ids (list of 18) (IDs that match between all tier 1 data)
