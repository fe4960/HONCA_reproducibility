cell=args[1]
dirInName=args[2] #"major_flt_cpm05_lib_fix_anc_correct_all_ancFix"  #args[2]
fn=args[3] #"/storage/chentemp/wangj/human_ret_anc/data/HRCA_obs_sample_corrected.txt"
dirIn=args[5]
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
res=read.table(paste0(cell,"_DEG_age"),sep="\t",header=T)

#logFC   AveExpr t       P.Value adj.P.Val       B       z.std   qval
#HES4    -0.0610147419753712     2.0949373039422 -8.31653480893321 
