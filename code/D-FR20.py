#########################D-FR20############################
#########################D-FR20###########################
import cmapPy
import sys
import pandas as pd
from cmapPy.pandasGEXpress.parse_gctx import parse
import gseapy as gp
df = pd.read_csv("fer_20_entrez_id.txt", sep="\t",header=None)
df.iloc[:,0]
up=df.iloc[0:8,1].tolist()
down=df.iloc[8:20,1].tolist()
up_list = [str(i) for i in up]
down_list = [str(i) for i in down]
geneset={'up' : up_list,'down': down_list}
tanzhen = parse("level5_beta_trt_cp_n720216x12328.gctx",cidx=[0]) #you can download from L1000


def GSEA(start,end):
    outname="cmap/output_all/"+str(start)+"_"+str(end)+".txt"
    for i in range(int(start), int(end)):
        tanzhen = parse("level5_beta_trt_cp_n720216x12328.gctx",cidx=[i])
        pre_res = gp.prerank(rnk=tanzhen.data_df, # or rnk = rnk,
                     gene_sets=geneset,
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )
        out=pre_res.res2d
        out['cmap_cellline_drug']=[tanzhen.data_df.columns.tolist(),tanzhen.data_df.columns.tolist()]
        out.to_csv(outname, sep="\t", index=False,header=False,mode="a")

if __name__=='__main__':
    GSEA(sys.argv[1],sys.argv[2])