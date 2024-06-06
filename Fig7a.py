import cmapPy
import sys
import pandas as pd
#help(cmapPy.pandasGEXpress)
from cmapPy.pandasGEXpress.parse_gctx import parse
import gseapy as gp
from gseapy.plot import gseaplot
df = pd.read_csv("fer_20_entrez_id.txt", sep="\t",header=None)
df.iloc[:,0]
up=df.iloc[0:8,1].tolist()
down=df.iloc[8:20,1].tolist()
up_list = [str(i) for i in up]
down_list = [str(i) for i in down]
geneset={'up' : up_list,'down': down_list}
data = parse("level5_beta_trt_cp_n720216x12328.gctx",ridx=[1])
data1=data.data_df
colname=data1.columns.tolist()
pattern="LCP001_MCF10A.WT_24H:B08"
match_index =[colname.index(pattern)]
tanzhen = parse("level5_beta_trt_cp_n720216x12328.gctx",cidx=[match_index[0]])
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
terms = pre_res.res2d.Term
tanzhen.data_df.to_csv("filgotinib_foldchange.txt",index=True)
gseaplot(rank_metric=pre_res.ranking, term=terms[0], ofname='fig7a-1.pdf', **pre_res.results[terms[0]])
gseaplot(rank_metric=pre_res.ranking, term=terms[1], ofname='fig7a-2.pdf', **pre_res.results[terms[1]])