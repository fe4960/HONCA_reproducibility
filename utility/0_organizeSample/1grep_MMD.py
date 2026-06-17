import pandas as pd
import numpy as np
dr="HCA/ON_ONH/data/sample_list"
lst=["ON_ONH_ATAC","ON_ONH"]
for lst1 in lst:
    mmd=pd.read_csv(f"{dr}/MMD_donor_ID",sep="\t",header=0)
    fl=pd.read_csv(f"{dr}/{lst1}.txt",sep="\t",header=0)
#fl=pd.read_csv(f"{dr}/ON_ONH.txt",sep="\t",header=0)
    fl=fl[fl["sampleid"].str.contains("MMD")]
    fl["donorid1"]=fl["sampleid"].str.replace("3v31_MMD_","").str.replace("_ON_bead","").str.replace("MMD_","").str.replace("_ONH","").str.replace("ATAC-MMD-","").str.replace("-ON-bead","").str.replace("-ONH","").str.replace("-ON","").str.replace("-","_")
    mmd["Tube Label"]=mmd["Tube Label"].str.replace("-","_")
    mmd=mmd.set_index("Tube Label")
    fl["donorid2"]=fl["donorid1"]
    fl=fl.set_index("sampleid")
    unk=[]
    for id1 in fl.index:
        if fl.loc[id1,"donorid1"] in mmd.index:
            fl.loc[id1,"donorid2"]=np.unique(mmd.loc[fl.loc[id1,"donorid1"],"Donor ID"])
        else:
            unk.append(id1)
    fl.to_csv(f"{dr}/{lst1}_corrected.txt",sep="\t",header=True)
    pd.DataFrame(unk).to_csv(f"{dr}/{lst1}_uncorrected.txt",sep="\t",header=True)

#id1=mmd.index.intersection(fl["donorid1"])
#fl["donorid2"]=mmd.loc[fl["donorid1"],"Donor ID"]
#fl["donorid2"]=mmd.loc[id1,"Donor ID"]
