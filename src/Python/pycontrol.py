
# make sure garbage has been cleaned out
import gc 
gc.collect()

def smilei_check():
    """Do checks over the script"""

    #count all the Smilei instances
    n_smilei=0
    for key,value in globals().items() :
        if isinstance(value,Smilei) :
            n_smilei += 1

    assert (n_smilei == 1),"Only one Smilei instance allowed : "+str(n_smilei)
    
    sim=get_smilei()
    
    
    ## do some checks on "sim" here
    
    

smilei_check()


