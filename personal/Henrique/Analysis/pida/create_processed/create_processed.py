from tqdm import tqdm
import polars as pl
import sys
import os
user_name = os.getlogin()
sys.path.append(f'/home/{user_name}/Documents/Atmos_Pandora/apc_atmo_repo/Anatree/')
sys.path.append(f'/home/henrique/Documents/Atmos_Pandora/apc_atmo_repo/personal/Henrique/Analysis/pida/')
from  ana_tools import *
from pida_functions import *
from itertools import product
import argparse


def write_complex(dfcomplex_en, basename="numu", lepname="mu", cheatlep=False, cheatpi=False, cheatpr=False):
    complexname = "complex"
    if cheatlep or cheatpi or cheatpr:
        complexname+="_cheat"
    if cheatlep and lepname!="":
        complexname+=f"_{lepname}"
    if cheatpi:
        complexname+="_pi"
    if cheatpr:
        complexname+="_pr"
    dfcomplex_en.write_parquet(f"../data/processed/{basename}/{complexname}.parquet")

    return
def process_data_nue(dfall:pl.DataFrame, cheatlep = False, cheatpi = False, cheatpr = False, minPNC=0.0):
    dfsimple = create_ecandidate(dfall, clear_low_en_as_shower=False, filter_protons=False)
    dfsimple_en = simple_energy_nue(dfsimple)
    dfe = create_ecandidate(dfall, cheat=cheatlep)
    dfprselected = create_proton_candidate(dfe, cheat=cheatpi, cheatMinPNC=minPNC, return_only_selected=False)
    dfprpiselected = create_pion_candidate(dfprselected, cheat=cheatpr,
                                           cheatMinPNC=minPNC,
                                           lower_threshold_len=20,
                                           return_only_selected=False)
    dfcomplex_en = complex_energy_nue(dfprpiselected, W="W")
    dfsimple_en = dfsimple_en.select(
        dfcomplex_en.drop('Kpi', 'Kpr', 'OtherPFPs', 'mass_to_add').columns
    )
    if not (cheatlep or cheatpi or cheatpr):
        dfsimple_en.write_parquet("../data/processed/nue/simple.parquet")
    write_complex(dfcomplex_en, basename="nue", lepname="e", cheatlep=cheatlep, cheatpi=cheatpi, cheatpr=cheatpr)

def process_data_nc(dfall:pl.DataFrame, cheatlep = False, cheatpi = False, cheatpr = False, minPNC=0.0):
    dfsimple_en = simple_energy_nc(dfall)
    dfe = dfall
    dfprselected = create_proton_candidate(dfe, cheat=cheatpi, cheatMinPNC=minPNC, return_only_selected=False)
    dfprpiselected = create_pion_candidate(dfprselected, cheat=cheatpr, cheatMinPNC=minPNC, lower_threshold_len=20, return_only_selected=False)
    dfcomplex_en = complex_energy_nc(dfprpiselected, W="W")

    dfsimple_en = dfsimple_en.select(
        dfcomplex_en.drop('Kpi', 'Kpr', 'OtherPFPs', 'mass_to_add').columns
    )
    if not (cheatlep or cheatpi or cheatpr):
        dfsimple_en.write_parquet("../data/processed/nc/simple.parquet")
    write_complex(dfcomplex_en, basename="nc", lepname="", cheatlep=cheatlep, cheatpi=cheatpi, cheatpr=cheatpr)

def process_data_numu(dfall:pl.DataFrame, cheatlep = False, cheatpi = False, cheatpr = False, minPNC=0.0):
    dfsimple = create_mucandidates(dfall)
    dfsimple_en = simple_energy(dfsimple)
    if not cheatlep:
        dfmus = create_mucandidates_2(dfall, fineselection=False, afterfine=False)
    else:
        dfmus = create_mucandidates(dfall, cheat=True)
    dfprselected = create_proton_candidate(dfmus, cheat=cheatpi, cheatMinPNC=minPNC)
    dfprselected = join_pi_candidate(dfprselected, dfmus, 'pr')
    dfprpiselected = create_pion_candidate(dfprselected, cheat=cheatpr, cheatMinPNC=minPNC)
    dfprpiselected = join_pi_candidate(dfprpiselected, dfprselected, 'pi')
    dfcomplex_en = complex_energy2(dfprpiselected, W="W", force_bigger_pr=False, do_not_force_bigger_pi=False)

    dfsimple_en = dfsimple_en.select(
        dfcomplex_en.drop('Kpi', 'Kpr', 'OtherPFPs', 'mass_to_add').columns
    )

    if not (cheatlep or cheatpi or cheatpr):
        dfsimple_en.write_parquet("../data/processed/numu/simple.parquet")
    write_complex(dfcomplex_en, basename="numu", lepname="mu", cheatlep=cheatlep, cheatpi=cheatpi, cheatpr=cheatpr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-minPNC", "--minPNC", type=float, default=0.0)
    parser.add_argument("-n", "--neutrinos", type=str, default="numu")
    args = vars(parser.parse_args())
    neutrinos = args["neutrinos"]
    minPNC=0.
    df = pl.scan_parquet("../data/pida_visen.parquet") 
    dftrkg4 = getbestof(df)

    
    if neutrinos == "nc":
        dfall = dftrkg4.filter(
            pl.col('ccnc_truth') == 1,
            )
    else:
        dfall = dftrkg4.filter(
            pl.col('ccnc_truth') == 0,
            )

    if neutrinos == "numu":
        dfall = dfall.filter(
           pl.col('nuPDG_truth').abs() == 14
           )
    elif neutrinos == "nue":
        dfall = dfall.filter(
           pl.col('nuPDG_truth').abs() == 12
           )
        
    if isinstance(dfall, pl.LazyFrame):
        dfall = dfall.collect()

    dfall = dfall.sort('trklen').with_columns(
           pnc = pl.col('trkpurity_planes_B')*pl.col('trkcompleteness_planes_B'),
           )

    possibles = [True, False]
    process_data = process_data_nue
    if neutrinos == "numu":
        process_data = process_data_numu
    elif neutrinos == "nue":
        process_data = process_data_nue
    elif neutrinos == "nc":
        process_data = process_data_nc

    repeat = 3 if neutrinos != "nc" else 2
    res = product(possibles, repeat=repeat)
    res = ([ False, False, False ],)
    # res = ([ False, False ],)

    for docheats in tqdm(res):
        c_lep, c_pi, c_pr = (docheats if neutrinos != "nc" else (False, *docheats))
        print(c_lep, c_pi, c_pr)

        process_data(dfall, cheatlep=c_lep, cheatpi=c_pi, cheatpr=c_pr, minPNC=minPNC)

