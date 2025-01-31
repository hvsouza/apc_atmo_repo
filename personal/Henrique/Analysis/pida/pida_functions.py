import polars as pl
import sys
import os
import particle
# using getlogin() returning username
user_name = os.getlogin()

sys.path.append(f'/home/{user_name}/Documents/Atmos_Pandora/apc_atmo_repo/Anatree/')
from  ana_tools import *


def getbestof(df:pl.DataFrame ) -> pl.DataFrame:
    df:pl.DataFrame
    dftmp = df.group_by(selection_events(), maintain_order=True).agg(
        pl.col('^trkcalo_.*$').sum(),
        pl.col('^trknhits_.*$').sum(),
        pl.col('^allcalo_.*$').last(),
        pl.col('^nallhits_.*$').last(),
    ).with_columns(
        otherallcalo_planes_U = pl.col('allcalo_planes_U') - pl.col('trkcalo_planes_U'),
        otherallcalo_planes_V = pl.col('allcalo_planes_V') - pl.col('trkcalo_planes_V'),
        otherallcalo_planes_W = pl.col('allcalo_planes_W') - pl.col('trkcalo_planes_W'),
        othernallhits_planes_U = pl.col('nallhits_planes_U')-pl.col('trknhits_planes_U'),
        othernallhits_planes_V = pl.col('nallhits_planes_V')-pl.col('trknhits_planes_V'),
        othernallhits_planes_W = pl.col('nallhits_planes_W')-pl.col('trknhits_planes_W'),
    ).select(
        pl.col('run'),
        pl.col('subrun'),
        pl.col('event'),
        pl.col('^other.*$')
    )
    df = df.join(dftmp, on=selection_events(), how='inner')
    allkeys = df.columns
    allkeys = [ k.strip("_U") for k in allkeys if "_U" in k]
    for k in allkeys:
        vartomax = "trknhits_planes"
        if "trkpidpida" in k:
            vartomax='trkpidndf'
        if "trkpidchi2" in k:
            vartomax='trkpidchi2ndf'
            if f"{vartomax}_W" not in allkeys:
                vartomax = 'trkpidndf'

        if k.startswith('allcalo'):
            vartomax='nallhits_planes'
        if k.startswith('other'):
            vartomax='othernallhits_planes'
        
        df = df.with_columns(
            pl.when(
                (pl.col(f'{vartomax}_W') >= pl.col(f'{vartomax}_V')) & 
                (pl.col(f'{vartomax}_W') >= pl.col(f'{vartomax}_U'))  
                ).then(
                    pl.col(f'{k}_W')
                ).otherwise(
                    pl.when(
                        (pl.col(f'{vartomax}_V') >= pl.col(f'{vartomax}_U'))  
                    ).then(
                        pl.col(f'{k}_V')
                    ).otherwise(
                        pl.col(f'{k}_U')
                    )
                ).alias(f"{k}_B")
        )
    return df

def ismu():
    return pl.col('trkg4pdg_planes_B').abs()==13


# I tried 1e3 cut, virtually no change...
def create_mucandidates(df, cheat=False):
    if not cheat:
        dfselected = df.sort('trklen').filter(
            pl.col('trkPFPIsTrack'),
            pl.col('trklen')>0
        ).group_by(selection_events(), maintain_order=True).agg(
            pl.col('trkId').last(),
        )
    else:
        dfselected = df.sort('pnc').filter(
            ismu(),
            pl.col('trklen')>0
        ).group_by(selection_events(), maintain_order=True).agg(
            pl.col('trkId').last()
        )

    dfselected = dfselected.with_columns(
            pl.lit(True).alias("selected_mu")
    )
    df = df.join(dfselected, on=selection_events('trkId'), how='left', coalesce=True).with_columns(
        pl.col('selected_mu').fill_null(False)
    ).group_by(selection_events(), maintain_order=True).agg(
        pl.len().alias('noptions'),
        pl.all()
    ).explode(
        pl.all().exclude(selection_events('noptions'))
    )
    return df

def get_options(df:pl.DataFrame):
    df = df.drop('noptions')
    df = df.group_by(selection_events(), maintain_order=True).agg(
        pl.len().alias('noptions'),
        pl.all()
    ).explode(
        pl.all().exclude(selection_events('noptions'))
    )
    return df

def get_biggest(df:pl.DataFrame,
                tailv:dict={'trklen':1,
                       'trkmomllhd':2,
                       'trkPFPScoreIsTrack':2,
                      }
                ):
    listsort = [ k for k in  tailv.keys()]
    listvars = [ f'big_{k}' for k in  tailv.keys()]
    df = df.drop(listvars)
    for var, varsort, tt in zip(listvars, listsort, tailv.values()):
        if varsort not in df.columns:
            continue
        dfsel = df.sort(varsort,).group_by(selection_events(), maintain_order=True).agg(
            pl.col('trkId').tail(tt),
        ).explode(
            'trkId'
        ).with_columns(
            pl.lit(True).alias(var)
        )
        df = df.join(dfsel, on=selection_events('trkId'), how='left', coalesce=True).with_columns(
            pl.col(var).fill_null(False)
        )
    return df
def get_rest(df:pl.DataFrame, restpidcut = 13):

    df = get_options(df)
    df = df.filter(~((pl.col('trkpidpida_B')>=restpidcut) & (pl.col('noptions')>1)))
    df = get_options(df)
    firsttails = {
       'trklen': 1,
       'trkmomllhd': 2,
       'trkPFPScoreIsTrack': 2,
    }
    df = get_biggest(df, firsttails)


    dffirst = df.with_columns(
        pl.when(((pl.col('big_trklen')) & (pl.col('big_trkPFPScoreIsTrack')) )).then(True).otherwise(False).alias('selected_mu')
    ).filter(
        pl.col('selected_mu')
    ).sort('trklen').group_by(selection_events(), maintain_order=True).agg(
        pl.col('trkId').last()
    )

    secondtails = {
       'trklen': 1,
       'trkmomllhd': 2,
       'trkPFPScoreIsTrack': 2,
    }
    df = df.join(dffirst, on=selection_events(), how='anti')

    df = get_biggest(df, secondtails)
    df = get_options(df)

    dfsecond = df.with_columns(
        pl.when(((pl.col('big_trklen')) & pl.col('big_trkPFPScoreIsTrack')  )).then(True).otherwise(False).alias('selected_mu')
    ).filter(
        pl.col('selected_mu')
    ).sort('trklen').group_by(selection_events(), maintain_order=True).agg(
        pl.col('trkId').last()
    )

    dfselected = pl.concat([dffirst, dfsecond], how='vertical', rechunk=False)

    return dfselected

def create_mucandidates_2(df:pl.DataFrame,
                          cut=1.7e3,
                          pidacut=13,
                          pidacut_low=10,
                          restpidcut=13,
                          filter_shower=True,
                          calocut=0, # was 3
                          calocutpida=0., # was 0.1
                          calocut_new=1,
                          fineselection=False,
                          firstselect={"trklen":1,
                                      "trkmomllhd":3,
                                      "trkPFPScoreIsTrack":3},
                          afterfine = False,
                          includeall = True,
                          ):

    calovariable = "allcalo_planes_B"
    dforiginal = df
    df = df.filter(
        pl.col('trklen')>0
    )

    df = get_options(df)
    df = df.filter(
        ~((pl.col('trklen') >= cut) & (pl.col('noptions')>1)),
        # pl.col('trkPFPIsTrack'),
    )
    df = get_options(df)
    df = df.filter(
        ~((pl.col(calovariable) >= calocutpida) & (pl.col('trkpidpida_B') > pidacut) & (pl.col('noptions')>1)),
        # pl.col('trkPFPIsTrack'),
    )
    if filter_shower:
    #     # Tested: there is only ~4% of tracks tagged as muon here with PxC > 0.36
        df = get_options(df)
        df = df.filter(
        ~((pl.col(calovariable).is_between(0.3,2)) & ~(pl.col('trkPFPIsTrack')) & (pl.col('noptions')>1))
        )
        df = get_options(df)
        df = df.filter(
        ~((pl.col('trkpidpida_B') < 5) & ~(pl.col('trkPFPIsTrack')) & (pl.col('noptions')>1))# & (pl.col("allcalo")>2))
        )

    df = get_options(df)
    dfallcontained = df.filter(
        (pl.col(calovariable) >= calocut) & ~(pl.col('trkIsContained')) & (pl.col('noptions')>1)
    ).group_by(selection_events()).agg(
    ).with_columns(
        pl.lit(False).alias('allcontained')
    )
    df = df.join(dfallcontained, on=selection_events(), how='left', coalesce=True).with_columns(
        pl.col('allcontained').fill_null(True)
    )
    df = df.filter(
        ~((pl.col(calovariable) >= calocut_new) & (pl.col('trkIsContained')) & ~(pl.col('allcontained'))),
    )
    df = get_options(df)
    df = df.filter(
        ~((pl.col(calovariable) >= calocut) & ~(pl.col('trkPFPIsTrack')) & (pl.col('noptions')>1)),
    )
    df = get_options(df)

    # df = df.filter(
    #     ~((pl.col(calovariable) >= calocut) & ~(pl.col('trkPFPIsTrack')) & (pl.col('trkpidpida_B')>=pidacut_low) & (pl.col('noptions')>1)),
    # )
    # df = get_options(df)

    df = df.filter(
        ~((pl.col(calovariable) >= calocut) & (pl.col('trkPFPIsTrack')) & (pl.col('noptions')>1) & (pl.col('trkpidpida_B')>pidacut_low)),
    )
    df = get_options(df)
    # df = df.filter(
    #     ~((pl.col(calovariable) >= calocut) & (pl.col('trkPFPIsTrack')) & (pl.col('noptions')>=2) & (pl.col('trkpidpida_B')>pidacut_low)),
    # )
    df = df.with_columns(
        pl.lit(True).alias('candidate')
    )
    dfcandidates = df
    df = get_options(df)
    dfselected = df.sort(
        'trklen',
    ).group_by(selection_events(), maintain_order=True).agg(
        pl.col('trkId').last(),
    )

    if afterfine:
        dfrest = df.join(dfselected, on=selection_events(), how='anti')
        dfrest = dfrest.filter(pl.col('candidate'))
        dfrestselected = get_rest(dfrest, restpidcut=restpidcut)
        dfrestselected = dfrestselected
        dfselected = pl.concat([dfselected, dfrestselected], how='vertical', rechunk=False)

    dfrest = dforiginal.join(dfselected, on=selection_events(), how='anti')
    if includeall:
        dfrest = dfrest.filter(pl.col('trklen')>0)
        dfrest = get_options(dfrest)
        dfrest = dfrest.filter(
            ~((pl.col('trklen') >= cut) & (pl.col('noptions')>1)),
        )
        dfrestselected = dfrest.sort('trklen').group_by(selection_events()).agg(
            pl.col('trkId').last()
        )
        dfselected = pl.concat([dfselected, dfrestselected], how='vertical', rechunk=False)


    dfselected = dfselected.with_columns(
        pl.lit(True).alias("candidate"),
        pl.lit(True).alias("selected_mu"),
    )
    # dforiginal = df
    dforiginal = dforiginal.join(dfselected, on=selection_events('trkId'), how='left', coalesce=True).with_columns(
        pl.col('candidate').fill_null(False),
        pl.col('selected_mu').fill_null(False),
    )
    dfcandidates = dfcandidates.with_columns(
        change_status=pl.col('candidate')
    ).select(selection_events('trkId','change_status'))
    dforiginal = dforiginal.join(dfcandidates, on=selection_events('trkId'), how='left', coalesce=True)
    dforiginal = dforiginal.with_columns(
        candidate = pl.when(pl.col('change_status')).then(True).otherwise(pl.col('candidate'))
    ).drop('change_status')

    return dforiginal


def isparticle(pdg):
    return pl.col('trkg4pdg_planes_B').abs()==pdg
def isproton():
    return isparticle(2212)
def onlyone(n:int=1):
    return pl.col('noptions')==n
def morethan(n:int=2):
    return pl.col('noptions')>=n

def get_all_the_rest(dffull:pl.DataFrame, dfcutted:pl.DataFrame):
    return dffull.join(dfcutted.group_by(selection_events()).agg(), on=selection_events(), how='anti')
def create_proton_candidate(df:pl.DataFrame,
                            cut_pida_easy=10,
                            cut_pida_easy_shower=13,
                            maxmom=1.5,
                            maxcalo=0.8,
                            apply_shower_filter=False,
                            nopt_shower_filter=3,
                            cheat=False,
                            cheatMinPNC=0.,
                            ):
    df = df.filter(
        ~pl.col('selected_mu'),
        pl.col('trklen')>0,
    )

    if not cheat:
        df = get_options(df)
        df = df.filter(
            ((pl.col('trkpidpida_B')>cut_pida_easy) & (pl.col('trkPFPIsTrack'))) |
            ((pl.col('trkpidpida_B')>cut_pida_easy_shower) & ~(pl.col('trkPFPIsTrack')))
            # | ((onlyone()) & (pl.col('trkpidpida_B')>6) & (pl.col('trkmomrange_pr')<0.6) & (pl.col('trkPFPIsTrack')))
        ).filter(
            pl.col('trkcalo_planes_B')<maxcalo,
            pl.col('trkmomrange_pr').is_between(0,maxmom)
        )
        if apply_shower_filter:
            df = df.filter(
                ~(~(pl.col('trkPFPIsTrack')) & morethan(nopt_shower_filter))
            )
        df = get_options(df)
    else:
        df = df.filter(
            isproton(),
            pl.col('pnc') > cheatMinPNC,
        )
        df = get_options(df)
    df = df.with_columns(
        pl.lit(True).alias('selected_pr')
    )
    return df


def create_pion_candidate(df:pl.DataFrame,
                          minimum_options_after_remove_mu=2,
                          return_full_efficiency=False,
                          do_not_apply_small_cut_pida=True,
                          do_not_apply_cut_en=True,
                          do_not_apply_containement=True,
                          othercuts=True,
                          cheat=False,
                          cheatMinPNC=0.,
                          ):
    df = df.filter(
        pl.col('trklen')>0,
    )
    df = df.filter(
        ~pl.col('selected_mu'),
    )
    df = get_options(df)
    dfwithpr = df.filter(pl.col('selected_pr')).group_by(selection_events()).agg(
        pl.col('selected_pr').count().alias('npr'),
    )
    df = df.join(dfwithpr, on=selection_events(), how='left', coalesce=True)
    df = df.with_columns(
        pl.col('npr').fill_null(0)
    )
    if return_full_efficiency or cheat:
        df = df.filter(
            ~pl.col('selected_pr'),
        )
        if cheat:
            df = df.filter(
                isparticle(211),
                pl.col('pnc') > cheatMinPNC
            )
        df = get_options(df)
        df = df.with_columns(
            pl.lit(True).alias('selected_pi')
        )
        return df

    df = df.filter(
        ~pl.col('selected_pr'),
        morethan(minimum_options_after_remove_mu), # Could not find any reliable way if there is muon and no proton
        pl.col('trkPFPIsTrack')
    )
    df = df.filter(
        (do_not_apply_small_cut_pida) | ((pl.col('trkpidpida_B').is_between(2.,7)) | (pl.col('trkpidpida_B').is_between(8.5,10))),
        (do_not_apply_cut_en) | (pl.col('trkmomrange_mu').is_between(0.15, 0.5) & pl.col('trkcalo_planes_B').is_between(0.05, 0.5)),
        # pl.col('npr')>=1,
        (do_not_apply_containement) | (pl.col('trkIsContained'))
    )

    if othercuts:
        df = df.filter(
            ((pl.col('trklen')).is_between(3.5, 350))
        ).filter(
            ~((~pl.col('trkIsContained')) & (pl.col('trklen')<10))
        )

    df = df.with_columns(
        pl.lit(True).alias('selected_pi')
    )

    return df


def K_mu_cont():
    return (pl.col('trkmomrange_mu')**2+0.1057**2).sqrt()-0.1057
def K_mu_notcont():
    return (pl.col('trkmomllhd')**2+0.1057**2).sqrt()-0.1057

def K_from_p(p:str, mass:float):
    return (pl.col(p)**2+mass**2).sqrt()-mass
def p_from_K(K:str, mass:float):
    return (pl.col(K)**2 + 2*pl.col(K)*mass).sqrt()

def K_from_p_v(p:float, mass:float):
    return (p**2+mass**2).sqrt()-mass
def p_from_K_v(K:float, mass:float):
    return (K**2 + 2*K*mass).sqrt()

def trkcalo(W:str):
    return pl.col(f'trkcalo_planes_{W}')

def simple_energy(df:pl.DataFrame,
                  forceaddmumass=True,
                  ):
    df = df.with_columns(
        lepmom = (pl.col('lepen')**2-0.1057**2).sqrt(),
        Pmu = pl.when((pl.col('trkIsContained')) & (pl.col('selected_mu'))).then('trkmomrange_mu').otherwise(
            pl.when(~(pl.col('trkIsContained')) & (pl.col('selected_mu') & (pl.col('trkmomllhd')>0))).then('trkmomllhd').otherwise(0)
        ),
        Kmu = pl.when((pl.col('trkIsContained')) & (pl.col('selected_mu'))).then(K_mu_cont()).otherwise(
            pl.when(~(pl.col('trkIsContained')) & (pl.col('selected_mu')) & (pl.col('trkmomllhd')>0)).then(K_mu_notcont()).otherwise(0)),
        # Kmu = pl.when((pl.col('trkIsContained')) & (pl.col('selected_mu'))).then(K_mu_cont()).otherwise(
        #     pl.when(~(pl.col('trkIsContained')) & (pl.col('selected_mu'))).then(pl.when(pl.col('trkmomllhd')>0).then(K_mu_notcont()).otherwise(pl.col('trkmomrange_mu'))).otherwise(0)),
    ).with_columns(
        Emu = pl.when(pl.col('Pmu')>0).then((pl.col('Pmu')**2 + 0.1057**2).sqrt()).otherwise(0),
        Kmucalo = pl.when((pl.col('selected_mu')) & (pl.col('Kmu')>0)).then('trkcalo_planes_W').otherwise(0),
    ).with_columns(
        Ehad = pl.col('allcalo_planes_W') - pl.col('Kmucalo'),
        Elep = pl.col('Kmu')+pl.when((pl.col('Kmu')>0) | (forceaddmumass)).then(0.1057).otherwise(0)
    ).with_columns(
        Etotal = pl.col("Elep") + pl.col("Ehad")
    )
    df = df.sort(selection_events('selected_mu')).group_by(selection_events()).agg(
        pl.all().last(),
        (pl.col('selected_mu').last()).alias('hasmu'),
        (pl.col('enu_truth').last()-pl.col('lepen').last()).alias('had_truth'),
    )
    return df


def complex_energy2(df:pl.DataFrame,
                   W="W",
                   force_bigger_pr=False,
                   do_not_force_bigger_pi=True,
                   return_full=False,
                   try_hard=False,
                   ):

    if try_hard:
        df = df.with_columns(
            pl.when(trkcalo(W)==0).then(trkcalo("B")).otherwise(trkcalo(W)).alias(f'trkcalo_planes_{W}')
        )

    df=df.with_columns(
        lepmom = (pl.col('lepen')**2-0.1057**2).sqrt(),
        trkg4mom = (pl.col(f'trkg4en_planes_{W}')**2 - pl.col(f'trkg4mass_planes_{W}')**2).sqrt(),
        trkg4K = (pl.col(f'trkg4en_planes_{W}') - pl.col(f'trkg4mass_planes_{W}')),
        Pmu = pl.when((pl.col('trkIsContained')) & (pl.col('selected_mu'))).then('trkmomrange_mu').otherwise(
            pl.when(~(pl.col('trkIsContained')) & (pl.col('selected_mu') & (pl.col('trkmomllhd')>0))).then('trkmomllhd').otherwise(0)
        ),
        Kmu = pl.when(
            pl.col('selected_mu')
            ).then(
                pl.when(
                    (pl.col('trkIsContained')) & (pl.col('trkmomrange_mu')>0)
                ).then(
                    K_mu_cont()
                ).otherwise(
                    pl.when(~(pl.col('trkIsContained'))).then(
                        pl.when( pl.col('trkmomllhd')>0).then( K_mu_notcont()).otherwise(trkcalo(W))
                    ).otherwise(trkcalo(W)))
            ).otherwise(0),
        Kpr = pl.when((pl.col('selected_pr'))).then(
            pl.when(pl.col('trkmomrange_pr') > 0).then(K_from_p('trkmomrange_pr', 0.9383)).otherwise(trkcalo(W))
        ).otherwise(0),
        Kpi = pl.when((pl.col('selected_pi'))).then(
            pl.when('trkIsContained').then(
                pl.when(
                    (pl.col('trkmomrange_mu') > 0) & (do_not_force_bigger_pi | (K_mu_cont() > trkcalo(W)))
                ).then(
                        K_mu_cont()
                ).otherwise(trkcalo(W))
            ).otherwise(
                pl.when(
                    (pl.col('trkmomllhd') > 0) & (do_not_force_bigger_pi | (K_mu_notcont() > trkcalo(W)))
                ).then(
                    K_mu_notcont()
                ).otherwise(trkcalo(W))
            )
        ).otherwise(0)
    )

    if force_bigger_pr:
        df = df.with_columns(
            Kpr = pl.when((pl.col('selected_pr'))).then(
                pl.when((pl.col('trkmomrange_pr') > 0) & (K_from_p('trkmomrange_pr', 0.9383) > trkcalo(W))).then(
                    K_from_p('trkmomrange_pr', 0.9383)
                ).otherwise(trkcalo(W))
            ).otherwise(0)
        )
    df = df.with_columns(
        Emu = pl.when(pl.col('Pmu')>0).then((pl.col('Pmu')**2 + 0.1057**2).sqrt()).otherwise(0),
        Kmu = pl.when((pl.col('selected_mu'))).then(
            pl.when(pl.col('Kmu')>0).then(pl.col('Kmu')).otherwise(trkcalo(W))
        ),
        Kpr = pl.when((pl.col('selected_pr'))).then(
            pl.when(pl.col('Kpr')>0).then(pl.col('Kpr')).otherwise(trkcalo(W))
        ),
        Kpi = pl.when((pl.col('selected_pi'))).then(
            pl.when(pl.col('Kpi')>0).then(pl.col('Kpi')).otherwise(trkcalo(W))
        ),
        OtherPFPs = pl.when( ~(pl.col('selected_mu')) & ~(pl.col('selected_pr')) & ~(pl.col('selected_pi'))).then( trkcalo(W) ).otherwise(0)
    )
    if return_full:
        return df

    df = df.sort('selected_mu').group_by(selection_events(),maintain_order=True).agg(
        pl.col('selected_mu').filter(pl.col('selected_mu')).any(),
        pl.col('enu_truth').last(),
        pl.col('lepen').last(),
        (pl.col('enu_truth').last()-pl.col('lepen').last()).alias('had_truth'),
        pl.col(f'otherallcalo_planes_{W}').last(),
        pl.col('Kmu').sum(),
        pl.col('Kpi').sum(),
        pl.col('Kpr').sum(),
        pl.col('OtherPFPs').sum(),
        pl.col('trkIsContained').filter(pl.col('selected_mu')).last(),
        (pl.col('npi').last()*0.1393).alias('mass_to_add'),
    ).with_columns(
        Elep = pl.col('Kmu')+0.1057,
        Ehad = pl.col('Kpr') + pl.col('Kpi') + pl.col('OtherPFPs') + pl.col('mass_to_add') + pl.col(f'otherallcalo_planes_{W}'),
    ).with_columns(
        Etotal = pl.col("Elep") + pl.col("Ehad")
        # Etotal = pl.col('Kmu') + pl.col('Kpr') + pl.col('Kpi') + pl.col('OtherPFPs') + pl.col('mass_to_add') + pl.col(f'otherallcalo_planes_{W}'),
    )
    return df


def get_n_prpi(df:pl.DataFrame):

    dfwith = df.filter(pl.col('selected_pr')).group_by(selection_events()).agg(
        pl.col('selected_pr').count().alias('npr'),
    )
    df = df.join(dfwith, on=selection_events(), how='left', coalesce=True)
    df = df.with_columns(
        pl.col('npr').fill_null(0)
    )
    if 'selected_pi' not in df.columns:
        return df
    dfwith = df.filter(pl.col('selected_pi')).group_by(selection_events()).agg(
        pl.col('selected_pi').count().alias('npi'),
    )
    df = df.join(dfwith, on=selection_events(), how='left', coalesce=True)
    df = df.with_columns(
        pl.col('npi').fill_null(0)
    )
    return df

def join_pi_candidate(df_pi:pl.DataFrame, df:pl.DataFrame, t='pi'):
    df = df.join(df_pi.select(selection_events('trkId', f'selected_{t}')), on=selection_events('trkId'), how='left', coalesce=True)
    df = df.with_columns(
        pl.col(f'selected_{t}').fill_null(False)
    )
    df = get_n_prpi(df)
    return df

