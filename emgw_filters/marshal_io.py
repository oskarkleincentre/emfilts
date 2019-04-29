"""
Basic io module using marshaltools to extract data from the marshal


Notes (This is also the file used for ZTF rates, and am simply moving here)
R. Biswas
- Changes to get_marshal_photometry to avoid trying to get redshifts
"""
__all__ = ['extract_data', 'extract_latest_marshal_table', 'clean_ztf_marshal',
           'class_aliases', 'get_summary', 'get_marshal_photometry']

import numpy as np
import pandas as pd
from tqdm import tqdm
import healpy as hp
from marshaltools import ProgramList

# define aliases
snia_norm = ['SN Ia', 'SN Ia ', 'SNIa', 'SN Ia-norm']
sn_cc = ['SN Ib', 'SN Ic', 'Ic-BL', 'SN Ic-BL', 'SN IIL', 'SN Ibn',
         'SN II-pec', 'SN Ic', 'SN IIb', 'SN II',
         'SN Ib/c', 'SN II ', 'SN IIn', 'SN IIP']
statics = ['galaxy', 'rock', 'Galaxy', 'star','star' ]
agns = ['AGN',  'NLS1', 'QSO', 'AGN', 'LINER', "CLAGN", "Q"]
other_vars = ['LBV', 'varstar','CV', "CV Candidate", 'Var',  'blazar']
other_trans = ['Gap I - Ca-rich', 'Gap', 'ILRT', 'SN']
tde = ['TDE']
undecideds = ['SN?', 'SN? ', 'SN Ia?', 'SN Ic?', "AGN?", "star?",
              'TDE? ', 'SN?', 'Star?','CV?', 'blazar?', 'AGN? ']
slsn = ['SLSN-I', 'SLSN-II', 'SLSN-R']
bads = ['unk', 'duplicate', 'None', None, 'old', 'unclassified', 'bogus']
others = set(sn_cc).union(statics).union(agns).union(other_vars).union(other_trans).union(tde).union(undecideds).union(slsn).union(bads)

class_aliases = (('snia_norm', snia_norm), 
                 ('sn_cc', sn_cc),
                 ('statics', statics),
                 ('agns', agns),
                 ('other_vars', other_vars),
                 ('other_trans', other_trans),
                 ('tde', tde),
                 ('undecideds', undecideds),
                 ('slsn', slsn),
                 ('bads', bads),)

cc_aliases = (('snibc',['SN Ib', 'SN Ic', 'Ic-BL', 'SN Ic-BL', 'SN Ibn', 'SN Ib/c']), 
              ('snii',['SN IIL', 'SN II-pec', 'SN II', 'SN II ', 'SN IIn', 'SN IIP']),
              ('sniib',['SN II b'])
             )

def snia_pec(seq):
    return list((set(seq) - set(snia_norm)) - others)

def extract_latest_marshal_table(fname,
                                 program_name='Redshift Completeness Factor'):
    """
    extract the latest summary (object name, ra, dec) from the marshal using
    `marshaltools` for a particular program

    Parameters
    ----------
    fname : string 
        absolute path of csv file to be output to
    program_name : string, defaults to `Redshift Completeness Factor`
        program name under which host is saved in

    Returns
    -------
    None 
    """
    pl = ProgramList(program_name)
    
    df  = pl.table.to_pandas()
    df.to_csv(fname, index=False)
    return df

def extract_data(name, sources):
    """
    extract the data to download from the marhsall sources

    Parameters
    ---------
    name : string
	name of the ZTF object in marshaltools

    sources : dict
        dictionary of sources with Marshal Information. Expected to
	be obtained by the `sources` attribute of the instance of
       	`marshaltools.ProgramList`

    Returns
    ------
    l : tuple of quantities for downloads, `None` if unsuccessful 
    headers : the names of the downloaded quantities, `None` if unsuccessful
    missing_name : The name of the sources if there is a problem leading
	to the downloads. `None` if succesfully downloaded.
    """
    try:
        data = sources[name]
        band, mag = data['mag']
        l = (name, data['classification'],
             data['ra'], data['dec'], data['redshift'],
             data['creationdate'], data['lastmodified'],
             data['iauname'], band, mag,)
        headers = ('name', 'classification', 'ra', 'dec', 'redshift',
                   'creationdate', 'lastmodified', 'iauname', 'band',
                   'mag')
    except KeyError:
        return (None, None, name)
    return l, headers, None


def get_summary(names, sources, healpix_NSIDE=None, cleanup=True, aliases=None,
                divide_cc=True):
    """ return the summary

    Parameters
    ----------
    names : 
    sources :
    healpix_NSIDE :
    cleanup : Bool
        if True, cleanup the classifications to a new class
    aliases:
        list of pairs (names, lists) of aliases.

    Returns
    -------
    summary : `pd.DataFrame` with summary of each object found
    missing : Object Ids that were not found.
    """
    tdas = list()
    missing = list()

    # iterate through the objects
    for name in tqdm(names):
        lst, headers, missing_name = extract_data(name, sources)  
        if lst is None:
            assert missing_name == name
            missing.append(name)
        else:
            tdas.append(lst)

    # Turn into a dataframe
    df = pd.DataFrame(tdas, columns=headers)
    
    # Add the healpix Column
    if healpix_NSIDE is not None:
        df['hp_{}'.format(healpix_NSIDE)] = hp.ang2pix(healpix_NSIDE,
						       df.ra, df.dec,
                                                       lonlat=True, nest=False)

    # Change names of columns for classes
    df.rename(columns=dict(name='tid', classification='mclass'), inplace=True)

    if aliases is None:
        aliases = list(class_aliases)

        aliases.append(('snia_pec', snia_pec(df.mclass.unique())))

    if cleanup:
        summary = clean_ztf_marshal(df, aliases)

    if divide_cc:
        summary_t = clean_ztf_marshal(df.query('tclass == "sn_cc"'),
                                      cc_aliases, 'cc_class')
        summary = df.join(summary_t.cc_class)

    return summary.set_index('tid'), missing
    
def clean_ztf_marshal(summary, aliases, classname='tclass'):
    """ clean up different classifications for marshalls
    
    Parameters
    ----------
    summary :
    aliases : 

    Returns
    -------
    summary : `pd.DataFrame`
    """
    vals = []

    # `None` gives a lot of trouble. Replace by `Unclassified`
    # This is the only change in the marshal classification we will do
    for v in summary['mclass'].values:
        if v is None:
            vals.append('Unclassified')
        else:
            vals.append(v)

    # Set up new classification for TDAs
    summary[classname] = vals
    
    for (k, v) in aliases:
        summary[classname] = summary[classname].replace(v,k)
    
    return summary 


def get_marshal_photometry(program, names, verbose=False):
    """download the marshal photometry for a particular program and a
    sequence of object names.


    Parameters
    ----------
    program: string, mandatory
	name of the program.
    names: Sequence of strings
	sequence of ZTF object names
    verbose: Bool
	Whether to provide verbose output

    Returns
    -------
    metadata
    photometry
    bad_tdas: 
    """
    bad_tdas = []
    
    summary = []
    phot = []
    
    for name in tqdm(names):
        try:
            lc = program.get_lightcurve(name)
            # summary.append([name, lc.ra, lc.dec, np.float(lc.redshift), lc.mwebv])
            pt = lc.table.to_pandas()
            if verbose:
               print(lc.table)
               print(lc.table.to_pandas())
            # if 'absmag' not in pt.columns:
            #    pt['absmag'] = None
            pt['name'] = name
            phot.append(pt)
            
            if verbose:
                print(len(pt))
        except:
            bad_tdas.append(name)

    # Now create the right dataframes 
    photometry = pd.concat(phot)
    photometry.rename(columns=dict(magpsf='mag', name='tid',
                                   sigmamagpsf='mag_err', limmag='m5',
                                   jdobs='jd'), inplace=True)
    # metadata = pd.DataFrame(summary, columns=('tid', 'ra', 'dec', 'redshift',
    #                                         'mwebv'))
    # metadata.set_index('tid', inplace=True)

    # It is possible to have bad objects, but these should match
    # assert photometry.tid.unique().size == len(metadata)
    assert len(names) == photometry.tid.unique().size + len(bad_tdas)

    return photometry, bad_tdas
