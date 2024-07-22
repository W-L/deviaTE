import logging
from pathlib import Path

import pandas as pd
import numpy as np
from numpy.typing import NDArray
import plotnine as p9
import gffutils




def _ref_allele(snp_index: int, refbases: NDArray, length: int) -> int:
    """
    Return the index of the reference allele in the long format dataframe.
    :param snp_index: Index of variant sites
    :param refbases: array of reference alleles
    :param length: length of the entire ref sequence
    :return: index of reference allele in long dataframe
    """
    if refbases[snp_index] == 'A':
        return snp_index
    elif refbases[snp_index] == 'C':
        return snp_index + length
    elif refbases[snp_index] == 'G':
        return snp_index + 2 * length
    else:
        return snp_index + 3 * length


def _snp_frame(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a long format data frame that houses the variants
    results in a DataFrame with 4 * length and 5 columns
    :param df: Input dataframe
    :return: Frame with just variants for visualisation
    """
    # concatenate the base coverage and create factor columns
    length = len(df)
    counts_col = np.concatenate([df['A'], df['C'], df['G'], df['T']])
    pos_col = np.tile(np.arange(0, length), 4)
    fam_col = np.repeat(df['#TEfam'], 4)
    sample_col = np.repeat(df['sample_id'], 4)
    base_col = np.repeat(['A', 'C', 'G', 'T'], length)
    # get snp positions
    snps = df[df['snp'] == True]['pos']
    polys = snps.values.astype(int)
    # remove reference allele from SNPs to visualise just the variant allele
    polys_ref_indices = [_ref_allele(idx, df['refbase'], length) for idx in polys]
    polys_ref_indices_sorted = sorted(polys_ref_indices)
    base_col[polys_ref_indices_sorted] = 'X'
    # create the new long format df
    lf = pd.DataFrame({'fam_col': fam_col, 'sample_col': sample_col,
                       'base_col': base_col, 'pos_col': pos_col, 'counts_col': counts_col})
    snps = lf[lf['pos_col'].isin(polys)]
    snps_noref = snps[snps['base_col'] != 'X']
    snps_noref['base_col'] = pd.Categorical(snps_noref['base_col'], categories=['A', 'C', 'G', 'T'])
    return snps_noref


def _load_annotation_db(gff_filename: str) -> gffutils.FeatureDB | None:
    """
    Load a database of annotations into memory from a gff file
    :param gff_filename: Filename of the gff annotations file
    :return: Database object of the gffutils package
    """
    if not Path(gff_filename).is_file():
        logging.info("gff file does not exist. Skipping loading of annotations")
        return
    db = gffutils.create_db(data=gff_filename, dbfn=":memory:")
    return db


def _get_annotations(frame: pd.DataFrame, gff_db: gffutils.FeatureDB) -> list :
    """
    Load annotations for a specific sequence
    :param frame: dataframe used just to get the sequence name
    :param gff_db: Feature databse from the gffutils package
    :return: list of annotated features for selected sequence
    """
    sid = frame['#TEfam'].iloc[0]
    length = len(frame)
    # first try if an annotation for this sequence is available
    try:
        s = gff_db[sid]
    except gffutils.exceptions.FeatureNotFoundError:
        logging.info(f"No annotation found for {sid}")
        return []
    # fill a list with all available annotations
    annotations = []
    for f in gff_db.region(seqid=sid, start=s.start, end=s.end):
        if f.start != 1 and f.end != length:
            annotations.append((f.featuretype, f.start - 1, f.end))
    return annotations


def _prep_annotation_df(frame: pd.DataFrame, annotations: str) -> pd.DataFrame:
    """
    Prepare a dataframe to visualise the annotations
    :param frame: Results dataframe
    :param annotations: Path to the annotation file
    :return: dataframe with coordinates to plot the annotations
    """
    # create a gff database to grab annotations from
    gff_db = _load_annotation_db(annotations)
    if not gff_db:
        return pd.DataFrame()
    ann = _get_annotations(frame, gff_db)
    # create a frame to visualise annotations
    ann_df = {'xmin': [], 'xmax': [], 'ymin': [], 'ymax': [], 'colour': []}
    ymax = -(max(frame['cov']) * 0.05)
    ymax_fac = -(max(frame['cov']) * 0.025)
    for tp, start, end in ann:
        ann_df['xmin'].append(start)
        ann_df['xmax'].append(end)
        ann_df['ymax'].append(ymax)
        ymax = ymax + (ymax_fac)
        ann_df['ymin'].append(ymax)
        ann_df['colour'].append(tp)
    adf = pd.DataFrame(ann_df)
    return adf


def visualise(input_file: str, annotations: str = None) -> None:
    """
    Visualise the results of one analysis. Optionally adds annotations
    :param input_file: Input results file from deviate analysis
    :param annotations: Optional file with annotations
    :return:
    """
    # logging.info('Loading data:', input_file)
    frame = pd.read_csv(input_file, sep=" ", skiprows=1)
    # create dataframe to visualise variants
    snp_f = _snp_frame(frame)
    # polygon total coverage
    pg = pd.DataFrame({'px': np.concatenate((frame['pos'], np.flip(frame['pos']))),
                       'py': np.concatenate((frame['cov'], np.zeros(len(frame['cov']))))})
    # polygon hq coverage
    pg_hq = pd.DataFrame({'px': np.concatenate((frame['pos'], np.flip(frame['pos']))),
                          'py': np.concatenate((frame['hq_cov'], np.zeros(len(frame['hq_cov']))))})
    # colours for the bases
    base_colors = ['#1b9e77', '#d95f02', '#7570b3', '#e6ab02']
    # put together the visualisation
    baseplot = (
        p9.ggplot() +
        p9.geom_polygon(data=pg,
                        mapping=p9.aes(x='px', y='py'),
                        fill='lightgrey',
                        color='lightgrey') +
        p9.geom_polygon(data=pg_hq,
                        mapping=p9.aes(x='px', y='py'),
                        fill='grey',
                        color='grey') +
        p9.geom_bar(data=snp_f,
                    mapping=p9.aes(x='pos_col', y='counts_col', fill='base_col'),
                    stat='identity',
                    width=1) +
        p9.labs(fill='') +
        p9.ylab("coverage") +
        p9.xlab("") +
        p9.ggtitle(f"{frame['#TEfam'].iloc[0]} {frame['sample_id'].iloc[0]}") +
        p9.scale_fill_manual(values=base_colors) +
        p9.theme_minimal() +
        p9.theme(legend_position='bottom')
    )
    # add rectangles for the annotations if any are given
    if annotations:
        adf = _prep_annotation_df(frame, annotations)
        if len(adf) > 0:
            baseplot += p9.geom_rect(data=adf,
                                     mapping=p9.aes(
                                         xmin='xmin',
                                         xmax='xmax',
                                         ymin='ymin',
                                         ymax='ymax',
                                         color='colour'
                                     ),
                                     fill='none',
                                     size=1)
            baseplot += p9.scale_color_gray()
            baseplot += p9.labs(colour="")
    # save visualisation to file
    baseplot.save(f"{Path(input_file).name}.pdf", width=10, height=5)

