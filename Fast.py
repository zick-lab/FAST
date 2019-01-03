'''
Fast is a tool for analyzing the structure and the characteristics of segments, based on CNVs and CVs
'''

__version__ = "$Revision: 4 $"
# $Source$

import ConfigParser
import argparse
import csv
import datetime
import os
import re
import shutil
import time
import pandas as pd
import numpy as np
from collections import OrderedDict


from focusedGenesProcessing import FL, segments2FocusedGenes
from excelUtils import write_excel


mut_types = ['ID', 'TR', 'DM', 'Other']

CP = ConfigParser.ConfigParser()

def readCNVData(cnvFile, remove_unmappble, unify_seg = True):
    '''
    read cnv file, clear and filter it
    :param cnvFile: cvn file
    :param remove_unmappble: flag to indicate whether to remove unmappble areas
    :return: cleared, filtered cnv dataframe
    '''

    print 'readCNVData'
    dfcnv = pd.read_csv(cnvFile, prefix='CNV', header=None, names=['chr','start','end', 'CopyNumber', 'VariationType'],sep=None)
    dfcnv['start'] = dfcnv['start'].astype('int32')
    dfcnv['end'] = dfcnv['end'].astype('int32')
    dfcnv['chr'] = dfcnv['chr'].astype('str')
    dfcnv = dfcnv[dfcnv['chr'].apply(lambda x: re.match('[^0-9xy]',x) is None)]

    # filtering

    # filter only the gain
    dfcnv = dfcnv[dfcnv['VariationType']=='gain']
    # filter gain above some threshold
    dfcnv = dfcnv[dfcnv['CopyNumber']>=CP.getint('PARAMS','CNV_MINIMAL_CN')]

    # unify consecutive segments
    if unify_seg:
        dfcnv['prev_chr'] = dfcnv.chr.shift()
        dfcnv['prev_end'] = dfcnv.end.shift()
        dfcnv['not_same_segment'] = ~((dfcnv.prev_chr==dfcnv.chr) & (dfcnv.start == dfcnv.prev_end))
        dfcnv.ix[0,'not_same_segment'] = True
        dfcnv['Segment ID'] = np.cumsum(dfcnv.not_same_segment)
        dfcnv['Segment ID'] = dfcnv['Segment ID'].astype('int32')
        dfcnv = dfcnv[:len(dfcnv)-1]
        grouped = dfcnv.groupby('Segment ID', as_index=False)
        dfcnv = grouped.agg(OrderedDict((('chr','first'),
                                           ('start','first'),
                                           ('end','last'),
                                           ('CopyNumber', lambda x: int(np.mean(x))),
                                           ('VariationType','first'))))
    dfcnv['start'] = dfcnv['start'].astype('int32')
    dfcnv['end'] = dfcnv['end'].astype('int32')

    dfcnv['segment_size']=dfcnv['end']-dfcnv['start']
    dfcnv['Segment ID'] = np.array(range(len(dfcnv)))

    # remove unmappble areas
    if remove_unmappble:
        dfcnv = removeUnmap(dfcnv)

    dfcnv.sort_values(by='Segment ID', ascending=True, inplace=True, na_position='last')

    return dfcnv


def removeUnmap(dfcnv):
    '''
    util function, to remove unmappble areas from cnv dataframe
    :param dfcnv:
    :return: dfcnv, after removing unmappble areas
    '''
    THRESHOLD = CP.getint('PARAMS','WINDOW')*CP.getint('PARAMS','WINDOW_MULT')

    dfunmap = pd.read_csv(CP.get('FILES','UNMAPPBLE_INP'))
    dfunmap['chr'] = dfunmap['chr'].apply(lambda x: x.replace('chr',''))
    dftemp = pd.merge(dfcnv,dfunmap,how='left',on='chr', sort = False)
    dftemp['isUnmap'] = (dftemp['start'] <= (dftemp['chromEnd']+THRESHOLD)) & (dftemp['end'] >= (dftemp['chromStart']-THRESHOLD))
    dftemp=dftemp[dftemp['isUnmap']==True]
    dfcnv = pd.merge(dfcnv,dftemp[['chr','start','end','isUnmap']],how='left',on=['chr','start','end'],sort=True)
    dfcnv = dfcnv[pd.isnull(dfcnv['isUnmap'])]
    del dfcnv['isUnmap']
    return dfcnv


def readSVData(svFile):
    '''
    Read the SV file, perform some filtering and return a dataframe
    :param svFile: sv file
    :return: cleared, filtered sv dataframe
    '''
    print 'readSVData'

    col_names = ['Chr1', 'Pos1', 'Orientation1', 'Chr2', 'Pos2', 'Orientation2', 'Type', 'Size', 'Score', 'num_Reads']
    dfsv = pd.read_csv(svFile, sep='\t', comment='#', header=None, usecols=range(len(col_names)), names=col_names)

    print 'dfsv read'

    # filtering records records
    dfsv = dfsv[dfsv['Score']>CP.getint('PARAMS','SV_MINIMAL_SCORE')]
    dfsv = dfsv[dfsv['num_Reads']>CP.getint('PARAMS','SV_MINIMAL_NUM_READS')]

    dfsv['Chr1'] = dfsv['Chr1'].astype('str')
    dfsv['Chr2'] = dfsv['Chr2'].astype('str')
    dfsv['Chr1'] = dfsv['Chr1'].apply(lambda x: x.replace('chr',''))
    dfsv['Chr2'] = dfsv['Chr2'].apply(lambda x: x.replace('chr',''))
    dfsv = dfsv[(dfsv['Chr1'].apply(lambda x: re.match('[^0-9xy]',x) is None)) & (dfsv['Chr2'].apply(lambda x: re.match('[^0-9xy]',x) is None))]

    # remove small deletions
    dfsv = dfsv[~(((dfsv['Type']=='DEL') | (dfsv['Type']=='ITX')) & (dfsv['Chr1']==dfsv['Chr2']) & (abs(dfsv['Pos2']-dfsv['Pos1'])<CP.getint('PARAMS','SV_MINIMAL_DELETIONS_REMOVE')))]

    # digitize
    win = CP.getint('PARAMS','WINDOW')
    dfsv['Pos1'] = dfsv['Pos1'].apply(lambda x: (x//win)*win)
    dfsv['Pos2'] = dfsv['Pos2'].apply(lambda x: (x//win)*win)

    # and one additional filtering
    dfsv=unifySvRecords(dfsv)
    return dfsv


def unifySvRecords (dfsv):
    '''
    unify close enough svs
    :param dfsv
    :return: dfsv with unified closed svs
    '''
    dfsv.sort_values(by=['Chr1','Pos1','Chr2','Pos2'],ascending=[True]*4,inplace=True)
    dfsv['unify'] = np.nan
    counter = 0
    row_it = dfsv.iterrows()
    i, last = row_it.next()
    dfsv.loc[i,'unify'] = counter

    for i, row in row_it:
        if not (row['Chr1']==last['Chr1'] and row['Chr2']==last['Chr2'] and (abs(row['Pos1']-last['Pos1']))<1000 and (abs(row['Pos2']-last['Pos2']))<1000 and row['Type']==last['Type']):
            counter += 1
        dfsv.loc[i,'unify'] = counter
        last = row

    grouped = dfsv.groupby('unify', as_index=False)
    dfsv = grouped.agg(OrderedDict((('Chr1','first'),
                        ('Pos1', lambda x: int(np.mean(x))),
                        ('Orientation1', 'first'),
                        ('Chr2', 'first'),
                        ('Pos2', lambda x: int(np.mean(x))),
                        ('Orientation2', 'first'),
                        ('Type', 'first'),
                        ('Size', 'sum'),
                        ('Score', lambda x: int(np.mean(x))),
                        ('num_Reads', 'sum'),
                                    )))
    del dfsv['unify']
    return dfsv


def assignSegmentsForSVs (dfcnv, dfsv):
    '''
    Assign segments to each <chr,location> in breakpoint
    :param dfcnv: df of cnv
    :param dfsv: df of sv
    :return: dfsv with segment for each <chr,location> in breakpoint
    '''
    THRESHOLD = CP.getint('PARAMS','WINDOW')*CP.getint('PARAMS','WINDOW_MULT')

    for j in [1, 2]:
        dfsv['Segment%d ID' % j] = None
        dfsv['edge%d' % j] = None

    grouped_segments = dfcnv.groupby('chr')

    for i, row in dfsv.iterrows():
        for j in [1,2]:
            if row['Chr%d' % j] in grouped_segments.groups:
                segments_in_chr = grouped_segments.get_group(row['Chr%d' % j])
                for seg_row, seg_record in segments_in_chr.iterrows():
                    if (seg_record['start']-THRESHOLD) <= row['Pos%d' % j] <= (seg_record['end'] + THRESHOLD):
                        dfsv.loc[i, 'Segment%d ID' % j] = seg_record['Segment ID']
                        if seg_record['start'] + THRESHOLD > row['Pos%d' % j]:
                            dfsv.loc[i, 'edge%d' % j] = 'left'
                        elif seg_record['end'] - THRESHOLD < row['Pos%d' % j]:
                            dfsv.loc[i, 'edge%d' % j] = 'right'
                        else:
                            dfsv.loc[i, 'edge%d' % j] = 'middle'
                        break

    return dfsv


def createAmpliconFromConnectedSegments(max_segment, dfsv, dfcnv):
    '''
    identify connected segments (= amplicons) based on breakpoints, and assign them an Amplicon ID
    :param max_segment: maximal number of segments
    :param dfsv: df of sv
    :param dfsegments: df of segments
    :return: dfsegments, dfamplicons
    '''

    def dfs(a, visited, u, val):
        """ DFS for finding connected segments"""
        for v, temp in enumerate(a[u]):
            if a[u][v] == 0:
                continue
            if visited.has_key(v) == False:
                visited[v] = val
                dfs(a, visited, v, val)

    relevant_dfsv = dfsv.dropna(subset=['Segment1 ID','Segment2 ID'])

    a = np.zeros((max_segment+1, max_segment+1))
    for i, row in relevant_dfsv.iterrows():
        a[row['Segment1 ID']][row['Segment2 ID']] = 1
        a[row['Segment2 ID']][row['Segment1 ID']] = 1

    con_comp = 0
    visited = {}
    for u in range(max_segment+1):
        if not visited.has_key(u):
            visited[u] = con_comp
            dfs(a,visited,u,con_comp)
            con_comp += 1

    # dfsegments['Amplicon ID']=dfsegments['Segment ID'].map(visited.get)
    dfcnv['Amplicon ID'] = dfcnv['Segment ID'].map(visited.get)

    dfamplicons = pd.DataFrame(dfcnv[['Amplicon ID','Segment ID','segment_size']].groupby('Amplicon ID').agg({'Segment ID':'count', 'segment_size':'sum',}))
    dfamplicons = dfamplicons[dfamplicons['segment_size']>0]
    dfamplicons.reset_index(inplace=True)
    dfamplicons.rename(columns={'Segment ID':'nsegments_in_amplicon','segment_size': 'amplicon_size'}, inplace=True)

    return dfcnv, dfamplicons


def updateAmpliconsInSVs(dfcnv, dfsv):
    '''
    Assign Amplicon ID for segments in dfcnv and dfsv
    :param dfcnv: df of cnv
    :param dfsv: df of sv
    :param dfsegments: df of segments
    :return: updated dfcnv and dfsv
    '''
    """ assign amplicons to segments and to breakpoints"""
    dfsv = pd.merge(dfsv,dfcnv[['Segment ID','Amplicon ID']],how='left',left_on = 'Segment1 ID',right_on = 'Segment ID',sort=False)
    del dfsv['Segment ID']

    dfsv = pd.merge(dfsv,dfcnv[['Segment ID','Amplicon ID']],how='left',left_on = 'Segment2 ID',right_on = 'Segment ID',sort=False, suffixes=['1','2'])
    del dfsv['Segment ID']

    dfsv = dfsv.rename(columns={'Amplicon ID1': 'Amplicon1 ID', 'Amplicon ID2': 'Amplicon2 ID'})
    dfsv.loc[dfsv['Amplicon1 ID'].isnull(),'Amplicon1 ID'] = dfsv['Amplicon2 ID']
    del dfsv['Amplicon2 ID']
    dfsv = dfsv.rename(columns={'Amplicon1 ID':'Amplicon ID'})
    dfsv.dropna(inplace=True, subset=['Amplicon ID'])
    dfsv['Amplicon ID'] = dfsv['Amplicon ID'].astype(int)

    # dfcnv = pd.merge(dfcnv,dfsegments[['Segment ID','Amplicon ID']],how='left',on='Segment ID',sort=False)

    return dfcnv, dfsv


def setStructureTypeToSV(dfsv):
    if dfsv.empty:
        return dfsv

    dfsv['Structure Type'] = None

    for i, row in dfsv.iterrows():
        structure_type = 'Other'
        if not (pd.isnull(row['Segment1 ID']) or pd.isnull(row['Segment2 ID'])):
            if row['Segment1 ID']==row['Segment2 ID']:
                if row['edge1']=='middle' and row['edge2']=='middle' and row['Type']== 'DEL':
                    structure_type = np.NaN
#michalin
#                elif row['edge1']==row['edge2'] and not row['edge1']=='middle':
                elif row['edge1']==row['edge2']:
                    if row['Type']=='INV':
                        structure_type = 'ID'
                elif row['Type']=='ITX' and not (row['edge1']=='middle' and row['edge2']=='middle'):
                    structure_type = 'TR'
            else:
                structure_type = 'DM'
        else:
            structure_type = np.NaN
        dfsv.loc[i, 'Structure Type'] = structure_type

    return dfsv

def assignAmpliconStructureType(dfamplicon, dfsv):
    if dfsv.empty:
        return pd.DataFrame()

    grouped_by_amp = dfsv.groupby(['Amplicon ID', 'Structure Type'], as_index = False)
    df_amp_counts = grouped_by_amp['num_Reads'].agg({'num_Reads': 'sum'})

    df_amp_counts = df_amp_counts.set_index(['Amplicon ID', 'Structure Type'])['num_Reads'].unstack().reset_index()
    df_amp_counts.columns = df_amp_counts.columns.tolist()

    # if len(df_amp_counts) == 0:
    #     return pd.DataFrame()

    for m in mut_types:
        if not m in df_amp_counts.columns:
            df_amp_counts[m] = np.NaN
            df_amp_counts[m] = df_amp_counts[m].apply(pd.to_numeric, errors='ignore')

    # find the winner structure type for each edge based on majority of supporting reads

    dfamplicon = pd.merge(dfamplicon, df_amp_counts[['Amplicon ID'] + mut_types], how='left', on='Amplicon ID', copy=False)
    dfamplicon['Winner_Amplicon_Structure'] = dfamplicon[mut_types].idxmax(axis=1)

    return dfamplicon

def assignSegmentStructureType(dfsv):
    if dfsv.empty:
        return pd.DataFrame()

    # create a dataframe where each row contains one edge of the BP
    df_bp1 = pd.DataFrame(data={'Segment ID': dfsv['Segment1 ID'], 'num_Reads': dfsv['num_Reads'], 'Structure Type': dfsv['Structure Type'], 'orig_ind':dfsv.index})
    df_bp2 = pd.DataFrame(data={'Segment ID': dfsv['Segment2 ID'], 'num_Reads': dfsv['num_Reads'], 'Structure Type': dfsv['Structure Type'],'orig_ind': dfsv.index})

    df_bps = pd.concat([df_bp1, df_bp2])
    # workarround what seems to be a bug in drop_duplicates
    df_bps["joined"] = df_bps["Segment ID"].map(str) + "" + df_bps["orig_ind"].map(str)
    df_bps.reset_index(drop=True, inplace=True)

    # if two sides of bp are the same segment, remove one of it
    #df_bps.drop_duplicates(subset=['Segment ID', 'orig_ind'], inplace=True)
    df_bps.drop_duplicates(subset='joined', inplace=True)

    grouped_by_segment = df_bps.groupby(['Segment ID', 'Structure Type'], as_index = False)
    dfsegment_counts = grouped_by_segment['num_Reads'].agg({'num_Reads': 'sum'})

    dfsegment_counts = dfsegment_counts.set_index(['Segment ID', 'Structure Type'])['num_Reads'].unstack().reset_index()
    dfsegment_counts.columns = dfsegment_counts.columns.tolist()

    # if len(df_amp_counts) == 0:
    #     return pd.DataFrame()

    for m in mut_types:
        if not m in dfsegment_counts.columns:
            dfsegment_counts[m] = np.NaN
            dfsegment_counts[m] = dfsegment_counts[m].apply(pd.to_numeric, errors='ignore')

    # find the winner structure type for each edge based on majority of supporting reads

    dfsegment = dfsegment_counts[['Segment ID'] + mut_types]
    dfsegment.loc[:,'Winner_Segment_Structure'] = dfsegment[mut_types].idxmax(axis=1)

    return dfsegment



def run(excel, remove_unmappble, unify_seg, config_file, sample_id, cnv_file, sv_file, fast_output_dir, focused_genes):
    '''
    main run
    :param excel: flag: should we output to excel file
    :param remove_unmappble: flag: should we remove unmappble regions
    :param config_file: configuration file
    :param sample_id: sample ID
    :param cnv_file: path to cnv file
    :param sv_file: path to sv file
    :param fast_output_dir: output directory
    :param focused_genes: focused genes list
    :return: --
    '''

    if not os.path.exists(fast_output_dir):
        os.makedirs(fast_output_dir)

    # Read Data
    dfcnv = readCNVData(cnv_file, remove_unmappble, unify_seg)
    dfsv = readSVData(sv_file)

    dfsv = assignSegmentsForSVs(dfcnv, dfsv)
    # dfsv, dfsegments = assignSegmentsForSVs(dfcnv, dfsv)

    # Amplicons Handing
    dfcnv, dfamplicons = createAmpliconFromConnectedSegments(dfcnv['Segment ID'].max(), dfsv, dfcnv)
    dfcnv, dfsv = updateAmpliconsInSVs(dfcnv, dfsv)

    # Define Structure Type
    dfsv = setStructureTypeToSV(dfsv)
    dfamplicons = assignAmpliconStructureType(dfamplicons, dfsv)

    dfsegments = assignSegmentStructureType(dfsv)

    # focused genes

    df_focused_genes = segments2FocusedGenes(dfcnv, dfamplicons, dfsegments, focused_genes)
    if len(dfamplicons)>0:
        df_amp_gene = pd.merge(dfamplicons, df_focused_genes[['Amplicon ID', 'Gene']], how='left', left_on='Amplicon ID', right_on='Amplicon ID', copy=False)
        temp = df_amp_gene.groupby('Amplicon ID')['Gene'].apply(lambda g: ', '.join([x for x in g if isinstance(x,str)])).to_frame().reset_index()
        dfamplicons = pd.merge(dfamplicons, temp, how='left', on='Amplicon ID', copy=False)


    dfsegment_gene = pd.merge(dfcnv, df_focused_genes[['Segment ID', 'Gene']], how='left', left_on='Segment ID', right_on='Segment ID', copy=False)
    temp = dfsegment_gene.groupby('Segment ID')['Gene'].apply(lambda g: ', '.join([x for x in g if isinstance(x,str)])).to_frame().reset_index()
    dfcnv = pd.merge(dfcnv, temp, how='left', on='Segment ID', copy=False)

    # write to csv files
    df_focused_genes.to_csv(open(os.path.join(fast_output_dir, sample_id + '_focused_genes_FAST.csv'), 'w'))
    dfamplicons.to_csv(open(os.path.join(fast_output_dir, sample_id + '_amplicons_FAST.csv'),'w'))
    dfsegments.to_csv(open(os.path.join(fast_output_dir, sample_id + '_segments_FAST.csv'), 'w'))
    dfcnv.to_csv(open(os.path.join(fast_output_dir, sample_id + '_cnv_FAST.csv'),'w'))
    dfsv.to_csv(open(os.path.join(fast_output_dir, sample_id + '_sv_FAST.csv'),'w'))
    shutil.copyfile(config_file,os.path.join(fast_output_dir, sample_id + '_config.txt'))

    if excel:
        write_excel(fast_output_dir, sample_id, dfcnv, dfsv, dfamplicons, dfsegments, df_focused_genes)

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Analyze amplicon structure base on CNV and SV data.')
    parser.add_argument('--config', '-c', type=str, required = True)
    parser.add_argument('--excel', '-e', action='store_true', default = False, help = 'generated colored excel file')
    parser.add_argument('--remove_unmappble', '-r', action='store_true', default = False, help = 'remove unmappable regions')
    parser.add_argument('--not_unify_seg', '-n', action='store_true', default=False, help='do not unify consecutive segments')

    args = parser.parse_args()

    CP.readfp(open(args.config))

    timestr = time.strftime("%Y%m%d_%H%M%S")
    log = 'Fast_log_%s.csv' % timestr
    f_log = open(os.path.join(CP.get('FILES', 'FAST_OUTPUT_DIR'), log),'w')

    focused_locuses = [FL(*fl) for fl in csv.reader(open(CP.get('FILES','GENES')))]

    sample_id = CP.get('FILES','SAMPLE_ID')
    run(args.excel, args.remove_unmappble, not (args.not_unify_seg), args.config, sample_id, CP.get('FILES', 'CNV_INP'),
        CP.get('FILES', 'SV_INP'), CP.get('FILES', 'FAST_OUTPUT_DIR'), focused_locuses)
    # try:
    #     run(args.excel, args.remove_unmappble, not (args.not_unify_seg), args.config, sample_id, CP.get('FILES', 'CNV_INP'), CP.get('FILES', 'SV_INP'), CP.get('FILES', 'FAST_OUTPUT_DIR'), focused_locuses)
    # except Exception as E:
    #     print >> f_log, '%s, %s, FAILED: %s' % (datetime.datetime.now(), sample_id, E)
    # else:
    #     print >> f_log, '%s, %s, SUCCESS' % (datetime.datetime.now(), sample_id)
    #
    # f_log.close()
