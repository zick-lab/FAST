import argparse
import ConfigParser

from excelUtils import *

class FL:
    def __init__(self,name, chrm, start,end):
        self.name = name
        self.chrm = chrm
        self.start = int(start)
        self.end = int(end)


def analyzeFileFocusedLocuses(sample_id, df_cnv, df_amplicon, focused_genes):
    '''
    Analyze whether and which segments contain specific genes
    :param sample_id: sample ID
    :param df_cnv: dataframe of cnvs
    :param df_segment: dataframe of segments
    :param df_amplicon: dataframe of amplicons
    :param focused_genes: list of focused genes
    :return:
    '''

    l_res = []
    for gene in focused_genes:
        df_cnv_gene = df_cnv[(df_cnv.chr==gene.chrm) & (df_cnv.start <= gene.end) & (df_cnv.end >= gene.start)]
        if len(df_cnv_gene)>0:
            segment_id = df_cnv_gene.iloc[0].segment
            copy_number = df_cnv_gene.iloc[0].CopyNumber
            # segment_data = df_segment[df_segment.segment==segment_id].iloc[0]
            amplicon = df_cnv_gene.iloc[0].amplicon
            # amplicon_fl = segment_data['amplicon']
            # amplicon_data = df_amplicon[df_amplicon.amplicon==amplicon_fg].iloc[0]
            # segment_var_type = segment_data['Winner_segment_structure']
            # amplicon_var_type = amplicon_data['Winner_amplicon_structure']
            amplicon_struct_type = df_amplicon[df_amplicon.amplicon==amplicon].iloc[0].Winner_Amplicon_Structure
            l_res.append([sample_id, gene.name, amplicon_struct_type, copy_number, amplicon, segment_id])
    if len(l_res)==0:
        print 'No segment in oncogenes: ', sample_id
    return l_res

def segments2FocusedGenes(s, df_cnv, df_amplicon, focused_locuses):
    '''
    :param s: sample ID
    :param df_cnv: dataframe of cnvs
    :param df_amplicon: dataframe of amplicons
    :param focused_locuses: list of focused genes
    :return:
    '''

    l_focusedGenes = analyzeFileFocusedLocuses(s, df_cnv, df_amplicon, focused_locuses)

    if len(l_focusedGenes)==0:
        df_focused_genes_analysis = pd.DataFrame(columns=['sample_id','Gene','Amplicon Structure', 'amplicon_id', 'segment_id'], data=[[s, None, None, None, None]])
    else:
        df_focused_genes_analysis = pd.DataFrame(l_focusedGenes, columns=['sample_id','Gene','Amplicon Structure','copy_number','amplicon_id', 'segment_id'])
        df_focused_genes_analysis.sort_values(by='Amplicon Structure',ascending=True, inplace=True, na_position='last')
        del df_focused_genes_analysis['sample_id']

    return  df_focused_genes_analysis

