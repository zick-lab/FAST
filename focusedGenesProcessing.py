import argparse
import ConfigParser

from excelUtils import *

class FL:
    def __init__(self,name, chrm, start,end):
        self.name = name
        self.chrm = chrm
        self.start = int(start)
        self.end = int(end)


def analyzeFileFocusedLocuses(dfcnv, dfamplicon, dfsegment, focused_genes):
    '''
    Analyze whether and which segments contain specific genes
    :param sample_id: sample ID
    :param dfcnv: dataframe of cnvs
    :param dfsegment: dataframe of segments
    :param dfamplicon: dataframe of amplicons
    :param focused_genes: list of focused genes
    :return:
    '''

    l_res = []

    for gene in focused_genes:
        dfcnv_gene = dfcnv[(dfcnv.chr==gene.chrm) & (dfcnv.start <= gene.end) & (dfcnv.end >= gene.start)]
        if len(dfcnv_gene)>0:
            copy_number = dfcnv_gene.iloc[0].CopyNumber
            amplicon = dfcnv_gene.iloc[0]['Amplicon ID']
            segment = dfcnv_gene.iloc[0]['Segment ID']
            if len(dfamplicon)==0:
                amplicon_struct_type = None
            else:
                amplicon_struct_type = dfamplicon[dfamplicon['Amplicon ID']==amplicon].iloc[0].Winner_Amplicon_Structure
            if len(dfsegment)==0:
                segment_struct_type = None
            else:
                df_seg = dfsegment[dfsegment['Segment ID']==segment]
                if len(df_seg)>0:
                    segment_struct_type = df_seg.iloc[0].Winner_Segment_Structure
                else:
                    segment_struct_type = None

            l_res.append([gene.name, copy_number, amplicon, amplicon_struct_type, segment, segment_struct_type])
    if len(l_res)==0:
        print 'No segment in oncogenes:'
    return l_res

def segments2FocusedGenes(dfcnv, dfamplicon, dfsegment, focused_locuses):
    '''
    :param s: sample ID
    :param dfcnv: dataframe of cnvs
    :param dfamplicon: dataframe of amplicons
    :param focused_locuses: list of focused genes
    :return:
    '''

    l_focusedGenes = analyzeFileFocusedLocuses(dfcnv, dfamplicon, dfsegment, focused_locuses)

    focused_genes_cols = ['Gene', 'copy_number', 'Amplicon ID', 'Amplicon Structure', 'Segment ID', 'Segment Structure']
    if len(l_focusedGenes)==0:
        df_focused_genes_analysis = pd.DataFrame(columns=focused_genes_cols, data=[[None]*len(focused_genes_cols)])
    else:
        df_focused_genes_analysis = pd.DataFrame(l_focusedGenes, columns=focused_genes_cols)
        df_focused_genes_analysis.sort_values(by='Amplicon Structure',ascending=True, inplace=True, na_position='last')

    return  df_focused_genes_analysis

