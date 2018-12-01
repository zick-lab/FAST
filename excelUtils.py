import os
import xlsxwriter
import pandas as pd
import numpy as np

COLORS = [('yellow','black'),
           ('blue', 'white'),
           ('brown','white'),
           ('cyan', 'black'),
           ('gray','black'),
           ('green','black'),
           ('lime', 'black'),
           ('orange','black'),
           ('pink','black'),
           ('purple','black'),
           ('red','black'),
           ('silver','black'),
           ('navy','white'),
           ('black','white')]



def setFormats(workbook):
    formats = {}

    title_frmt = workbook.add_format()
    title_frmt.set_bold()
    title_frmt.set_font_size(16)
    title_frmt.set_underline()
    formats['title'] = title_frmt

    header = workbook.add_format()
    header.set_bold()
    formats['header'] = header

    colored = []
    bordered = []
    for c in COLORS:
        frmt = workbook.add_format()
        frmt.set_pattern(1)
        frmt.set_bg_color(c[0])
        frmt.set_font_color(c[1])
        colored.append(frmt)

        frmt = workbook.add_format()
        frmt.set_pattern(1)
        frmt.set_bg_color(c[0])
        frmt.set_font_color(c[1])
        frmt.set_border(3)
        frmt.set_pattern(15)
        frmt.set_bold()
        bordered.append(frmt)

    formats['colored'] = colored
    formats['bordered'] = bordered
    return formats

def writeHeader(title, l,worksheet, OFFSET, formats, comments=None):
    worksheet.write(OFFSET, 0, title, formats['title'])
    OFFSET += 1
    if comments:
        for comment in comments:
            print comment
            worksheet.write(OFFSET, 0, comment + '\n')
            OFFSET+=1
    if len(l)==0:
        worksheet.write(OFFSET, 0, 'No data for this table', formats['header'])
    for i,item in enumerate(l):
        if str(item)=='nan':
            worksheet.write_blank(OFFSET,i,None)
        else:
            worksheet.write(OFFSET, i, l[i], formats['header'])
    OFFSET += 1
    return OFFSET

def writeEmptyHeader(title, message, worksheet, OFFSET, formats):
    worksheet.write(OFFSET, 0, title, formats['title'])
    OFFSET += 1
    worksheet.write(OFFSET, 0, message)
    OFFSET += 1
    return OFFSET

def writeEmptyRows(worksheet, OFFSET, n):
    for i in range(n):
        worksheet.write_blank(OFFSET+i,0,None)
    OFFSET += (i+1)
    return OFFSET

def setColors(l):
    return dict(zip(l,range(len(l)))), len(l)

def write_simple_table(worksheet, name, df, OFFSET,formats, index=False, comments=None):
    if df is not None and len(df)>0:
        cols = list(df.columns.values)
        cols_header = [df.index.name]+cols if index else cols
    else:
        cols = []
        cols_header = []
    OFFSET = writeHeader(name, cols_header, worksheet, OFFSET, formats,comments)
    if df is not None:
        iir = 0
        for ir, row in df.iterrows():
            if index:
                worksheet.write(OFFSET+iir, 0, ir)
                index_offset = 1
            else:
                index_offset = 0
            for ic in range(len(cols)):
                if not pd.isnull(row[ic]):
                    worksheet.write(OFFSET+iir, ic+index_offset, row[ic])
                else:
                    worksheet.write_blank(OFFSET+iir,ic+index_offset,None)
            iir+=1
        OFFSET += iir
    OFFSET = writeEmptyRows(worksheet, OFFSET, 3)
    return OFFSET


def write_excel(fast_output_dir, sample_id, dfcnv, dfsv, dfamplicons, dfsegments, df_focused_genes):
    outputFileXls = os.path.join(fast_output_dir, sample_id + '_FAST.xlsx')
    workbook   = xlsxwriter.Workbook(outputFileXls)
    worksheet1 = workbook.add_worksheet('FAST')
    formats = setFormats(workbook)
    used_colors, ci = setColors(dfcnv['Segment ID'].unique())
    used_amplicons_colors , ci1= setColors(dfcnv['Amplicon ID'].unique())

    OFFSET = 0
    #
    # focused genes
    OFFSET = write_simple_table(worksheet1, 'Focused Genes', df_focused_genes, OFFSET, formats)

    OFFSET = writeEmptyRows(worksheet1, OFFSET, 3)

    # write dfamplicons data
    cols = [] if dfamplicons is None else list(dfamplicons.columns.values)
    OFFSET = writeHeader('Amplicon Analysis ', cols, worksheet1, OFFSET, formats)
    iir = 0
    if dfamplicons is not None:
        for ir, row in dfamplicons.iterrows():
            ##            if np.isfinite(row['Amplicon ID']):
            frmt = formats['bordered'][used_amplicons_colors[row['Amplicon ID']]%len(COLORS)]
            worksheet1.write(OFFSET+iir, 0, row['Amplicon ID'],frmt)
            for ic in range(1,len(cols)):
                if not pd.isnull(row[ic]):
                    worksheet1.write(OFFSET+iir, ic, row[ic])
                else:
                    worksheet1.write_blank(OFFSET+iir,ic,None)
            iir+=1

        OFFSET += iir

    # write empty rows
    OFFSET = writeEmptyRows(worksheet1, OFFSET, 3)

    # write dfsegments data
    cols = [] if dfsegments is None else list(dfsegments.columns.values)
    OFFSET = writeHeader('Segment Analysis ', cols, worksheet1, OFFSET, formats)
    iir = 0
    if dfsegments is not None:
        for ir, row in dfsegments.iterrows():
            ##            if np.isfinite(row['Amplicon ID']):
            frmt = formats['colored'][used_colors[row['Segment ID']] % len(COLORS)]
            worksheet1.write(OFFSET+iir, 0, row['Segment ID'],frmt)
            for ic in range(1,len(cols)):
                if not pd.isnull(row[ic]):
                    worksheet1.write(OFFSET+iir, ic, row[ic])
                else:
                    worksheet1.write_blank(OFFSET+iir,ic,None)
            iir+=1

        OFFSET += iir

    # write empty rows
    OFFSET = writeEmptyRows(worksheet1, OFFSET, 3)


    # write CNV data with corresponding colors
    cols = [] if dfcnv is None else list(dfcnv.columns.values)
    OFFSET = writeHeader('CNV Data', cols, worksheet1, OFFSET, formats)
    iir = 0
    if dfcnv is not None:
        for ir, row in dfcnv.iterrows():
            for ic, col in enumerate(dfcnv.columns):
                if col=='Segment ID':
                    if np.isfinite(row['Segment ID']):
                        frmt = formats['colored'][used_colors[row['Segment ID']] % len(COLORS)]
                    else:
                        frmt = None
                elif col=='Amplicon ID':
                    if np.isfinite(row['Amplicon ID']):
                        frmt = formats['bordered'][used_amplicons_colors[row['Amplicon ID']] % len(COLORS)]
                    else:
                        frmt = None
                else:
                    frmt = None

                if not pd.isnull(row[ic]):
                    if frmt:
                        worksheet1.write(OFFSET+iir, ic, row[ic],frmt)
                    else:
                        worksheet1.write(OFFSET + iir, ic, row[ic])
                else:
                    worksheet1.write_blank(OFFSET+iir,ic,None)
            iir+=1
    OFFSET += iir


    # write empty rows
    OFFSET = writeEmptyRows(worksheet1, OFFSET, 3)

    # write SV Data
    cols = [] if dfsv is None else list(dfsv.columns.values)
    OFFSET = writeHeader('SV Data', cols, worksheet1, OFFSET, formats)
    iir = 0

    # write dfsv
    if dfsv is not None:
        for ir, row in dfsv.iterrows():
            if pd.isnull(row['Segment1 ID']) and pd.isnull(row['Segment2 ID']):
                continue
            frmt_list = [None, None]
            for i in range(2):
                if not pd.isnull(row['Segment%d ID' % (i+1)]):
                    frmt_list[i] = formats['colored'][used_colors[row['Segment%d ID'%(i+1)]]%len(COLORS)]

            for ic, col in enumerate(dfsv.columns):
                # one of the following columns: Chr1/2	Pos1/2	Orientation1/2 segment1/2   edge1/2
                if col.find('1') > 0:
                    frmt = frmt_list[0]
                elif col.find('2') > 0:
                    frmt = frmt_list[1]
                else:
                    if col=='Amplicon ID':
                        if np.isfinite(row['Amplicon ID']):
                            frmt = formats['bordered'][used_amplicons_colors[row['Amplicon ID']] % len(COLORS)]
                            worksheet1.write(OFFSET + iir, ic, row['Amplicon ID'], frmt)
                            continue
                    else: frmt = None

                if not pd.isnull(row[col]):
                    if frmt is None:
                        worksheet1.write(OFFSET + iir, ic, row[col])
                    else:
                        worksheet1.write(OFFSET + iir, ic, row[col], frmt)
                else:
                    worksheet1.write_blank(OFFSET + iir, ic, None)
            iir+=1

    workbook.close()
