#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 10:06:15 2025

@author: keskusa2
"""


from collections import defaultdict, Counter
import networkx as nx
import plotly
import plotly.graph_objects as go
import numpy as np
import os

def get_clusters(cnas):
    LEN_THR = 5000000
    CLUST_LEN = 5000000
    CLUST_FREQ = 1000000
    MIN_SV = 2
    
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]
    
    by_sv = defaultdict(list)
    g = nx.MultiGraph()
    node_ids = {}
    id_to_kmers = {}
    id_to_cn = {}
    node_to_sv = defaultdict(list)
    for cns in cnas.values():
        for cn in cns:
            node_id = node_to_id(cn.get_name())
            id_to_cn[node_id] = cn
            g.add_node(node_id, _cn = cn.cn)
            if cn.sv1:
                by_sv[cn.sv1].append(node_id)
            if cn.sv2:
                by_sv[cn.sv2].append(node_id)
    clusters = []
    for cns in cnas.values():
        cur_cluster = []
        for cn in cns:
            if cn.pos_2 - cn.pos_1 > LEN_THR:
                continue
            if cur_cluster and cn.pos_1 - cur_cluster[-1].pos_2 > CLUST_LEN:
                clust_span = cur_cluster[-1].pos_2 - cur_cluster[0].pos_1
                if len(cur_cluster) > 3 and clust_span/len(cur_cluster)  <= CLUST_FREQ:
                    clusters.append(cur_cluster)
                cur_cluster = [cn]
            else:
                cur_cluster.append(cn)
        if cur_cluster:
            clust_span = cur_cluster[-1].pos_2 - cur_cluster[0].pos_1
            if len(cur_cluster) > 3 and clust_span/len(cur_cluster)  <= CLUST_FREQ:
                clusters.append(cur_cluster)
    for cl in clusters:
        for a,b in zip(cl[:-1], cl[1:]):
            node_id1 = node_to_id(a.get_name())
            id_to_cn[node_id1] = a
            node_id2 = node_to_id(b.get_name())
            g.add_edge(node_id1, node_id2, _type = 'adj', _supp = 0,  _jun_type = '' , _dtype ='', _sv = '')
            id_to_cn[node_id2] = b
    svs_to_cl = defaultdict(list)
    posls = []
    for i, cc in enumerate(nx.connected_components(g)):
        posls.append((id_to_cn[list(cc)[0]].ref_id, id_to_cn[list(cc)[0]].haplotype))
        if len(cc) < 2:
            continue
        for sv, nodes in by_sv.items():
            if len(nodes) < 2:
                continue
            ncs = [x for x in nodes if x in cc]
            for n in ncs:
                svs_to_cl[sv].append(i)
    cc = Counter([(c[0], c[1]) for c in svs_to_cl.values() if len(c) > 1])
    svls = []
    for key, val in cc.items():
        if val < 1:
            continue
        if posls[key[0]][0] == posls[key[1]][0] and not posls[key[0]][1] == posls[key[1]][1]:
            continue
        svls += [sv for sv, cl in svs_to_cl.items() if cl == [key[0], key[1]] or cl == [key[1], key[0]]]
    for sv, nodes in by_sv.items():
        if len(nodes) < 2:
            continue
        if not sv in svls:
            continue
        g.add_edge(nodes[0], nodes[1], _type = 'junction', _supp = sv.supp, _jun_type = sv.get_jun(), _dtype = sv.detailed_type, _sv = sv)
        node_to_sv[nodes[0]].append(sv)
        node_to_sv[nodes[1]].append(sv)
    all_segments = []
    for i, cc in enumerate(nx.connected_components(g)):
        if len(cc) < 2:
            continue
        segments = []
        cnls = []
        ref_ids = []
        juntype = []
        d_type = []
        jun_svs = []
        junctions = defaultdict(list)
        for c in cc:
            segments.append([id_to_cn[c].ref_id, id_to_cn[c].pos_1, id_to_cn[c].pos_2, id_to_cn[c].cn, c, id_to_cn[c].haplotype])
            cnls.append(id_to_cn[c].cn)
            ref_ids.append(id_to_cn[c].ref_id)
        for u, v, key, value in g.edges(cc, keys=True, data = True):
            juntype.append(value['_jun_type'])
            if value['_type'] == 'junction':
                junctions[(u,v)].append(value)
                d_type.append(value['_dtype'])
                jun_svs.append(value['_sv'])
        c_juntype = Counter(juntype)
        c_dtype = Counter(d_type)
        c_refids = Counter(ref_ids)
        if len(junctions) >= MIN_SV:
            all_segments.append([segments, junctions,c_juntype, c_dtype, c_refids])
    return all_segments
     
def html_plot(segments, junctions, subgr_num, out_dir):

    AXIS_OFF = 500000
    TOL=1000000
    segment_list = defaultdict(list)
    for seg in segments:
        segment_list[seg[0]].append(seg[1:]) 
    segment_list = dict(sorted(segment_list.items()))
    seg_bound = TOL
    chr_bound = [0]
    xticktext = list(segment_list.keys())
    for seg_ls in segment_list.values():
        seg_ls.sort(key = lambda x:x[1])
        st =  seg_bound + TOL - seg_ls[0][0]
        for seg in seg_ls:
            seg.append(seg[0] + st)
            seg.append(seg[1] + st)
        seg_bound = seg_bound = seg_ls[-1][6] + TOL
        chr_bound +=[seg_bound -TOL, seg_bound]
    seg_dat = defaultdict(list)
    seg_dat[1] = [[],[],[]]
    seg_dat[2] = [[],[],[]]
    xtick = [int(np.median(chr_bound[i-1:i+1])) for i, v in enumerate(chr_bound) if i%2]  
    for ref , seg_ls in segment_list.items():
        for seg in seg_ls:
            hoverdata = f"{ref}:{seg[0]}-{seg[1]}<br>Length:{seg[1]-seg[0]}"
            seg_dat[seg[4]][0] += [seg[5], seg[6], None]
            seg_dat[seg[4]][1] += [seg[2], seg[2], None]
            seg_dat[seg[4]][2]  += [hoverdata,hoverdata,None]
    colors = {'HH':"#256676",'TT':"#256676",'Foldback': "#a20655",'TH': "#4ea6dc",'HT': "#f19724", 'Interchr':'#cdcc50'}
    cols = {1:"#79AA3B", 2:"#6C3BAA"}
    x_nan = [s for xx in seg_dat.values() for s in xx[0] if s]
    x_limit = [min(x_nan) - AXIS_OFF, max(x_nan) + AXIS_OFF]
    y_limit = max([s for xx in seg_dat.values() for s in xx[1] if s])

    fig = go.Figure()
    add_chr_boxes(fig, chr_bound, segment_list.keys(),  y_limit)
    add_segments(fig, seg_dat, cols)
    #add_legend(fig,x[1][0], colors)
    add_dbs(fig, segment_list, junctions, colors, y_limit)
    plots_layout_settings(fig, list(segment_list.keys()), x_limit, subgr_num, xtick, xticktext)
    fig.write_html(out_dir + '/plots/padfoot_' + str(subgr_num) + ".html")



def add_legend(fig, x0, colors):
    fig.add_trace(go.Scatter(x=[x0], y=[1], legendgroup="Foldback", mode = 'lines',yaxis="y5",  
                             line = dict(shape = 'spline', color = colors['Foldback'], width= 7, dash = 'solid'),
                             name="Foldback"))
    fig.add_trace(go.Scatter(x=[x0], y=[1], legendgroup="TT/HH", mode = 'lines',  yaxis="y5",
                             line = dict(shape = 'spline', color = colors['TT'], width= 7, dash = 'solid'),
                             name="TT/HH"))
    fig.add_trace(go.Scatter(x=[x0], y=[1], legendgroup="TH", mode = 'lines', yaxis="y5", 
                             line = dict(shape = 'spline', color = colors['TH'], width= 7, dash = 'solid'),
                             name="TH"))
    fig.add_trace(go.Scatter(x=[x0], y=[1], legendgroup="HT", mode = 'lines', yaxis="y5", 
                             line = dict(shape = 'spline', color = colors['HT'], width= 7, dash = 'solid'),
                             name="HT"))
    fig.add_trace(go.Scatter(x=[x0], y=[1], legendgroup="Interchr", mode = 'lines', yaxis="y5", 
                             line = dict(shape = 'spline', color = colors['0-0'], width= 7, dash = 'solid'),
                             name="Interchr"))
    

def add_segments(fig, seg_dat, cols):
    for hp ,(x,y,hoverdata) in seg_dat.items():
        col = cols[hp]
        fig.add_trace(go.Scatter(
        x=x,
        y=y,
        name= 'HP' + str(hp),
        yaxis="y5",
        line = dict(shape = 'spline', color = col, width= 3, dash = 'solid'),
        mode='lines',
        opacity=0.9,
        text=hoverdata,
        hoverinfo="text"))
    
def add_chr_boxes(fig, chr_bound, chr_list, ylim):
    for i, ref in enumerate(chr_list):
        if i%2:
            x0 = chr_bound[i*2]
            x1 = chr_bound[i*2+1]
            fig.add_shape(type="rect",
                          x0=x0, y0=0, x1=x1, y1=ylim,
                          line=dict(
                              color="RoyalBlue",
                              width=0),
                          fillcolor="LightSalmon", opacity=0.2)
  
def add_dbs(fig, segment_list, junctions, colors, ylim):

    segments = [s for seg in segment_list.values() for s in seg]
    x = defaultdict(list)
    y = defaultdict(list)
    labs = defaultdict(list)
    for nodes, jun in junctions.items():
        segls = [seg for seg in segments if seg[3] in nodes]
        sv = jun[0]['_sv']
        lab = sv.get_name() +  '<br>supp:' + str(sv.supp)
        
        c1 = 'Interchr'
        if sv.bp_1[0] == sv.bp_2[0]:
            c1 = jun[0]['_jun_type']
            
        if 'foldback' in jun[0]['_dtype']:
            lab += '<br>Foldback'
            c1 = 'Foldback'
            
        seg1 = [seg for seg in segls if abs(seg[0] - sv.bp_1[1]) < 2 or abs(seg[1] - sv.bp_1[1]) < 2][0]
        seg2 = [seg for seg in segls if abs(seg[0] - sv.bp_2[1]) < 2 or abs(seg[1] - sv.bp_2[1]) < 2][0]
        x0 = seg1[5] if abs(seg1[1] - sv.bp_1[1]) < 2 else seg1[6]
        x2 = seg2[5] if abs(seg1[1] - sv.bp_1[1]) < 2 else seg2[6]
        x[c1] += [x0,x0, None,x2,x2, None, x0,int(np.mean([x0,x2])), x2, None]
        y[c1] += [0, ylim, None,0, ylim, None, ylim, ylim+1, ylim, None]
        labs[c1] += [lab, lab, None,lab, lab, None,lab, lab, lab, None]
    for c1, x1 in x.items():
        y1 = y[c1]
        lab = labs[c1]
        col = colors[c1] 
        fig.add_trace(go.Scatter(
        x=x1,
        y=y1,
        name=c1,
        yaxis="y5",
        line = dict(shape = 'spline', color = col, width= 1, dash = 'solid'),
        mode='lines',
        opacity=0.9,
        text=lab,
        hoverinfo="text"))


def plots_layout_settings(fig, chr_list, x_limit, cluster_ind, xtick, xticktext):
    fig.update_layout(
        xaxis=dict(
            type="linear",
            showline=True,
            zeroline=True,
            tickmode = 'array',
            tickvals = xtick,
            ticktext = xticktext,
            linecolor = "dimgray",
            range=x_limit,
            tickfont={"color": "black", 'size':15}
        ),
        yaxis5=dict(
            linecolor="dimgray",
            tickmode = 'array',
            side="left",
            tickfont={"color": "black", 'size':15},
            ticks="outside",
            title=dict(text="", font={"color": "dimgray"}),
            type="linear",
            showline=True,
            zeroline=True,
        ))
    
    fig.update_layout(
        template="plotly_white",
        font_family="Helvetica"
    )
    
    fig.update_layout(legend=dict(
        orientation = 'h', xanchor = "center", x = 0.45, y= 1.2))
    
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    fig.update_xaxes(tick0=0.0, rangemode="nonnegative")
    fig.update_layout(legend={'itemsizing': 'constant'})
    fig.update_layout(font_family= "Helvetica")
    
    fig.update_layout(
        title={
            'text': 'subcluster - ' + str(cluster_ind),
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family = "Helvetica",
        font_color = "black",
        font_size = 15,
        title_font_family = "Helvetica",
        title_font_color = "black",
        legend_font_size = 15
    )
    
    if len(chr_list) < 3:
        height = 400
        width = 800
    elif len(chr_list) < 7:
        height = 600
        width = 1200
    else:
        height = 800
        width = 1200
    
    xlen = x_limit[1]-x_limit[0]
    if xlen > 100000000:
        width += 200
    elif xlen > 150000000:
        width +=400
     
    fig.update_layout(
        width=width,
        height=height,
       )
            
def out_findings(all_segments, out_dir):
    all_segments.sort(key = lambda s:len(s[0]), reverse = True)
    out_folderhtml = os.path.join(out_dir , 'plots')
    if not os.path.isdir(out_folderhtml):
        os.mkdir(out_folderhtml)
    #out_summary = open(out_dir + '/summary_clusters.tsv')
    for subgr_num, (segments, junctions,c_juntype, c_dtype, c_refids) in enumerate(all_segments):
        html_plot(segments, junctions, subgr_num, out_dir)
        
        
def cluster_cn(cnas, out_dir):
    all_segments = get_clusters(cnas)
    out_findings(all_segments, out_dir)
      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
        
    
    

        
