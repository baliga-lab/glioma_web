#!/usr/bin/env python
import sys
import traceback as tb
import logging
import json
import os
from functools import wraps

import MySQLdb
from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask

import gzip
import pandas
import numpy as np

convert = {'Evading apoptosis':'cellDeath.gif', 'Evading immune detection':'avoidImmuneDestruction.gif', 'Genome instability and mutation':'genomicInstability.gif', 'Insensitivity to antigrowth signals':'evadeGrowthSuppressors.gif', 'Limitless replicative potential':'immortality.gif', 'Reprogramming energy metabolism':'cellularEnergetics.gif', 'Self sufficiency in growth signals':'sustainedProliferativeSignalling.gif', 'Sustained angiogenesis':'angiogenesis.gif', 'Tissue invasion and metastasis':'invasion.gif', 'Tumor promoting inflammation':'promotingInflammation.gif'}

app = Flask(__name__)
app.config.from_envvar('GLIOMA_SETTINGS')

######################################################################
#### General helpers
######################################################################

def dbconn():
    return MySQLdb.connect(host=app.config['HOST'], user=app.config['USER'],
                           passwd=app.config['PASS'], db=app.config['DB'])


def read_exps():
    with gzip.open(app.config['GENE_EXPR_FILE'], 'rb') as f:
        return pandas.read_csv(f, sep=',', index_col=0, header=0)

######################################################################
#### Boxplot functionality
######################################################################

def submat_data(submat, col_indexes):
    """given a sub matrix and a list of column indexes
    that specify the columns, of the matrix, return a list
    of (col_idx, median, min, max, lower_quartile, upper_quartile)
    tuples
    """
    col_medians = np.median(submat, axis=0)
    col_mins = np.min(submat, axis=0)
    col_maxs = np.max(submat, axis=0)
    col_upper_quarts = np.percentile(submat, q=75.0, axis=0)
    col_lower_quarts = np.percentile(submat, q=25.0, axis=0)
    data = [[idx,
             col_mins[i],
             col_lower_quarts[i],
             col_medians[i],
             col_upper_quarts[i],
             col_maxs[i]
             ]
            for i, idx in enumerate(col_indexes)]
    return sorted(data, key=lambda x: x[1])


def cluster_data(cursor, cluster_id, df):
    patient_map = {name: index
                   for index, name in enumerate(df.columns.values)}
    gene_map = {name: index for index, name in enumerate(df.index)}

    cursor.execute("""select g.symbol from bic_gene bg
join gene g on bg.gene_id=g.id where bicluster_id=%s""", [cluster_id])
    genes = [row[0] for row in cursor.fetchall()]

    cursor.execute("""select name from bic_pat bp
join patient p on bp.patient_id=p.id where bicluster_id=%s""",
                   [cluster_id])
    patients = [row[0] for row in cursor.fetchall()]

    cursor.execute("""select name from patient where id not in
(select patient_id from bic_pat where bicluster_id=%s)""",
                   [cluster_id])
    excluded_patients = [row[0] for row in cursor.fetchall()]


    gene_indexes = sorted([gene_map[g] for g in genes])
    patient_indexes = sorted([patient_map[p] for p in patients])
    ex_patient_indexes = sorted([patient_map[p] for p in excluded_patients])

    submat = df.values[np.ix_(gene_indexes, patient_indexes)]
    in_data = submat_data(submat, patient_indexes)

    ex_submat = df.values[:,ex_patient_indexes]
    out_data = submat_data(ex_submat, ex_patient_indexes)
    return in_data, out_data

BOXPLOT_COLOR_MAP = {
    'control': 'green',
    'classical': 'black',
    'neural': 'lightBlue',
    'NA': 'grey',
    'g_cimp': 'brown',
    'proneural': 'red',
    'mesenchymal': 'orange'
}
######################################################################
#### Available application paths
######################################################################

@app.errorhandler(Exception)
def unhandled_exception(e):
    app.logger.exception(e)
    return render_template('unknown_error.html')

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/gene')
def gene():
    return render_template('gene.html')

@app.route('/patient')
def patient():
    return render_template('patient.html')

@app.route('/bicluster/<bicluster>')
def bicluster(bicluster=None):
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id,name,var_exp_fpc,var_exp_fpc_p_value,survival,survival_p_value
FROM bicluster WHERE name=%s""", [bicluster])
    bc_pk, bc_name, bc_varexp_fpc, bc_varexp_fpc_pval, bc_survival, bc_survival_pval = c.fetchone()
    bic_info = {
        'pk': bc_pk,
        'name': bc_name,
        'varexp_fpc': bc_varexp_fpc,
        'varexp_fpc_pval': bc_varexp_fpc_pval,
        'survival': bc_survival,
        'survival_pval': bc_survival_pval,
        'varexp_flag': bc_varexp_fpc_pval <= 0.05,
        'survival_flag': bc_survival_pval <= 0.05
        }

    c.execute("""SELECT g.id, g.symbol, g.entrez FROM bic_gene bg join gene g on bg.gene_id=g.id where bg.bicluster_id=%s order by g.symbol""", [bc_pk])
    genes = list(c.fetchall())
    c.execute("""SELECT p.id, p.name FROM bic_pat bp join patient p on p.id=bp.patient_id where bp.bicluster_id=%s order by p.name""", [bc_pk])
    tumors = list(c.fetchall())
    # Replication
    c.execute("""SELECT * FROM replication WHERE bicluster_id=%s""", [bc_pk])
    tmp = list(c.fetchall())
    repConvert = {'French':'Gravendeel, et al. 2009','REMBRANDT':'Madhavan, et al. 2009','GSE7696':'Murat, et al. 2008'}
    repPubmed = {'French':'19920198','REMBRANDT':'19208739','GSE7696':'18565887'}
    replication = []

    replicated = [0, 0]
    bic_info['repl_coexp'] = False
    bic_info['repl_survival'] = False

    for i in tmp:
        tmp1 = [0,0]
        if bic_info['varexp_fpc_pval'] <= 0.05 and float(i[4])<=0.05:
            tmp1[0] = 1
            replicated[0] = 1
            bic_info['repl_coexp'] = True

        if (( bic_info['survival'] > 0 and (float(i[5]) > 0 and float(i[6])<=0.05)) or
            ( bic_info['survival'] < 0 and (float(i[5]) < 0 and float(i[6])<=0.05))):
            tmp1[1] = 1
            replicated[1] = 1
            bic_info['repl_survival'] = True

        replication.append(list(i)+[repConvert[i[2]], repPubmed[i[2]]]+tmp1)

    # Regulators
    regulators = []
    c.execute("""SELECT gene.id, gene.symbol, tf_regulator.action FROM tf_regulator, gene WHERE tf_regulator.bicluster_id=%s AND gene.id=tf_regulator.gene_id""",
              [bc_pk])
    tfs = list(c.fetchall())
    tfList = []
    for tf in tfs:
        known = 'No'
        c.execute("""SELECT * FROM tf_crispr WHERE gene_id=%s""", [tf[0]])
        for crispr in c.fetchall():
            if float(crispr[4])<=0.05:
                known = 'Yes'
        regulators.append(['TF', tf[0], tf[1], tf[2].capitalize(), known])
        tfList.append(tf[1])
    c.execute("""SELECT mirna.id, mirna.name, mirna.mir2disease, mirna.hmdd FROM mirna_regulator, mirna WHERE mirna_regulator.bicluster_id=%s AND mirna.id=mirna_regulator.mirna_id""", [bc_pk])
    mirnas = list(c.fetchall())
    mirnaList = []
    for mirna in mirnas:
        if not mirna[0] in mirnaList:
            known = 'No'
            if (not mirna[2]=='no') or (not mirna[3]==0):
                known = 'Yes'
            regulators.append(['miRNA', mirna[0], mirna[1], 'Repressor', known])
            mirnaList.append(mirna[1])
    regulators = sorted(regulators, key=lambda name: name[1])
    # Get causal flows with bicluster
    c.execute("""SELECT * FROM causal_flow WHERE bicluster_id=%s""", [bc_pk])
    tmp_cf = c.fetchall()
    causalFlows = []
    for cf1 in tmp_cf:
        if cf1[3]=='tf':
            c.execute("""SELECT symbol FROM gene WHERE id=%s""", [cf1[2]])
            g1 = c.fetchall()[0][0]
        else:
            c.execute("""SELECT name FROM mirna WHERE id=%s""", [cf1[2]])
            g1 = c.fetchall()[0][0]
        if (cf1[3]=='tf' and g1 in tfList) or (cf1[3]=='mirna' and g1 in mirnaList):
            c.execute("""SELECT * FROM somatic_mutation WHERE id=%s""", [cf1[1]])
            m1 = c.fetchall()[0]
            if m1[2]=='gene':
                c.execute("""SELECT symbol FROM gene WHERE id=%s""", [m1[1]])
                mut = c.fetchall()[0][0]
            elif m1[2]=='pathway':
                c.execute("""SELECT name FROM nci_nature_pathway WHERE id=%s""", [m1[1]])
                mut = c.fetchall()[0][0]
            causalFlows.append([mut, g1])

    causalFlows = sorted(causalFlows, key=lambda mutation: mutation[0])

    # Hallmarks of Cancer
    c.execute("""SELECT hallmark.name FROM hallmark, bic_hal WHERE bic_hal.bicluster_id=%s AND hallmark.id=bic_hal.hallmark_id""", [bc_pk])
    h1 = list(c.fetchall())
    h2 = [[i[0],convert[i[0]]] for i in h1]

    # GO
    c.execute("""SELECT go_bp.id, go_bp.go_id, go_bp.name FROM bic_go, go_bp WHERE bic_go.bicluster_id=%s AND go_bp.id=bic_go.go_bp_id""", [bc_pk])
    tmps = list(c.fetchall())
    gobps = []
    for gobp in tmps:
        c.execute("""SELECT gene.symbol FROM go_gene, gene, bic_gene WHERE go_gene.go_bp_id=%s AND bic_gene.bicluster_id=%s AND go_gene.gene_id=gene.id AND gene.id=bic_gene.gene_id""", [gobp[0], bc_pk])
        tmp = c.fetchall()
        gobps.append(list(gobp)+[sorted(list(set([i[0] for i in tmp])))])

    exp_data = read_exps()
    in_data, out_data = cluster_data(c, bc_pk, exp_data)    
    ratios_mean = np.mean(exp_data.values)
    print "# in_data: ", len(in_data)
    hallmarks = h2
    all_boxplot_data = in_data + out_data
    patients = [exp_data.columns.values[item[0]] for item in all_boxplot_data]
    c.execute("""select pt.name from patient p join phenotypes pt on p.phenotype_id=pt.id where p.name in %s""",
              [patients])
    phenotypes = [row[0] for row in c.fetchall()]
    boxplot_colors = [BOXPLOT_COLOR_MAP[pt] for pt in phenotypes]
    js_boxplot_data = [item[1:] for item in all_boxplot_data]
    perc20 = len(in_data) / 5
    quintiles = [perc20 * i for i in range(1, 6)]
    return render_template('bicluster.html', **locals())


@app.route('/search')
def search():
    gene = request.args.get('gene')
    db = dbconn()
    c = db.cursor()
    type ='gene'
    if not gene:
        return render_template('index.html')
    if gene.find('hsa-')==-1:
        c.execute("""SELECT * FROM gene WHERE symbol=%s""", [gene])
        geneData = c.fetchall()
    else:
        c.execute("""SELECT * FROM mirna WHERE name=%s""", [gene])
        geneData = c.fetchall()
        type = 'mirna'
    if len(geneData)==0:
        return render_template('index.html')
    else:
        # Get causal flows downstream of mutation in gene
        geneData = geneData[0]
        muts = {}
        c.execute("""SELECT * FROM somatic_mutation WHERE mutation_type='gene' AND ext_id=%s""", [geneData[0]])
        tmp_muts = c.fetchall()
        if len(tmp_muts)==1:
            muts['name'] = gene
            c.execute("""SELECT * FROM causal_flow WHERE somatic_mutation_id=%s""", [tmp_muts[0][0]])
            tmp_cf = c.fetchall()
            muts['flows'] = 0
            muts['regs'] = []
            muts['tfs'] = []
            muts['miRNAs'] = []
            muts['biclusters'] = []
            muts['data'] = []
            for cf1 in tmp_cf:
                g1 = ''
                if cf1[3]=='tf':
                    c.execute("""SELECT * FROM tf_regulator WHERE gene_id=%s AND bicluster_id=%s""", [cf1[2], cf1[4]])
                    if len(c.fetchall())>0:
                        c.execute("""SELECT symbol FROM gene WHERE id=%s""", [cf1[2]])
                        g1 = c.fetchall()[0][0]
                        if not g1 in muts['regs']:
                            muts['regs'].append(g1)
                            muts['tfs'].append(g1)
                else:
                    c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s AND bicluster_id=%s""", [cf1[2], cf1[4]])
                    if len(c.fetchall())>0:
                        c.execute("""SELECT name FROM mirna WHERE id=%s""", [cf1[2]])
                        g1 = c.fetchall()[0][0]
                        if not g1 in muts['regs']:
                            muts['regs'].append(g1)
                            muts['miRNAs'].append(g1)
                if not g1=='':
                    c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [cf1[4]])
                    b1 = c.fetchall()[0]
                    if not b1 in muts['biclusters']:
                        muts['biclusters'].append(b1[0])
                    c.execute("""SELECT hallmark.name FROM hallmark, bic_hal WHERE bic_hal.bicluster_id=%s AND hallmark.id=bic_hal.hallmark_id""", [cf1[4]])
                    tmp1 = c.fetchall()
                    h1 = list(set([convert[i[0]] for i in tmp1]))
                    h2 = [[i[0],convert[i[0]]] for i in tmp1]
                    muts['data'].append([gene, g1, b1[0], b1[1], b1[2], h2])

        # Get biclusters regulated by gene
        regs = {}
        if type=='gene':
            c.execute("""SELECT * FROM tf_regulator WHERE gene_id=%s""", [geneData[0]])
        else:
            c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s""", [geneData[0]])
        tmp_regs = c.fetchall()
        if len(tmp_regs)>0:
            regs['name'] = gene
            regs['biclusters'] = len(set([i[1] for i in tmp_regs]))
            regs['data'] = []
            # Collect all biclusters downstream regulated by TF or miRNA
            for reg in tmp_regs:
                action = 'Rep.'
                if type=='gene' and reg[3]=='activator':
                    action = 'Act.'
                c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [reg[1]])
                b1 = c.fetchall()[0]
                c.execute("""SELECT hallmark.name FROM hallmark, bic_hal WHERE bic_hal.bicluster_id=%s AND hallmark.id=bic_hal.hallmark_id""", [reg[1]])
                tmp1 = c.fetchall()
                h1 = list(set([convert[i[0]] for i in tmp1]))
                h2 = [[i[0],convert[i[0]]] for i in tmp1]
                regs['data'].append([gene, action, b1[0], b1[1], b1[2], h2])

        # Get biclusters that gene resides
        bics = {}
        if type=='gene':
            c.execute("""SELECT * FROM bic_gene, bicluster WHERE bic_gene.gene_id=%s AND bic_gene.bicluster_id=bicluster.id""", [geneData[0]])
            tmp_bics = c.fetchall()
            if len(tmp_bics)>0:
                bics['name'] = gene
                bics['biclusters'] = len(tmp_bics)
                bics['data'] = []
                for bic1 in tmp_bics:
                    c.execute("""SELECT hallmark.name FROM bic_hal, hallmark WHERE bic_hal.bicluster_id=%s AND bic_hal.hallmark_id=hallmark.id""", [bic1[3]])
                    tmp1 = c.fetchall()
                    h1 = list(set([convert[i[0]] for i in tmp1]))
                    h2 = [[i[0],convert[i[0]]] for i in tmp1]
                    bics['data'].append([bic1[4], bic1[5], bic1[6], bic1[7], bic1[8], h2])
        return render_template('search.html', gene=gene, muts=muts, regs=regs, bics=bics)

@app.route('/network')
def network():
    return render_template('network.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/download')
def download():
    return render_template('download.html')

@app.route('/citation')
def citation():
    return render_template('citation.html')

@app.route('/genecompletions')
def genecompletions():
    term = request.args.get('term')
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT symbol FROM gene WHERE symbol LIKE %s""", [str(term)+'%'])
    json1 = json.dumps([i[0] for i in c.fetchall()])
    return Response(response=json1, status=200, mimetype='application/json')

if __name__ == '__main__':
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    app.debug = True
    app.secret_key = 'supercalifragilistic'
    app.logger.addHandler(handler)
    app.run(host='0.0.0.0', debug=True)

