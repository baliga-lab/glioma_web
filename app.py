#!/usr/bin/env python
import sys
import traceback as tb
import logging
import json
import os
import csv
from functools import wraps

import MySQLdb
from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask

import math
import itertools
import gzip
import pandas
import numpy as np
import rpy2.robjects as robjects


NUM_PARTS = 5

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
#### Graph/Visualization functionality
######################################################################

GRAPH_COLOR_MAP = {
    'control': '#6abd45',
    'classical': 'black',
    'neural': '#32689b',
    'NA': 'grey',
    'g_cimp': '#8a171a',
    'proneural': '#ed2024',
    'mesenchymal': '#faa41a'
}

# order in which the enrichment phenotypes are ordered
ENRICHMENT_PHENOTYPES = [
    'g_cimp', 'proneural', 'neural', 'classical', 'mesenchymal', 'control'
]


def phyper(q, m, n, k, lower_tail=False):
    """calls the R function phyper"""
    r_phyper = robjects.r['phyper']
    kwargs = {'lower.tail': lower_tail}
    return float(r_phyper(float(q), float(m), float(n), float(k), **kwargs)[0])


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
    return sorted(data, key=lambda x: x[3])


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

    ex_submat = df.values[np.ix_(gene_indexes, ex_patient_indexes)]
    out_data = submat_data(ex_submat, ex_patient_indexes)
    return in_data, out_data


def subtype_enrichment(cursor, cluster_id, df):
    patient_map = {name: index
                   for index, name in enumerate(df.columns.values)}
    gene_map = {name: index for index, name in enumerate(df.index)}

    cursor.execute("""select g.symbol from bic_gene bg
join gene g on bg.gene_id=g.id where bicluster_id=%s""", [cluster_id])
    genes = [row[0] for row in cursor.fetchall()]

    cursor.execute("""select name from bic_pat bp
join patient p on bp.patient_id=p.id where bicluster_id=%s""",
                   [cluster_id])
    included_patients = [row[0] for row in cursor.fetchall()]

    cursor.execute("""select name from patient where id not in
(select patient_id from bic_pat where bicluster_id=%s)""",
                   [cluster_id])
    excluded_patients = [row[0] for row in cursor.fetchall()]

    gene_indexes = sorted([gene_map[g] for g in genes])

    # above is exactly like boxplot
    # now make a phenotype map
    cursor.execute("""select p.name, pt.name from patient p join phenotypes pt on p.phenotype_id=pt.id
where pt.name <> 'NA'""")
    ptmap = {patient: phenotype for patient, phenotype in cursor.fetchall()}
    all_patients = {patient for patient in ptmap.keys()}
    phenotypes = {phenotype for phenotype in ptmap.values()}

    in_patient_indexes = sorted([patient_map[p] for p in included_patients if p in all_patients])
    ex_patient_indexes = sorted([patient_map[p] for p in excluded_patients if p in all_patients])

    # we use the submat_data function to sort our patients
    in_submat = df.values[np.ix_(gene_indexes, in_patient_indexes)]
    in_data = submat_data(in_submat, in_patient_indexes)
    sorted_in_indexes = [row[0] for row in in_data]
    ex_submat = df.values[np.ix_(gene_indexes, ex_patient_indexes)]
    ex_data = submat_data(ex_submat, ex_patient_indexes)
    sorted_ex_indexes = [row[0] for row in ex_data]

    # sorted by median pValue
    sorted_patient_indexes = sorted_in_indexes + sorted_ex_indexes
    
    # group patients into phenotype groups.
    # NOTE: the ptmap items need to be sorted, otherwise groupby fails to group correctly
    pt_patients = itertools.groupby(sorted(ptmap.items(), key=lambda pair: pair[1]),
                                    key=lambda pair: pair[1])
    pt_patients = {phenotype: set(map(lambda p: p[0], patients))
                   for phenotype, patients in pt_patients}

    num_columns = len(all_patients)
    cols_per_part = int(math.floor(num_columns / NUM_PARTS))
    #print "# columns: %d # cols/part: %d" % (num_columns, cols_per_part)

    pvalues = []
    min_pvalue = 100.0
    max_pvalue = -100.0

    for i in range(NUM_PARTS):
        part_pvalues = {}
        start = cols_per_part * i
        end = (cols_per_part * (i + 1)) - 1

        # adjust end for the last part
        if i == (NUM_PARTS - 1) and end != (num_columns - 1):
            end = num_columns - 1
        cur_patients = [df.columns.values[p_i] for p_i in sorted_patient_indexes[start:end + 1]]
        #print "Part %d, %d-%d, # current patients: %d" % (i, start, end, len(cur_patients))
        for phenotype in phenotypes:
            q = len([p for p in cur_patients if p in pt_patients[phenotype]])
            k = len(cur_patients)
            m = len(pt_patients[phenotype])
            n = num_columns - m
            pvalue = phyper(q, m, n, k)
            #print "part %d, phenotype: %s, q=%d, k=%d, m=%d, n=%d" % (i, phenotype, q, k, m, n)
            if pvalue == 0.0:
                pvalue = 10e-10

            if pvalue <= 0.5:
                pvalue = -math.log10(2 * pvalue)
            else:
                pvalue = math.log10(2 * (1.0 - pvalue))

            if math.isinf(pvalue):
                signum = -1 if pvalue < 0 else 1
                pvalue = signum * -math.log10(10e-10)
            
            if pvalue < min_pvalue:
                min_pvalue = pvalue
            if pvalue > max_pvalue:
                max_pvalue = pvalue

            part_pvalues[phenotype] = pvalue
        pvalues.append(part_pvalues)

    return pvalues, min_pvalue, max_pvalue


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
    elements = []
    elements.append({'data': { 'id': 'bc%d' % bc_pk, 'name': bc_name}, 'classes': 'bicluster' })

    regulators = []
    c.execute("""SELECT g.id, g.symbol, tfr.action FROM tf_regulator tfr join gene g on tfr.gene_id=g.id
WHERE tfr.bicluster_id=%s""", [bc_pk])
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
        elements.append({'data': { 'id': 'reg%d' % tf[0], 'name': tf[1] }, 'classes': 'tf' })
        elements.append({'data': { 'id': 'tfbc%d' % tf[0], 'source': 'reg%d' % tf[0], 'target': 'bc%d' % bc_pk }, 'classes': tf[2] })

    c.execute("""SELECT mirna.id, mirna.name, mirna.mir2disease, mirna.hmdd
FROM mirna_regulator mr join mirna on mirna.id=mr.mirna_id WHERE mr.bicluster_id=%s""", [bc_pk])
    mirnas = list(c.fetchall())

    mirnaList = []
    for mirna in mirnas:
        if not mirna[0] in mirnaList:
            known = 'No'
            if (not mirna[2]=='no') or (not mirna[3]==0):
                known = 'Yes'
            regulators.append(['miRNA', mirna[0], mirna[1], 'Repressor', known])
            mirnaList.append(mirna[1])
            elements.append({'data': { 'id': 'reg%d' % mirna[0], 'name': mirna[1]}, 'classes': 'mirna' })
            elements.append({'data': { 'id': 'mirnabc%d' % mirna[0], 'source': 'reg%d' % mirna[0], 'target': 'bc%d' % bc_pk }, 'classes': 'repressor' })

    regulators = sorted(regulators, key=lambda name: name[1])

    # Get causal flows with bicluster
    c.execute("""SELECT id,somatic_mutation_id,regulator_id,regulator_type,bicluster_id,leo_nb_atob,mlogp_m_atob
FROM causal_flow WHERE bicluster_id=%s""", [bc_pk])
    tmp_cf = c.fetchall()
    causalFlows = []
    

    for cf_pk, cf_som_mut_id, cf_reg_id, cf_reg_type, cf_bc_id, cf_leo, cf_mlogp in tmp_cf:
        if cf_reg_type == 'tf':
            c.execute("""SELECT symbol FROM gene WHERE id=%s""", [cf_reg_id])
            g1 = c.fetchone()[0]
        else:
            c.execute("""SELECT name FROM mirna WHERE id=%s""", [cf_reg_id])
            g1 = c.fetchone()[0]

        if (cf_reg_type == 'tf' and g1 in tfList) or (cf_reg_type == 'mirna' and g1 in mirnaList):
            c.execute("""SELECT * FROM somatic_mutation WHERE id=%s""", [cf_som_mut_id])
            m1 = c.fetchone()
            if m1[2]=='gene':
                c.execute("""SELECT symbol FROM gene WHERE id=%s""", [m1[1]])
                mut = c.fetchone()[0]
            elif m1[2]=='pathway':
                c.execute("""SELECT name FROM nci_nature_pathway WHERE id=%s""", [m1[1]])
                mut = c.fetchone()[0]
            causalFlows.append([mut, g1])
            elements.append({'data': { 'id': 'mut%d' % cf_som_mut_id, 'name': mut}, 'classes': 'genotype' })
            elements.append({'data': { 'id': 'cf%d' % cf_pk, 'source': 'mut%d' % cf_som_mut_id, 'target': 'reg%d' % cf_reg_id } })

    causalFlows = sorted(causalFlows, key=lambda mutation: mutation[0])

    # Hallmarks of Cancer
    c.execute("""SELECT hm.id,hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id
WHERE bh.bicluster_id=%s""", [bc_pk])
    hallmarks = []
    for hm_id, hm_name in c.fetchall():
        hallmarks.append([hm_name, convert[hm_name] ])
        elements.append({'data': { 'id': 'hm%d' % hm_id, 'name': hm_name}, 'classes': 'hallmark' })
        elements.append({'data': { 'id': 'bchm%d' % hm_id, 'source': 'bc%d' % bc_pk, 'target': 'hm%d' % hm_id } })

    # GO
    c.execute("""SELECT go_bp.id, go_bp.go_id, go_bp.name FROM bic_go join go_bp on go_bp.id=bic_go.go_bp_id
WHERE bic_go.bicluster_id=%s""", [bc_pk])
    tmps = list(c.fetchall())
    gobps = []
    for gobp in tmps:
        c.execute("""SELECT distinct gene.symbol FROM go_gene, gene, bic_gene WHERE go_gene.go_bp_id=%s AND bic_gene.bicluster_id=%s AND go_gene.gene_id=gene.id AND gene.id=bic_gene.gene_id order by gene.symbol""", [gobp[0], bc_pk])
        gobps.append(list(gobp) + [[row[0] for row in c.fetchall()]])

    # Prepare graph plotting data
    exp_data = read_exps()
    in_data, out_data = cluster_data(c, bc_pk, exp_data)
    enrichment_pvalues, min_enrichment_pvalue, max_enrichment_pvalue = subtype_enrichment(c, bc_pk, exp_data)
    js_enrichment_data = []
    js_enrichment_colors = []
    for part in enrichment_pvalues:
        for phenotype in ENRICHMENT_PHENOTYPES:
            js_enrichment_data.append([phenotype, part[phenotype]])
            js_enrichment_colors.append(GRAPH_COLOR_MAP[phenotype])
    enrichment_upper = -math.log10(0.05/30.0)
    enrichment_lower = math.log10(0.05/30.0)
    enrich_perc20 = len(js_enrichment_data) / 5
    enrich_quintiles = [enrich_perc20 * i for i in range(1, 6)]

    ratios_mean = np.mean(exp_data.values)
    all_boxplot_data = in_data + out_data
    patients = [exp_data.columns.values[item[0]] for item in all_boxplot_data]
    c.execute("""select p.name, pt.name from patient p join phenotypes pt on p.phenotype_id=pt.id where p.name in %s""",
              [patients])
    ptmap = {patient: phenotype for patient, phenotype in c.fetchall()}
    phenotypes = [ptmap[patient] for patient in patients]
    boxplot_colors = [GRAPH_COLOR_MAP[pt] for pt in phenotypes]
    js_boxplot_data = [[patients[i]] + item[1:] for i, item in enumerate(all_boxplot_data)]
    perc20 = len(in_data) / 5
    quintiles = [perc20 * i for i in range(1, 6)]
    db.close()

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
        c.execute("""SELECT symbol FROM gene WHERE symbol=%s""", [gene])
        geneData = c.fetchall()
    else:
        c.execute("""SELECT name FROM mirna WHERE name=%s""", [gene])
        geneData = c.fetchall()
        type = 'mirna'
    db.close()

    if len(geneData)==0:
        return render_template('index.html')

    symbol = geneData[0][0]
    if type == 'gene':
        return redirect(url_for('gene', symbol=symbol))
    else:
        return redirect(url_for('mirna', symbol=symbol))        


def __get_muts(c, gene_pk, symbol):
    # Get causal flows downstream of mutation in gene
    muts = {}
    c.execute("""SELECT * FROM somatic_mutation WHERE mutation_type='gene' AND ext_id=%s""", [gene_pk])
    tmp_muts = c.fetchall()
    if len(tmp_muts)==1:
        muts['name'] = symbol
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
                c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id  WHERE bh.bicluster_id=%s",
                          [cf1[4]])

                tmp1 = c.fetchall()
                h1 = list(set([convert[i[0]] for i in tmp1]))
                h2 = [[i[0],convert[i[0]]] for i in tmp1]
                muts['data'].append([symbol, g1, b1[0], b1[1], b1[2], h2])
    return muts


def __get_regulators(c, symbol):
    regs = {}
    tmp_regs = c.fetchall()
    if len(tmp_regs)>0:
        regs['name'] = symbol
        regs['biclusters'] = len(set([i[1] for i in tmp_regs]))
        regs['data'] = []
        # Collect all biclusters downstream regulated by TF or miRNA
        for reg in tmp_regs:
            action = 'Rep.'
            if type=='gene' and reg[3] == 'activator':
                action = 'Act.'
            c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [reg[1]])
            b1 = c.fetchall()[0]
            c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id WHERE bh.bicluster_id=%s",
                      [reg[1]])
            tmp1 = c.fetchall()
            h1 = list(set([convert[i[0]] for i in tmp1]))
            h2 = [[i[0],convert[i[0]]] for i in tmp1]
            regs['data'].append([symbol, action, b1[0], b1[1], b1[2], h2])
    return regs


@app.route('/mirna')
@app.route('/mirna/<symbol>')
def mirna(symbol=None, defaults={'symbol': None}):
    # Get biclusters regulated by mirna
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id FROM mirna WHERE name=%s""", [symbol])
    mirna_pk = c.fetchone()[0]
    muts = __get_muts(c, mirna_pk, symbol)
    c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s""", [mirna_pk])
    regs = __get_regulators(c, symbol)
    db.close()
    return render_template('search.html', gene=symbol, muts={}, regs=regs, bics={})


@app.route('/gene')
@app.route('/gene/<symbol>')
def gene(symbol=None, defaults={'symbol': None}):
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id FROM gene WHERE symbol=%s""", [symbol])
    gene_pk = c.fetchone()[0]
    muts = __get_muts(c, gene_pk, symbol)
    c.execute("""SELECT * FROM tf_regulator WHERE gene_id=%s""", [gene_pk])
    regs = __get_regulators(c, symbol)

    # Get biclusters that gene resides
    bics = {}
    c.execute("SELECT * FROM bic_gene bg join bicluster b on bg.bicluster_id=b.id where bg.gene_id=%s", [gene_pk])
    tmp_bics = c.fetchall()
    if len(tmp_bics) > 0:
        bics['name'] = gene
        bics['biclusters'] = len(tmp_bics)
        bics['data'] = []
        for bic1 in tmp_bics:
            c.execute("SELECT hm.name FROM bic_hal bh join hallmark hm on bh.hallmark_id=hm.id WHERE bh.bicluster_id=%s",
                      [bic1[3]])
            tmp1 = c.fetchall()
            h1 = list(set([convert[i[0]] for i in tmp1]))
            h2 = [[i[0],convert[i[0]]] for i in tmp1]
            bics['data'].append([bic1[4], bic1[5], bic1[6], bic1[7], bic1[8], h2])
    db.close()
    return render_template('search.html', gene=symbol, muts=muts, regs=regs, bics=bics)


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
    try:
        c = db.cursor()
        c.execute("""SELECT symbol FROM gene WHERE symbol LIKE %s""", [str(term)+'%'])
        json1 = json.dumps([i[0] for i in c.fetchall()])
    finally:
        db.close()
    return Response(response=json1, status=200, mimetype='application/json')


@app.route('/combinatorial_network')
def combinatorial_network():
    with open(app.config['NODES_FILE'], 'r') as infile:
        csvreader = csv.reader(infile, delimiter=',')
        csvreader.next()
        nodes = {node_id: {'id': node_id, 'tf_ko': tf_ko, 'in_gbm': in_gbm}
                 for node_id, tf_ko, in_gbm in csvreader}

    with open(app.config['EDGES_FILE'], 'r') as infile:
        csvreader = csv.reader(infile, delimiter=',')
        csvreader.next()
        edges = []
        for edge, sig_coocc in csvreader:
            source, edge_type, target = edge.split()
            edges.append({'source': source, 'target': target, 'type': edge_type,
                          'sig_coocc': sig_coocc})

    graph_data = []
    for node_id, node_data in nodes.items():
        classes = []
        if node_id.startswith('hsa-miR'):
            classes.append('mirna')
        else:
            classes.append('gene')

        if node_data['tf_ko'] == 'Yes':
            classes.append('crispr')
        if node_data['in_gbm'] == 'Yes':
            classes.append('in_gbm')
        if 'in_gbm' in classes and 'crispr' in classes:
            classes.append('crispr_gbm')
        graph_data.append({ 'data': { 'id': node_id }, 'classes': ' '.join(classes) })

    for i, edge in enumerate(edges):
        if edge['sig_coocc'] == 'Yes':
            graph_data.append({ 'data': { 'id': 'e%d' % i, 'source': edge['source'], 'target': edge['target'] }, 'classes': 'sigcoocc' })
        else:
            graph_data.append({ 'data': { 'id': 'e%d' % i, 'source': edge['source'], 'target': edge['target'] } })

    return render_template('combinatorial_network.html', **locals())
    

if __name__ == '__main__':
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    app.debug = True
    app.secret_key = 'supercalifragilistic'
    app.logger.addHandler(handler)
    app.run(host='0.0.0.0', debug=True)

