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

convert = {'Evading apoptosis':'cellDeath.gif', 'Evading immune detection':'avoidImmuneDestruction.gif', 'Genome instability and mutation':'genomicInstability.gif', 'Insensitivity to antigrowth signals':'evadeGrowthSuppressors.gif', 'Limitless replicative potential':'immortality.gif', 'Reprogramming energy metabolism':'cellularEnergetics.gif', 'Self sufficiency in growth signals':'sustainedProliferativeSignalling.gif', 'Sustained angiogenesis':'angiogenesis.gif', 'Tissue invasion and metastasis':'invasion.gif', 'Tumor promoting inflammation':'promotingInflammation.gif'}

app = Flask(__name__)
app.config.from_envvar('GLIOMA_SETTINGS')

######################################################################
#### General helpers
######################################################################

def dbconn():
    return MySQLdb.connect(host=app.config['HOST'], user=app.config['USER'],
                           passwd=app.config['PASS'], db=app.config['DB'])


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

@app.route('/bicluster')
def bicluster():
    bicluster = request.args.get('bicluster')
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT * FROM bicluster WHERE name=%s""", [bicluster])
    bicInfo = c.fetchall()[0]
    tmp = [0,0]
    if float(bicInfo[3])<=0.05:
        tmp[0] = 1
    if float(bicInfo[5])<=0.05:
        tmp[1] = 1
    bicInfo = list(bicInfo) + tmp
    c.execute("""SELECT gene.id, gene.symbol, gene.entrez FROM bic_gene, gene WHERE bic_gene.bicluster_id=%s AND gene.id=bic_gene.gene_id""", [bicInfo[0]])
    genes = sorted(c.fetchall(), key=lambda symbol: symbol[1])
    c.execute("""SELECT patient.id, patient.name FROM bic_pat, patient WHERE bic_pat.bicluster_id=%s AND patient.id=bic_pat.patient_id""", [bicInfo[0]])
    tumors = sorted(c.fetchall(), key=lambda name: name[1])
    # Replication
    c.execute("""SELECT * FROM replication WHERE bicluster_id=%s""", [bicInfo[0]])
    tmp = list(c.fetchall())
    repConvert = {'French':'Gravendeel, et al. 2009','REMBRANDT':'Madhavan, et al. 2009','GSE7696':'Murat, et al. 2008'}
    repPubmed = {'French':'19920198','REMBRANDT':'19208739','GSE7696':'18565887'}
    replication = []
    replicated = [0, 0]
    for i in tmp:
        tmp1 = [0,0]
        if float(bicInfo[3])<=0.05 and float(i[4])<=0.05:
            tmp1[0] = 1
            replicated[0] = 1
        if (float(bicInfo[4])>0 and (float(i[5])>0 and float(i[6])<=0.05)) or (float(bicInfo[4])<0 and (float(i[5])<0 and float(i[6])<=0.05)):
            tmp1[1] = 1
            replicated[1] = 1
        replication.append(list(i)+[repConvert[i[2]], repPubmed[i[2]]]+tmp1)
    bicInfo = bicInfo + replicated
    # Regulators
    regulators = []
    c.execute("""SELECT gene.id, gene.symbol, tf_regulator.action FROM tf_regulator, gene WHERE tf_regulator.bicluster_id=%s AND gene.id=tf_regulator.gene_id""", [bicInfo[0]])
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
    c.execute("""SELECT mirna.id, mirna.name, mirna.mir2disease, mirna.hmdd FROM mirna_regulator, mirna WHERE mirna_regulator.bicluster_id=%s AND mirna.id=mirna_regulator.mirna_id""", [bicInfo[0]])
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
    c.execute("""SELECT * FROM causal_flow WHERE bicluster_id=%s""", [bicInfo[0]])
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
    c.execute("""SELECT hallmark.name FROM hallmark, bic_hal WHERE bic_hal.bicluster_id=%s AND hallmark.id=bic_hal.hallmark_id""", [bicInfo[0]])
    h1 = list(c.fetchall())
    h2 = [[i[0],convert[i[0]]] for i in h1]
    # GO
    c.execute("""SELECT go_bp.id, go_bp.go_id, go_bp.name FROM bic_go, go_bp WHERE bic_go.bicluster_id=%s AND go_bp.id=bic_go.go_bp_id""", [bicInfo[0]])
    tmps = list(c.fetchall())
    gobps = []
    for gobp in tmps:
        c.execute("""SELECT gene.symbol FROM go_gene, gene, bic_gene WHERE go_gene.go_bp_id=%s AND bic_gene.bicluster_id=%s AND go_gene.gene_id=gene.id AND gene.id=bic_gene.gene_id""", [gobp[0], bicInfo[0]])
        tmp = c.fetchall()
        gobps.append(list(gobp)+[sorted(list(set([i[0] for i in tmp])))])
    
    return render_template('bicluster.html', bicluster=bicluster, genes=genes, tumors=tumors, bicInfo=bicInfo, replication=replication, regulators=regulators, hallmarks=h2, gobps=gobps, causalFlows=causalFlows)

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

