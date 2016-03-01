import MySQLdb
import re
db = MySQLdb.connect(host='localhost', user='glioma', passwd='glioma', db='glioma')

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        if a[-3:]=='-3p':
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(a+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(b):
            return 1
    else:
        if b[-3:]=='-3p':
            re1 = re.compile(b+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(b+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(a):
            return 1
    return 0

def miRNAInDict(miRNA, dict1):
    retMe = []
    for i in dict1.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe

### SHOW TABLES ###
c = db.cursor()
c.execute("""SHOW TABLES""")
print c.fetchall()

### GENE ###
# Load in gene information
genes = []
inFile = open('data/Homo_sapiens.gene_info','r')
inFile.readline() # Skip header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split('\t')
    genes.append([splitUp[1],splitUp[2]])
inFile.close()

# Enter in gene information into the database
#(('id', 'int(10) unsigned', 'NO', 'PRI', None, 'auto_increment'), ('symbol', 'varchar(100)', 'YES', 'MUL', None, ''), ('entrez', 'int(11)', 'YES', '', None, ''))
c.execute("""SHOW COLUMNS FROM gene""")
print c.fetchall()

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genesInDb = [i[1] for i in tmp]

#for gene in genes:
#    if not gene[1] in genesInDb:
#        c.execute("""INSERT INTO gene (symbol, entrez) VALUES (%s,%s)""", [gene[1],gene[0]])
#c.connection.commit()

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
print len(genes), len(tmp)


### PATIENTS ###
# Load in patient information
patients = []
inFile = open('data/conditions.txt','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    patients.append(line.strip())
inFile.close()

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM patient""")
print c.fetchall()

c.execute("""SELECT * FROM patient""")
tmp = c.fetchall()
patsInDb = [i[1] for i in tmp]

#for patient in patients:
#    if not patient in patsInDb:
#        c.execute("""INSERT INTO patient (name) VALUES (%s)""", [patient])
#c.connection.commit()

c.execute("""SELECT * FROM patient""")
tmp = c.fetchall()
print len(patients), len(tmp)


### HALLMARKS ###
hallmarks = ['Sustained angiogenesis', 'Insensitivity to antigrowth signals', 'Evading apoptosis', 'Limitless replicative potential', 'Evading immune detection', 'Tissue invasion and metastasis', 'Self sufficiency in growth signals', 'Tumor promoting inflammation', 'Reprogramming energy metabolism', 'Genome instability and mutation']

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM hallmark""")
print c.fetchall()

c.execute("""SELECT * FROM hallmark""")
tmp = c.fetchall()
hmInDb = [i[1] for i in tmp]

#for hallmark in hallmarks:
#    if not hallmark in hmInDb:
#        c.execute("""INSERT INTO hallmark (name) VALUES (%s)""", [hallmark])
#c.connection.commit()

c.execute("""SELECT * FROM hallmark""")
tmp = c.fetchall()
print len(hallmarks), len(tmp)


### GO_BP ###
geneOntology = []
inFile = open('data/gene_ontology_ext.obo','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    
    if line.strip()=='[Term]':
        goIds = []
        id1 = inFile.readline().strip().split('id: ')[1]
        name = inFile.readline().strip().split('name: ')[1]
        category = inFile.readline().strip().split('namespace: ')[1]
        if category=='biological_process':
            geneOntology.append([id1, name])
    
    if line.strip()[0:6]=='alt_id' and category=='biological_process':
        geneOntology.append([line.strip().split('alt_id: ')[1], name])

inFile.close()

geneOntology.append(['GO:0036017','response to erythropoietin'])
geneOntology.append(['GO:0036018','cellular response to erythropoietin'])
geneOntology.append(['GO:0038084','vascular endothelial growth factor signaling pathway'])
geneOntology.append(['GO:0044380','protein localization to cytoskeleton'])
geneOntology.append(['GO:1900046','regulation of hemostasis'])
geneOntology.append(['GO:2001141','regulation of RNA biosynthetic process'])
geneOntology.append(['GO:2001198','regulation of dendritic cell differentiation'])

# Enter in GO BP information into the database
c.execute("""SHOW COLUMNS FROM go_bp""")
print c.fetchall()

c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
gobpInDb = [i[1] for i in tmp]

#for go in geneOntology:
#    if not go[0] in gobpInDb:
#        c.execute("""INSERT INTO go_bp (go_id, name, description) VALUES (%s, %s, %s)""", [go[0],go[1],'NA'])
#c.connection.commit()

c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
print len(geneOntology), len(tmp)


### GO_BP TO GENE ###
c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
goTerms = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([str(i[2]) for i in tmp],[i[0] for i in tmp]))

go2gene = []
inFile = open('data/gene2go.hsa','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split('\t')
    if splitUp[1] in genes and splitUp[2] in goTerms:
        go2gene.append([goTerms[splitUp[2]], genes[splitUp[1]]])
inFile.close()

# Enter in GO BP information into the database
c.execute("""SHOW COLUMNS FROM go_gene""")
print c.fetchall()

c.execute("""SELECT * FROM go_gene""")
tmp = c.fetchall()
goGeneInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

#for goGene in go2gene:
#    if not str(goGene[0])+'_'+str(goGene[1]) in goGeneInDb:
#        c.execute("""INSERT INTO go_gene (go_bp_id, gene_id) VALUES (%s, %s)""", goGene)
#c.connection.commit()

c.execute("""SELECT * FROM go_gene""")
tmp = c.fetchall()
print len(go2gene), len(tmp)


### NCI Nature Pathways ###
# Load in NCI Nature Pathways information
nciNature = []
inFile = open('data/NCI_Nature_Pathways.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    nciNature.append(line.strip())
inFile.close()

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM nci_nature_pathway""")
print c.fetchall()

c.execute("""SELECT * FROM nci_nature_pathway""")
tmp = c.fetchall()
nciNatureInDb = [i[1] for i in tmp]

#for path in nciNature:
#    if not path in nciNatureInDb:
#        c.execute("""INSERT INTO nci_nature_pathway (name) VALUES (%s)""", [path])
#c.connection.commit()

c.execute("""SELECT * FROM nci_nature_pathway""")
tmp = c.fetchall()
print len(nciNature), len(tmp)


### miRNAs ###
# Load in miRNA information
miRNAs = []
miRNAIDs = {}
rev_miRNAIDs = {}
inFile = open('data/hsa.mature.fa','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    tmp = line.strip().split(' ')
    tmp[0] = re.sub('-5p','',tmp[0])
    miRNAs.append(tmp)
    miRNAIDs[tmp[0].lower()] = tmp[1]
    rev_miRNAIDs[tmp[1]] = tmp[0].lower()
inFile.close()

"""
# Load in miR2disease
miR2disease = {}
inFile = open('data/miR2Disease_5_2011.csv','r')
inFile.readline()
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    miRs = miRNAInDict(splitUp[0].lower(), miRNAIDs)
    if splitUp[0].lower()=='hsa-mir-17-92':
        miRs = [miRNAIDs['hsa-mir-17'], miRNAIDs['hsa-mir-18a'], miRNAIDs['hsa-mir-19a'], miRNAIDs['hsa-mir-20a'], miRNAIDs['hsa-mir-92a']]
    if len(miRs)==0:
        print 'ARRGGGHHH!!(miR2disease): ' + splitUp[0].lower()
    else:
        for miR in miRs:
            miR2disease[miR] = splitUp[1].lower()
inFile.close()

# Load in hmdd
hmdd = []
inFile = open('data/hmdd_9_9_2012.csv','r')
inFile.readline()
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    miRs = miRNAInDict(splitUp[0].lower(), miRNAIDs)
    if len(miRs)==0:
        print 'ARRGGGHHH!!(hmdd): ' + splitUp[0].lower()
    else:
        for miR in miRs:
            hmdd.append(miR)
inFile.close()
"""
# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM mirna""")
print c.fetchall()

c.execute("""SELECT * FROM mirna""")
tmp = c.fetchall()
miRNAInDb = [i[1] for i in tmp]

#for miRNA in miRNAs:
#    if not miRNA[1] in miRNAInDb:
#        m2d1 = 'no'
#        if miRNA[1] in miR2disease:
#            m2d1 = miR2disease[miRNA[1]]
#        hmdd1 = 0
#        if miRNA[1] in hmdd:
#            hmdd1 = 1
#        c.execute("""INSERT INTO mirna (mature_seq_id, name, mir2disease, hmdd) VALUES (%s, %s, %s, %s)""", [miRNA[1], miRNA[0].lower(), m2d1, hmdd1])
        #if not m2d1=='no' or hmdd1>0:
        #    print [miRNA[1], miRNA[0].lower(), m2d1, hmdd1]
#c.connection.commit()

c.execute("""SELECT * FROM mirna""")
tmp = c.fetchall()
print len(miRNAs), len(tmp)


### BICLUSTER ###
# Bicluster,Number_Of_Genes,Genes,Number_of_Tumors_in_Bicluster,Tumors_in_Bicluster,TCGA_Var_Exp_First_PC,TCGA_Var_Exp_First_PC_Perm_P_Value,TCGA_Survival,TCGA_Survival_P_Value,French_Var_Exp_First_PC,French_Average_Random_Var_Exp_First_PC,French_Var_Exp_First_PC_Perm_P_Value,French_Survival,French_Survival_P_Value,REMBRANDT_Var_Exp_First_PC,REMBRANDT_Average_Random_Var_Exp_First_PC,REMBRANDT_Var_Exp_First_PC_Perm_P_Value,REMBRANDT_Survival,REMBRANDT_Survival_P_Value,GSE7696_Var_Exp_First_PC,GSE7696_Average_Random_Var_Exp_First_PC,GSE7696_Var_Exp_First_PC_Perm_P_Value,GSE7696_Survival,GSE7696_Survival_P_Value,Enriched_GO_BP_Terms,Sustained_Angiogenesis,Insensitivity_To_Antigrowth_Signals,Evading_Apoptosis,Limitless_Replicative_Potential,Evading_Immune_Detection,Tissue_Invasion_And_Metastasis,Self_Sufficiency_In_Growth_Signals,Tumor_Promoting_Inflammation,Reprogramming_Energy_Metabolism,Genome_Instability_And_Mutation
# Load in bicluster information
biclusters = []
inFile = open('data/biclusters.csv','r')
header = inFile.readline().strip().split(',')
while 1:
    line = inFile.readline()
    if not line:
        break
    biclusters.append(dict(zip(header,line.strip().split(','))))
inFile.close()

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM bicluster""")
print c.fetchall()

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
biclustersInDb = [i[1] for i in tmp]

#for bic1 in biclusters:
#    if not bic1 in biclustersInDb:
#        c.execute("""INSERT INTO bicluster (name, var_exp_fpc, var_exp_fpc_p_value, survival, survival_p_value) VALUES (%s, %s, %s, %s, %s)""", [bic1['Bicluster'], bic1['TCGA_Var_Exp_First_PC'], bic1['TCGA_Var_Exp_First_PC_Perm_P_Value'], bic1['TCGA_Survival'], bic1['TCGA_Survival_P_Value']])
#c.connection.commit()

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
print len(biclusters), len(tmp)


### BICLUSTER GENES ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([str(i[1]).upper() for i in tmp],[i[0] for i in tmp]))

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM bic_gene""")
print c.fetchall()

c.execute("""SELECT * FROM bic_gene""")
tmp = c.fetchall()
bicGeneInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_genes = []
missing = []
for bic1 in biclusters:
    for gene in bic1['Genes'].split(' '):
        if gene in genes:
            bic_genes.append([bics[bic1['Bicluster']],genes[gene]])
        else:
            missing.append(gene)
print 'Missing genes = ',len(missing)

#outFile = open('missingGenes.csv','w')
#outFile.write('\n'.join(missing))
#outFile.close()

#for bg1 in bic_genes:
#    if not str(bg1[0])+'_'+str(bg1[1]) in bicGeneInDb:
#        c.execute("""INSERT INTO bic_gene (bicluster_id, gene_id) VALUES (%s, %s)""", bg1)
#c.connection.commit()

c.execute("""SELECT * FROM bic_gene""")
tmp = c.fetchall()
print len(bic_genes), len(tmp)


### BICLUSTER PATIENTS ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM patient""")
tmp = c.fetchall()
patients = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM bic_pat""")
print c.fetchall()

c.execute("""SELECT * FROM bic_pat""")
tmp = c.fetchall()
bicPatInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_pat = []
for bic1 in biclusters:
    for patient in bic1['Tumors_in_Bicluster'].split(' '):
        if patient in patients:
            bic_pat.append([bics[bic1['Bicluster']],patients[patient]])

#for bp1 in bic_pat:
#    if not str(bp1[0])+'_'+str(bp1[1]) in bicPatInDb:
#        c.execute("""INSERT INTO bic_pat (bicluster_id, patient_id) VALUES (%s, %s)""", bp1)
#c.connection.commit()

c.execute("""SELECT * FROM bic_pat""")
tmp = c.fetchall()
print len(bic_pat), len(tmp)


### BICLSUTER REPLICATION ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM replication""")
print c.fetchall()

c.execute("""SELECT * FROM replication""")
tmp = c.fetchall()
repInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
replications = []
for bic1 in biclusters:
    for study in ['French','REMBRANDT','GSE7696']:
        replications.append([bics[bic1['Bicluster']], study, bic1[study+'_'+'Var_Exp_First_PC'], bic1[study+'_'+'Var_Exp_First_PC_Perm_P_Value'], bic1[study+'_'+'Survival'], bic1[study+'_'+'Survival_P_Value']])

#for rep1 in replications:
#    if not str(rep1[0])+'_'+str(rep1[1]) in repInDb:
#        c.execute("""INSERT INTO replication (bicluster_id, study, var_exp_fpc, var_exp_fpc_p_value, survival, survival_p_value) VALUES (%s, %s, %s, %s, %s, %s)""", rep1)
#c.connection.commit()

c.execute("""SELECT * FROM replication""")
tmp = c.fetchall()
print len(replications), len(tmp)


### BICLUSTER GO_BP ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM go_bp""")
tmp = c.fetchall()
terms = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM bic_go""")
print c.fetchall()

c.execute("""SELECT * FROM bic_go""")
tmp = c.fetchall()
bicGoInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_go = []
#missing = []
for bic1 in biclusters:
    if not bic1['Enriched_GO_BP_Terms']=='NA':
        for go in bic1['Enriched_GO_BP_Terms'].split(';'):
            if go in terms:
                bic_go.append([bics[bic1['Bicluster']], terms[go]])
            else:
                print bic1['Bicluster'], go
                #missing.append(go)

#outFile = open('missingTerms.csv','w')
#outFile.write('\n'.join(missing))
#outFile.close()

#for bg1 in bic_go:
#    if not str(bg1[0])+'_'+str(bg1[1]) in bicGoInDb:
#        c.execute("""INSERT INTO bic_go (bicluster_id, go_bp_id) VALUES (%s, %s)""", bg1)
#c.connection.commit()

c.execute("""SELECT * FROM bic_go""")
tmp = c.fetchall()
print len(bic_go), len(tmp)


### BICLUSTER HALLMARKS ###
c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

#c.execute("""SELECT * FROM hallmark""")
#tmp = c.fetchall()
#hallmarks = dict(zip([str(i[1]) for i in tmp],[i[0] for i in tmp]))

hallmarks = {'Self_Sufficiency_In_Growth_Signals': 7L, 'Genome_Instability_And_Mutation': 10L, 'Insensitivity_To_Antigrowth_Signals': 2L, 'Sustained_Angiogenesis': 1L, 'Tumor_Promoting_Inflammation': 8L, 'Limitless_Replicative_Potential': 4L, 'Reprogramming_Energy_Metabolism': 9L, 'Evading_Apoptosis': 3L, 'Evading_Immune_Detection': 5L, 'Tissue_Invasion_And_Metastasis': 6L}

# Enter in patient information into the database
c.execute("""SHOW COLUMNS FROM bic_hal""")
print c.fetchall()

c.execute("""SELECT * FROM bic_hal""")
tmp = c.fetchall()
bicHalInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

# Create list of bicluster to gene mappings
bic_hal = []
for bic1 in biclusters:
    for hal1 in hallmarks:
        if not bic1[hal1]=='NA' and float(bic1[hal1])>=0.8:
            bic_hal.append([bics[bic1['Bicluster']], hallmarks[hal1]])

#for bh1 in bic_hal:
#    if not str(bh1[0])+'_'+str(bh1[1]) in bicHalInDb:
#        c.execute("""INSERT INTO bic_hal (bicluster_id, hallmark_id) VALUES (%s, %s)""", bh1)
#c.connection.commit()

c.execute("""SELECT * FROM bic_hal""")
tmp = c.fetchall()
print len(bic_hal), len(tmp)


### TF REGULATORS ###
c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

tfRegulators = []
inFile = open('data/tfRegulators_V2.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    tfs = splitUp[1].split(' ')
    for tf in tfs:
        if tf in genes:
            tfRegulators.append([bics[splitUp[0]], genes[tf], splitUp[2]])
        else:
            print tf, splitUp[0]

inFile.close()

# Enter in TF reglator information into the database
c.execute("""SHOW COLUMNS FROM tf_regulator""")
print c.fetchall()

c.execute("""SELECT * FROM tf_regulator""")
tmp = c.fetchall()
tfRegInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

#for tf in tfRegulators:
#    if not str(tf[0])+'_'+str(tf[1]) in tfRegInDb:
#        c.execute("""INSERT INTO tf_regulator (bicluster_id, gene_id, action) VALUES (%s, %s, %s)""", tf)
#c.connection.commit()

c.execute("""SELECT * FROM tf_regulator""")
tmp = c.fetchall()
print len(tfRegulators), len(tmp)


### miRNA REGULATORS ###
c.execute("""SELECT * FROM mirna""")
tmp = c.fetchall()
mirnas = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

mirnaRegulators = []
"""inFile = open('data/miRNARegulators.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    miRs = splitUp[1].lower().split(' ')
    for miR in miRs:
        miRs2 = miRNAInDict(miR, miRNAIDs)
        if not len(miRs2)==0:
            for miR2 in miRs2:
                mirnaRegulators.append([bics[splitUp[0]], mirnas[miR2]])
        else:
            print miR
inFile.close()
"""

# Enter in miRNA reglator information into the database
c.execute("""SHOW COLUMNS FROM mirna_regulator""")
print c.fetchall()

c.execute("""SELECT * FROM mirna_regulator""")
tmp = c.fetchall()
mirnaRegInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

#for mirna in mirnaRegulators:
#    if not str(miRNA[0])+'_'+str(miRNA[1]) in mirnaRegInDb:
#        c.execute("""INSERT INTO mirna_regulator (bicluster_id, mirna_id) VALUES (%s, %s)""", mirna)
#c.connection.commit()

c.execute("""SELECT * FROM mirna_regulator""")
tmp = c.fetchall()
print len(mirnaRegulators), len(tmp)


### CAUSAL FLOWS ###
c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM mirna""")
tmp = c.fetchall()
mirnas = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM nci_nature_pathway""")
tmp = c.fetchall()
nci_nat = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

c.execute("""SELECT * FROM bicluster""")
tmp = c.fetchall()
bics = dict(zip([i[1] for i in tmp],[i[0] for i in tmp]))

somaticMutations = []
causalFlows = []
inFile = open('data/causality.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    somMut = ''
    if splitUp[1]=='gene':
        somaticMutations.append([genes[splitUp[0]],splitUp[1]])
        somMut = str(genes[splitUp[0]])+'_'+splitUp[1]
    elif splitUp[1]=='pathway':
        somaticMutations.append([nci_nat[splitUp[0]],splitUp[1]])
        somMut = str(nci_nat[splitUp[0]])+'_'+splitUp[1]
    if splitUp[3]=='TF':
        causalFlows.append([somMut, genes[splitUp[2]], 'tf', bics[splitUp[4]], splitUp[5], splitUp[6]])
    if splitUp[3]=='miRNA':
        miR = miRNAInDict(splitUp[2].lower(), miRNAIDs)[0]
        if not miR=='':
            causalFlows.append([somMut, mirnas[miR], 'mirna', bics[splitUp[4]], splitUp[5], splitUp[6]])
inFile.close()

# Somatic mutations (genes or pathways)
c.execute("""SHOW COLUMNS FROM somatic_mutation""")
print c.fetchall()

c.execute("""SELECT * FROM somatic_mutation""")
tmp = c.fetchall()
somMutInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

somaticMutations2 = [i.split('_') for i in list(set([str(i[0])+'_'+str(i[1]) for i in somaticMutations]))]
#for somMut in somaticMutations2:
#    if not str(somMut[0])+'_'+str(somMut[1]) in somMutInDb:
#        c.execute("""INSERT INTO somatic_mutation (ext_id, mutation_type) VALUES (%s, %s)""", somMut)
#c.connection.commit()

c.execute("""SELECT * FROM somatic_mutation""")
tmp = c.fetchall()
print len(somaticMutations2), len(tmp)
somMuts = dict(zip([str(i[1])+'_'+str(i[2]) for i in tmp],[i[0] for i in tmp]))

# Causal flows
c.execute("""SHOW COLUMNS FROM causal_flow""")
print c.fetchall()

c.execute("""SELECT * FROM causal_flow""")
tmp = c.fetchall()
causalFlowsInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

#for causalFlow in causalFlows:
#    if not str(causalFlow[0])+'_'+str(causalFlow[1]) in causalFlowsInDb:
#        sm1 = somMuts[causalFlow.pop(0)]
#        c.execute("""INSERT INTO causal_flow (somatic_mutation_id, regulator_id, regulator_type, bicluster_id, leo_nb_atob, mlogp_m_atob) VALUES (%s, %s, %s, %s, %s, %s)""", [sm1]+causalFlow)
#c.connection.commit()

c.execute("""SELECT * FROM causal_flow""")
tmp = c.fetchall()
print len(causalFlows), len(tmp)


### TF CRISPR-Cas9 ###
c.execute("""SELECT * FROM gene""")
tmp = c.fetchall()
genes = dict(zip([str(i[2]) for i in tmp],[i[0] for i in tmp]))

# TF,1502.logFC,1502.FDR,U5.logFC,U5.FDR,131.logFC,131.FDR,CB660.logFC,CB660.FDR,827.logFC,827.FDR
tfCrispr = []
inFile = open('data/sgRNA_data.csv','r')
header = inFile.readline().strip().split(',')
cellLines = ['131','827','CB660','U5']
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = dict(zip(header,line.strip().split(',')))
    # Add each cell line data into tfCrispr
    # TODO
    for cellLine in cellLines:
        tfCrispr.append([genes[splitUp['TF']], cellLine, splitUp[cellLine+'.logFC'], splitUp[cellLine+'.FDR']])

inFile.close()

# Enter in TF reglator information into the database
c.execute("""SHOW COLUMNS FROM tf_crispr""")
print c.fetchall()

c.execute("""SELECT * FROM tf_crispr""")
tmp = c.fetchall()
tcInDb = [str(i[1])+'_'+str(i[2]) for i in tmp]

#for tc in tfCrispr:
#    if not str(tc[0])+'_'+str(tf[1]) in tcInDb:
#        c.execute("""INSERT INTO tf_crispr (gene_id, cell_line, log2fc, fdr) VALUES (%s, %s, %s, %s)""", tc)
#c.connection.commit()

c.execute("""SELECT * FROM tf_crispr""")
tmp = c.fetchall()
print len(tfCrispr), len(tmp)

