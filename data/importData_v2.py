import re, cPickle, os
from glioma.models import *
from copy import deepcopy
from cMonkeyWrapper import cMonkeyWrapper
from bicluster import bicluster
from pssm import pssm

#######################################################################
## Create a dictionary to convert the miRNAs to there respective ids ##
#######################################################################
inFile = open('data/hsa.mature.fa','r')
miRNAIDs = {}
miRNAIDs_rev = {}
mature_sequence_ids = []
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.split(' ')
    mature_sequence_ids.append(splitUp[1])
    if not splitUp[1] in miRNAIDs_rev:
        miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()
    if not splitUp[0].lower() in miRNAIDs:
        miRNAIDs[splitUp[0].lower()] = splitUp[1]
    else:
        print 'Uh oh!',splitUp

def miRNAInDict(miRNA, dict1):
    retMe = []
    for i in dict1.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        re1 = re.compile(a+'[a-z]$')
        if re1.match(b):
            return 1
    else:
        re1 = re.compile(b+'[a-z]$')
        if re1.match(a):
            return 1
    return 0

# Load up the miRNAs
miRNAs = {}
m_all = MiRNA.objects.all()
for m1 in m_all:
    miRNAs[m1.mature_sequence_id] = deepcopy(m1)

# Load up the biclusters
biclusters = {}
b_all = Bicluster.objects.all()
for b1 in b_all:
    biclusters[b1.name] = deepcopy(b1)

# Load up the datasets
datasets = {}
datasets['Rembrandt'] = Dataset.objects.get(name='rembrandt')
datasets['French'] = Dataset.objects.get(name='french')
datasets['GSE7696'] = Dataset.objects.get(name='gse7696')

# Load up categories
categories = {}
c_all = Category.objects.all()
for c1 in c_all:
    categories[c1.name] = deepcopy(c1)

# Load up categories
genes = {}
g_all = Gene.objects.all()
for g1 in g_all:
    genes[g1.probe] = deepcopy(g1)

# Load up categories
patients = {}
p_all = Patient.objects.all()
for p1 in p_all:
    patients[p1.name] = deepcopy(p1)

# Load up the patient categories
patient_categories = {}
inFile = open('data/rug.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    patient_categories[splitUp[0]] = categories[splitUp[1]]

inFile.close()

#############################################
## Load up the bicluster data for all runs ##
#############################################
# Load up p3utr results - read in pickled object of cMonkey results
postProcessed = {}
#genes = {}
#patients = {}
biclusters = {}
runs = ['all','upstream_p3utr']
run = 'upstream_p3utr'
pklFile = open('data/combined/'+str(run)+'/c1.pkl','rb')
c1_p3utr = cPickle.load(pklFile)
pklFile.close()
# Read in results from postProcessedVFinal.csv
postProcessed[run] = {}
# id,k.rows,k.cols,resid,resid.norm,motif1.E,motif1.consensus,motif1.matches,motif1.permutedEV<=10,motif1.permPV,motif2.E,motif2.consensus,motif2.matches,motif2.permutedEV<=10,motif2.permPV,3pUTRmotif1.weederScore,3pUTRmotif1.localPermP,3pUTRmotif1.allPermP,3pUTRmotif1.consensus,3pUTRmotif1.miRNAs,3pUTRmotif1.model,3pUTR_pita.miRNAs,3pUTR_pita.percTargets,3pUTR_pita.pValue,3pUTR_targetScan.miRNAs,3pUTR_targetScan.percTargets,3pUTR_targetScan.pValue,SEX.bi,SEX.bi.p,AGE,AGE.p,Survival,Survival.p,Survival.AGE,Survival.AGE.p,Survival.var,Survival.var.p,Survival.var.AGE,Survival.var.AGE.p,French_new.resid.norm,French_avg.resid.norm,French_norm.perm.p,French_survival,French_survival.p,French_survival.age,French_survival.age.p,GSE7696_new.resid.norm,GSE7696_avg.resid.norm,GSE7696_norm.perm.p,GSE7696_survival,GSE7696_survival.p,GSE7696_survival.age,GSE7696_survival.age.p,NON_TUMOR,ASTROCYTOMA,MIXED,OLIGODENDROGLIOMA,GBM
#headers = ['k.rows','k.cols','resid','resid.norm','motif1.E','motif1.consensus','motif1.matches','motif1.permutedEV<=10','motif1.permPV','motif2.E','motif2.consensus','motif2.matches','motif2.permutedEV<=10','motif2.permPV','3pUTRmotif1.weederScore','3pUTRmotif1.localPermP','3pUTRmotif1.allPermP','3pUTRmotif1.consensus','3pUTRmotif1.miRNAs','3pUTRmotif1.model','3pUTR_pita.miRNAs','3pUTR_pita.percTargets','3pUTR_pita.pValue','3pUTR_targetScan.miRNAs','3pUTR_targetScan.percTargets','3pUTR_targetScan.pValue','SEX.bi','SEX.bi.p','AGE','AGE.p','Survival','Survival.p','Survival.AGE','Survival.AGE.p','Survival.var','Survival.var.p','Survival.var.AGE','Survival.var.AGE.p','French_new.resid.norm','French_avg.resid.norm','French_norm.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p','NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']
headers = ['k.rows','k.cols','resid','resid.norm','motif1.E','motif1.consensus','motif1.matches','motif1.permutedEV<=10','motif1.permPV','motif2.E','motif2.consensus','motif2.matches','motif2.permutedEV<=10','motif2.permPV','3pUTRmotif1.weederScore','3pUTRmotif1.localPermP','3pUTRmotif1.allPermP','3pUTRmotif1.consensus','3pUTRmotif1.miRNAs','3pUTRmotif1.model','SEX.bi','SEX.bi.p','AGE','AGE.p','Survival','Survival.p','Survival.AGE','Survival.AGE.p','Survival.var','Survival.var.p','Survival.var.AGE','Survival.var.AGE.p','French_new.resid.norm','French_avg.resid.norm','French_norm.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p','NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']
inFile = open('data/combined/'+run+'/postProcessedVFinal.csv','r')
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    bicluster = int(splitUp.pop(0))
    postProcessed[run][bicluster] = dict(zip(headers,splitUp))

inFile.close()
# Get the biclusters from the cMonkey object
all_biclusters = c1_p3utr.getBiclusters()
for bicluster in all_biclusters:
    # Initialize the bicluster
    inThere = 0
    print bicluster
    b1 = Bicluster(name=str(run)+'.'+str(bicluster), residual=float(postProcessed[run][bicluster]['resid']), residual_normalized=float(postProcessed[run][bicluster]['resid.norm']), residual_normalized_permuted_p=float(1), survival_beta=float(postProcessed[run][bicluster]['Survival']), survival_p=float(postProcessed[run][bicluster]['Survival.p']), survival_age_beta=float(postProcessed[run][bicluster]['Survival.AGE']), survival_age_p=float(postProcessed[run][bicluster]['Survival.AGE.p']) )
    b1.save()
    # Put genes into database if not already there
    print 'Adding genes...'
    genes_in_bicluster = all_biclusters[bicluster].getGenes()
    for gene in genes_in_bicluster:
        if not gene in genes:
            genes[gene] = Gene(probe=gene)
            print 'Adding gene!'
            genes[gene].save()
        b1.genes.add(genes[gene])
    # Put patients into database if not already there
    print 'Adding patients...'
    patients_in_bicluster = all_biclusters[bicluster].getConditions()
    for patient in patients_in_bicluster:
        if not patient in patients:
            patients[patient] = Patient(name=patient, category = patient_categories[patient])
            print 'Adding patient!'
            patients[patient].save()
            datasets['Rembrandt'].patients.add(patients[patient])
        b1.patients.add(patients[patient])
    
    biclusters[str(run)+'.'+str(bicluster)] = b1
    
    print 'Adding replication...'
    # Load up the replication data for all runs
    repsets = ['French','GSE7696']
    for repset in repsets:
        if not postProcessed[run][bicluster][repset+'_norm.perm.p']=='NA':
            print repset,'Replication!!'
            r1 = Replication(dataset = datasets[repset], bicluster = b1, residual = 1, residual_normalized = float(postProcessed[run][bicluster][repset+'_new.resid.norm']), residual_normalized_permuted_p = float(postProcessed[run][bicluster][repset+'_norm.perm.p']), survival_beta = float(postProcessed[run][bicluster][repset+'_survival']), survival_p = float(postProcessed[run][bicluster][repset+'_survival.p']), survival_age_beta = float(postProcessed[run][bicluster][repset+'_survival.age']), survival_age_p = float(postProcessed[run][bicluster][repset+'_survival.age.p']) )
            r1.save()
    
    print 'Adding enrichment...'
    # Load enrichment in or out of the bicluster
    enrichment_categories = ['NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']
    for cat in enrichment_categories:
        if not postProcessed[run][bicluster][cat]=='-Inf':
            print cat, 'Enrichment!'
            e1 = Enrichment(bicluster = b1, category = categories[cat], negative_log_p = postProcessed[run][bicluster][cat])
            e1.save()
    
    # Load regulatory factor information
    # 'motif1.E','motif1.consensus','motif1.matches','motif1.permutedEV<=10','motif1.permPV'
    print 'Adding upstream factor...'
    if (not postProcessed[run][bicluster]['motif1.E']=='NA') and (not postProcessed[run][bicluster]['motif1.permPV']=='NA'):
        # Load upstream motif 1
        u1 = Upstream(bicluster=b1, consensus = postProcessed[run][bicluster]['motif1.consensus'], e_value = float(postProcessed[run][bicluster]['motif1.E']), motif_matches = postProcessed[run][bicluster]['motif1.matches'], permuted_p = float(postProcessed[run][bicluster]['motif1.permPV']) )
        print 'Upstream Motif 1!'
        u1.save()
    
    if (not postProcessed[run][bicluster]['motif2.E']=='NA') and (not postProcessed[run][bicluster]['motif2.permPV']=='NA'):
        # Load upstream motif 2
        u1 = Upstream(bicluster=b1, consensus = postProcessed[run][bicluster]['motif2.consensus'], e_value = float(postProcessed[run][bicluster]['motif2.E']), motif_matches = postProcessed[run][bicluster]['motif2.matches'], permuted_p = float(postProcessed[run][bicluster]['motif2.permPV']) )
        print 'Upstream Motif 2!'
        u1.save()
    
    # Load p3utr motif
    # '3pUTRmotif1.weederScore','3pUTRmotif1.localPermP','3pUTRmotif1.allPermP','3pUTRmotif1.consensus','3pUTRmotif1.miRNAs','3pUTRmotif1.model'
    print 'Adding p3utr factor...'
    # Load p3utr motifs
    if not postProcessed[run][bicluster]['3pUTRmotif1.weederScore']=='NA':
        # Initialize object
        p1 = P3utr(bicluster = b1, consensus = postProcessed[run][bicluster]['3pUTRmotif1.consensus'], model = postProcessed[run][bicluster]['3pUTRmotif1.model'] )
        print 'P3UTR Motif!'
        p1.save()
        # Add miRNA motif matches
        if not postProcessed[run][bicluster]['3pUTRmotif1.miRNAs']=='NA':
            splitUp = postProcessed[run][bicluster]['3pUTRmotif1.miRNAs'].strip().split('_')
            for mirna in splitUp:
                if miRNA[-3:]=='-5p':
                    miRNA = miRNA[:-3]
                    miRs = miRNAInDict(miRNA.lower(), miRNAIDs)
                    if len(miRs)>0:
                        for miRNA in miRs:
                            p1.mirnas.add(miRNAs[miRNA])
    
    # Load pita predictions
    print 'Adding pita factor...'
    if not postProcessed[run][bicluster]['3pUTR_pita.miRNAs']=='NA':
        # Initialize object
        pt1 = TGPD(bicluster=b1, database = 'pita', percent_targets = float(postProcessed[run][bicluster]['3pUTR_pita.percTargets']), p = float(postProcessed[run][bicluster]['3pUTR_pita.pValue']) )
        print 'PITA!'
        pt1.save()
        # Add miRNA motif matches
        splitUp = postProcessed[run][bicluster]['3pUTR_pita.miRNAs'].strip().split(';')
        for miRNA in splitUp:
            if miRNA[-3:]=='-5p':
                miRNA = miRNA[:-3]
            miRs = miRNAInDict(miRNA.lower(), miRNAIDs)
            if len(miRs)>0:
                print miRs
                for miRNA in miRs:
                    pt1.mirnas.add(miRNAs[miRNA])
            else:
                print miRNA.lower(), miRs
    
    # Load targetscan predictions
    # '3pUTR_targetScan.miRNAs','3pUTR_targetScan.percTargets','3pUTR_targetScan.pValue'
    print 'Adding targetscan factor...'
    if not postProcessed[run][bicluster]['3pUTR_targetScan.miRNAs']=='NA':
        # Initialize object
        t1 = TGPD(bicluster=b1, database = 'targetscan', percent_targets = float(postProcessed[run][bicluster]['3pUTR_targetScan.percTargets']), p = float(postProcessed[run][bicluster]['3pUTR_targetScan.pValue']) )
        print 'TARGETSCAN!'
        t1.save()
        # Add miRNA motif matches
        splitUp = postProcessed[run][bicluster]['3pUTR_targetScan.miRNAs'].strip().split(';')
        for miRNA in splitUp:
            if miRNA[-3:]=='-5p':
                miRNA = miRNA[:-3]
            miRs = miRNAInDict(miRNA.lower(), miRNAIDs)
            if len(miRs)>0:
                print miRs
                for miRNA in miRs:
                    t1.mirnas.add(miRNAs[miRNA])
                
            else:
                print miRNA.lower(), miRs



Replication.objects.filter(bicluster=b1).delete()
Enrichment.objects.filter(bicluster=b1).delete()
b1.delete()

#######################################
## Load GO BP functional enrichments ##
#######################################
# Load up gene ontology
geneOntology = {}
go_all = Gene_Ontology.objects.all()
for go1 in go_all:
    geneOntology[go1.go_id] = deepcopy(go1)

runs = ['all','upstream_p3utr']
for run in runs:
    print run
    inFile = open('data/combined/'+run+'/biclusterEnrichment_GOBP.csv','r')
    inFile.readline() # Get rid of header
    # "","Top10.Terms.BP","BH.sig.GO.Ids.BP"
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        if not splitUp[2]=='':
            go_ids = splitUp[2].split('|')
            for go_id in go_ids:
                fe1 = Functional_Enrichment(bicluster = biclusters[run+'.'+splitUp[0]], gene_ontology = geneOntology[go_id])
                fe1.save()

inFile.close()


#####################################################
## Load up Hallmarks of Cancer Semantic Similarity ##
#####################################################
hallmarks = ['Self Sufficiency in Growth Signals','Insensitivity to Antigrowth Signals','Evading Apoptosis','Limitless Replicative Potential','Sustained Angiogenesis','Tissue Invasion and Metastasis','Genome Instability and Mutation','Tumor Promoting Inflammation','Reprogramming Energy Metabolism','Evading Immune Detection']
for run in runs:
    print run
    inFile = open('data/combined/'+run+'/jiangConrath_hallmarks.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        if not splitUp[1]=='NA':
            cluster_id= splitUp[0].strip('"')
            for i in range(len(hallmarks)):
                hoc1 = Hallmarks_Of_Cancer(bicluster=biclusters[run+'.'+cluster_id],hallmark_of_cancer=hallmarks[i],semantic_similarity=float(splitUp[i+1]))
                hoc1.save()
    
    inFile.close()

# Now read in permuted p-value files for each run
for run in runs:
    print run
    inFile = open('data/combined/'+run+'/residualPermutedPvalues_permAll.csv','r')
    inFile.readline() # Get rid of header
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        if not splitUp[1]=='NA':
            biclusters[run+'.'+splitUp[0].strip('"')].residual_normalized_permuted_p = float(splitUp[10])
            biclusters[run+'.'+splitUp[0].strip('"')].save()
    
    inFile.close()


#############################################
## Load up the bicluster data for all runs ##
#############################################
import re, cPickle, os
from glioma.models import *
from copy import deepcopy
from cMonkeyWrapper import cMonkeyWrapper
from bicluster import bicluster
from pssm import pssm

#######################################################################
## Create a dictionary to convert the miRNAs to there respective ids ##
#######################################################################
inFile = open('data/hsa.mature.fa','r')
miRNAIDs = {}
miRNAIDs_rev = {}
mature_sequence_ids = []
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.split(' ')
    mature_sequence_ids.append(splitUp[1])
    if not splitUp[1] in miRNAIDs_rev:
        miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()
    if not splitUp[0].lower() in miRNAIDs:
        miRNAIDs[splitUp[0].lower()] = splitUp[1]
    else:
        print 'Uh oh!',splitUp

def miRNAInDict(miRNA, dict1):
    retMe = []
    for i in dict1.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        re1 = re.compile(a+'[a-z]$')
        if re1.match(b):
            return 1
    else:
        re1 = re.compile(b+'[a-z]$')
        if re1.match(a):
            return 1
    return 0

# Load up the miRNAs
miRNAs = {}
m_all = MiRNA.objects.all()
for m1 in m_all:
    miRNAs[m1.mature_sequence_id] = deepcopy(m1)

# Load up the biclusters
biclusters = {}
b_all = Bicluster.objects.all()
for b1 in b_all:
    biclusters[b1.name] = deepcopy(b1)

# Load up the datasets
datasets = {}
datasets['Rembrandt'] = Dataset.objects.get(name='rembrandt')
datasets['French'] = Dataset.objects.get(name='french')
datasets['GSE7696'] = Dataset.objects.get(name='gse7696')

# Load up categories
categories = {}
c_all = Category.objects.all()
for c1 in c_all:
    categories[c1.name] = deepcopy(c1)

# Load up categories
genes = {}
g_all = Gene.objects.all()
for g1 in g_all:
    genes[g1.probe] = deepcopy(g1)

# Load up categories
patients = {}
p_all = Patient.objects.all()
for p1 in p_all:
    patients[p1.name] = deepcopy(p1)

# Load up the patient categories
patient_categories = {}
inFile = open('data/rug.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    patient_categories[splitUp[0]] = categories[splitUp[1]]

inFile.close()

# Load up p3utr results - read in pickled object of cMonkey results
postProcessed = {}
#genes = {}
#patients = {}
biclusters = {}
runs = ['all','upstream_p3utr']
#run = 'upstream_p3utr'
run = 'all'
pklFile = open('data/combined/'+str(run)+'/c1.pkl','rb')
c1_p3utr = cPickle.load(pklFile)
pklFile.close()
# Read in results from postProcessedVFinal.csv
postProcessed[run] = {}
headers = ['k.rows','k.cols','resid','resid.norm','motif1.E','motif1.consensus','motif1.matches','motif1.permutedEV<=10','motif1.permPV','motif2.E','motif2.consensus','motif2.matches','motif2.permutedEV<=10','motif2.permPV','3pUTRmotif1.weederScore','3pUTRmotif1.localPermP','3pUTRmotif1.allPermP','3pUTRmotif1.consensus','3pUTRmotif1.miRNAs','3pUTRmotif1.model','3pUTR_pita.miRNAs','3pUTR_pita.percTargets','3pUTR_pita.pValue','3pUTR_targetScan.miRNAs','3pUTR_targetScan.percTargets','3pUTR_targetScan.pValue','SEX.bi','SEX.bi.p','AGE','AGE.p','Survival','Survival.p','Survival.AGE','Survival.AGE.p','Survival.var','Survival.var.p','Survival.var.AGE','Survival.var.AGE.p','French_new.resid.norm','French_avg.resid.norm','French_norm.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p','NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']
#headers = ['k.rows','k.cols','resid','resid.norm','motif1.E','motif1.consensus','motif1.matches','motif1.permutedEV<=10','motif1.permPV','motif2.E','motif2.consensus','motif2.matches','motif2.permutedEV<=10','motif2.permPV','3pUTRmotif1.weederScore','3pUTRmotif1.localPermP','3pUTRmotif1.allPermP','3pUTRmotif1.consensus','3pUTRmotif1.miRNAs','3pUTRmotif1.model','SEX.bi','SEX.bi.p','AGE','AGE.p','Survival','Survival.p','Survival.AGE','Survival.AGE.p','Survival.var','Survival.var.p','Survival.var.AGE','Survival.var.AGE.p','French_new.resid.norm','French_avg.resid.norm','French_norm.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p','NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']
inFile = open('data/combined/'+run+'/postProcessedVFinal.csv','r')
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    bicluster = int(splitUp.pop(0))
    postProcessed[run][bicluster] = dict(zip(headers,splitUp))

inFile.close()
# Get the biclusters from the cMonkey object
all_biclusters = c1_p3utr.getBiclusters()
for bicluster in all_biclusters:
    # Initialize the bicluster
    b1 = Bicluster.objects.get(name=str(run)+'.'+str(bicluster))
    print bicluster
    biclusters[str(run)+'.'+str(bicluster)] = b1
    
    # Load p3utr motif
    # '3pUTRmotif1.weederScore','3pUTRmotif1.localPermP','3pUTRmotif1.allPermP','3pUTRmotif1.consensus','3pUTRmotif1.miRNAs','3pUTRmotif1.model'
    print 'Adding p3utr factor...'
    # Load p3utr motifs
    if not postProcessed[run][bicluster]['3pUTRmotif1.weederScore']=='NA':
        # Initialize object
        p1 = P3utr.objects.get(bicluster=b1)
        m_names = [m1.name for m1 in p1.mirnas.all()]
        # Add miRNA motif matches
        if not postProcessed[run][bicluster]['3pUTRmotif1.miRNAs']=='NA':
            print 'P3UTR Motif!'
            splitUp = postProcessed[run][bicluster]['3pUTRmotif1.miRNAs'].strip().split('_')
            if not 'hsa-miR-96-3p' in splitUp and 'hsa-mir-96-3p' in m_names:
                p1.mirnas.clear()
            for miRNA in splitUp:
                if miRNA[-3:]=='-5p':
                    miRNA = miRNA[:-3]
                miRs = miRNAInDict(miRNA.lower(), miRNAIDs)
                if len(miRs)>0:
                    for miRNA in miRs:
                        if not miRNAs[miRNA].name in m_names:
                            print 'Adding!'
                            p1.mirnas.add(miRNAs[miRNA])



###################
## Load GO terms ##
###################
geneOntology = {}
geneOntologyDag = {}
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
            geneOntology[id1] = Gene_Ontology.objects.get(go_id=id1, term=name, category=category)
            geneOntology[id1].save()
            goIds.append(id1)
            geneOntologyDag[id1] = []
    
    if line.strip()[0:6]=='alt_id' and category=='biological_process':
        geneOntology[line.strip().split('alt_id: ')[1]] = Gene_Ontology.objects.get(go_id=line.strip().split('alt_id: ')[1], term=name, category=category)
        goIds.append(line.strip().split('alt_id: ')[1])
        geneOntologyDag[line.strip().split('alt_id: ')[1]] = []
    
    if line.strip()[0:4]=='is_a' and category=='biological_process':
        parent = (line.strip().split('is_a: ')[1]).split(' ')[0]
        for goId in goIds:
            geneOntologyDag[goId].append(parent)

inFile.close()

##################################################
## Load GO terms and genes associated with them ##
##################################################
def getParents(go_id, geneOntologyDag):
    retMe = []
    if not len(geneOntologyDag[go_id])==0:
        for parent in geneOntologyDag[go_id]:
            retMe += getParents(parent, geneOntologyDag)
        return [go_id]+retMe
    else:
        return [go_id]

def uniquify(seq):
    # not order preserving
    set = {}
    map(set.__setitem__, seq, [])
    return set.keys()

inFile = open('data/gene2go.hsa','r')
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split('\t')
    if splitUp[2] in geneOntology:
        ga_all = Gene_Annotation.objects.filter(category='entrez',annotation=splitUp[1])
        if len(ga_all) > 0:
            for ga1 in ga_all:
                #geneOntology[splitUp[2]].annotated_genes.add(genes[ga1.gene.probe])
                parents = []
                for p1 in geneOntologyDag[splitUp[2]]:
                    parents += getParents(p1, geneOntologyDag)
                parents = uniquify(parents)
                for p1 in parents:
                    geneOntology[p1].annotated_genes.add(genes[ga1.gene.probe])

inFile.close()



########################################################
## Load up the cell-line differential expression data ##
########################################################
import re, cPickle, os
from glioma.models import *
from copy import deepcopy
from cMonkeyWrapper import cMonkeyWrapper
from bicluster import bicluster
from pssm import pssm

#######################################################################
## Create a dictionary to convert the miRNAs to there respective ids ##
#######################################################################
inFile = open('data/hsa.mature.fa','r')
miRNAIDs = {}
miRNAIDs_rev = {}
mature_sequence_ids = []
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.split(' ')
    mature_sequence_ids.append(splitUp[1])
    if not splitUp[1] in miRNAIDs_rev:
        miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()
    if not splitUp[0].lower() in miRNAIDs:
        miRNAIDs[splitUp[0].lower()] = splitUp[1]
    else:
        print 'Uh oh!',splitUp

def miRNAInDict(miRNA, dict1):
    retMe = []
    for i in dict1.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        if a[-3:]=='-3p':
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(a+'[a-oq-z]?(-\d)?$')
        if re1.match(b):
            return 1
    else:
        if b[-3:]=='-3p':
            re1 = re.compile(b+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(b+'[a-oq-z]?(-\d)?$')
        if re1.match(a):
            return 1
    return 0

# Load up the miRNAs
miRNAs = {}
m_all = MiRNA.objects.all()
for m1 in m_all:
    miRNAs[m1.mature_sequence_id] = deepcopy(m1)

# Read in file and store cell-line associations
# ,fc.T98G,t.p.T98G,bh.p.T98G,fc.A172,t.p.A172,bh.p.A172,fc.U87,t.p.U87,bh.p.U87,fc.AU,t.p.AU,bh.p.AU,fc.All,t.p.All,bh.p.All
inFile = open('data/output.csv','r')
inFile.readline()
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.strip().split(',')
    miRNA = splitUp.pop(0)
    if miRNA[-3:]=='-5p':
        miRNA = miRNA[:-3]
    miRs = miRNAInDict(miRNA.lower(), miRNAIDs)
    print miRNA, [miRNAs[m1].name for m1 in miRs]
    if len(miRs)>0:
        for miRNA in miRs:
            cl1 = Cell_Line(mirna=miRNAs[miRNA], t98g_fc = splitUp[0], t98g_p = splitUp[1], t98g_bh_p = splitUp[2], a172_fc = splitUp[3], a172_p = splitUp[4], a172_bh_p = splitUp[5], u87_fc = splitUp[6], u87_p = splitUp[7], u87_bh_p = splitUp[8], au_fc = splitUp[9], au_p = splitUp[10], au_bh_p = splitUp[11], all_fc = splitUp[12], all_p = splitUp[13], all_bh_p = splitUp[14])
            cl1.save()
    else:
        print miRNA

inFile.close()

