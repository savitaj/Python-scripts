#Python script ot create a cohort of patients who recieved 'imatinib' and/or have 'leukemia'
%load_ext autoreload
%autoreload 2
%matplotlib inline
sdkPath = '../'
import sys
sys.path.append(sdkPath)
from TestCohort import *

#connect to the sdkpath
rec = RecordsAPIWrapper(path=sdkPath, debug=True)

#Define a cohort using RecordAPIWrapper.makeCohort function. The first argument is the cohort name, second argument is the cohort Specifier
#selecting patients that receive imatinib
query = inQuery('med_drugs', 'imatinib')
new_cohort = rec.makeCohort('my_new_cohort', cohortSpecifier=query)

#choose specific fields and columns to be retrieved in a pandas dataframe
new_cohort.initDump(cohortProjector='patient_id,timestamp,order_drugs,disease') # Fields are separated by commas
df = new_cohort.dumpData()

# Instead of previous step, iterate over all chunks using getDF (prevents error with dumpData function and time taken)
# new_cohort.advanceDF() returns True if there is another df to advance to, and False otherwise
new_cohort.initDump(cohortProjector='patient_id,timestamp,order_drugs,disease')
while new_cohort.advanceDF():
    df = new_cohort.getDF() # Gets a new chunk each time called

#gets a list of diagnosis codes for a particular disease
codes = rec.getDiagnosisCodesFromDescription('leukemia')


#or get a list of diagnosis codes associated with 'diagnosis_code' field
query = inQuery('diagnosis_code', codes)
ilc = rec.makeCohort('imatinib_leukemia_cohort', cohortSpecifier=query)

#select patients who have imatinib in 'order_drugs' and 'med_durgs'  and patient has recieved a leukemia 'diagnosis_code'
query = andQuery(inQuery('order_drugs, med_drugs', 'imatinib'), inQuery('diagnosis_code', codes))
ilc = rec.makeCohort('imatinib_leukemia_cohort', cohortSpecifier=query)

#patients who recieved 'imatinib' or were assigned leukemia 'diagnosis_code'
query = orQuery(inQuery('order_drugs, med_drugs', 'imatinib'), inQuery('diagnosis_code', codes))
ilc = rec.makeCohort('imatinib_leukemia_cohort', cohortSpecifier=query)

#range query function builds a cohort of patients that have a value appearing in a selected range
query = andQuery(
    inQuery('diagnosis_code', codes),
    rangeQuery('Hemoglobin [Mass/volume] in Blood', rangeStart=50, rangeStop=100)
)

#create a cohort
sca = rec.makeCohort('imatinib_leukemia_cohort', cohortSpecifier=query)
