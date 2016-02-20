## Usage
# sh step7-regMatVsDERs.sh snyder run3-v1.0.10
# sh step7-regMatVsDERs.sh hippo run3-v1.0.10

# Define variables
EXPERIMENT=$1
PREFIX=$2
SHORT="regVsDERs-${EXPERIMENT}"
ncore=5
cores="${ncore}cores"

# Directories
ROOTDIR=/dcs01/ajaffe/Brain/derRuns/derCountSupp
MAINDIR=${ROOTDIR}/${EXPERIMENT}

# Construct shell files
sname="${SHORT}.${PREFIX}"
echo "Creating script ${sname}"

if [[ "${EXPERIMENT}" == "snyder" ]]
then
    CUTOFF=5
elif [[ "${EXPERIMENT}" == "hippo" ]]
then
    CUTOFF=3
else
    echo "Specify a valid experiment: snyder or hippo"
fi

WDIR=${MAINDIR}/regionMatrix-vs-DERs/cut${CUTOFF}-vs-${PREFIX}

cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=200G,h_vmem=400G,h_fsize=30G
#$ -N ${sname}
#$ -hold_jid regMat-${EXPERIMENT},derM-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Compare DERs vs regionMatrix
cd ${WDIR}
module load R/3.2.x
Rscript -e "analysisPath <- '${WDIR}'; load('${MAINDIR}/regionMatrix/regionMat-cut${CUTOFF}.Rdata'); proc.time(); load('${MAINDIR}/derAnalysis/${PREFIX}/fullRegions.Rdata'); proc.time(); library(rmarkdown); library(knitrBootstrap); render('${ROOTDIR}/step7-regMatVsDERs.Rmd', output_file='${WDIR}/step7-regMatVsDERs.html')"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
