#!/bin/sh


## Usage
# sh step1-fullCoverage.sh snyder
# sh step1-fullCoverage.sh hippo

# Define variables
EXPERIMENT=$1
SHORT="fullCov-${EXPERIMENT}"
CORES=10

# Directories
ROOTDIR=/dcs01/ajaffe/Brain/derRuns/derCountSupp
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/CoverageInfo

if [[ "${EXPERIMENT}" == "snyder" ]]
then
    DATADIR=/dcs01/ajaffe/Snyder/RNAseq/TopHat
    CUTOFF=5
elif [[ "${EXPERIMENT}" == "hippo" ]]
then
    DATADIR=/dcs01/ajaffe/Hippo/TopHat
    CUTOFF=3
else
    echo "Specify a valid experiment: snyder or hippo"
fi



# Construct shell file
echo 'Creating script for loading the Coverage data'
cat > ${ROOTDIR}/.${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=20G,h_vmem=45G,h_fsize=40G
#$ -N ${SHORT}
#$ -pe local ${CORES}

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load the data, save the coverage without filtering, then save each file separately
cd ${WDIR}
module load R/3.2.x
module load samtools/1.1
Rscript ${ROOTDIR}/step1-fullCoverage.R -d "${DATADIR}" -p "out$" -c "${CUTOFF}" -m ${CORES}

## Move log files into the logs directory
mv ${ROOTDIR}/${SHORT}.* ${WDIR}/logs/

echo '**** Job ends ****'
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
