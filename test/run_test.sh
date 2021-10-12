OUTDIR=`pwd -P`

BSIM_DIR=$(echo "$OUTDIR" | sed -e 's|/test||')
EXPORTPATH="${OUTDIR}/test_output"
FOLDER='benchmarking'
RUNDIR="${BSIM_DIR}/out/production/bsim-bristol/"
cd $RUNDIR
JAVACMD="java -cp .:${BSIM_DIR}/lib/* BasicSimulation2D.BasicSimulation2D"
ARGS='-simt 1.0 -simdt 0.05 -export_time 0.05'
DATAPATHS="-input_data ${BSIM_DIR}/run/initialization_data/onecell-1800by1800.csv -export_path ${EXPORTPATH}/${FOLDER}"

cmd="${JAVACMD} ${ARGS} ${DATAPATHS}"

echo $cmd

$cmd
