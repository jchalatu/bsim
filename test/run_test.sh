OUTDIR=`pwd -P`
BSIMDIR=$(echo "$OUTDIR" | sed -e 's|/test||')
EXPORTPATH="${OUTDIR}/test_output"
FOLDER='benchmarking'
RUNDIR="${BSIMDIR}/out/production/bsim-bristol/"
cd $RUNDIR
JAVACMD="java -cp .:${BSIMDIR}/lib/* BasicSimulation2D.BasicSimulation2D"
ARGS='-simt 1.0 -simdt 0.05 -export_time 0.05'
DATAPATHS="-input_data ${BSIMDIR}/run/initialization_data/onecell-1800by1800.csv -export_path ${EXPORTPATH}/${FOLDER}"

cmd="${JAVACMD} ${ARGS} ${DATAPATHS}"

echo $cmd

$cmd
