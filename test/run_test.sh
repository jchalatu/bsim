OUTDIR=`pwd -P`
EXPORTPATH="${OUTDIR}/test_output"
FOLDER='benchmarking'
RUNDIR='/scratch/jchalatu/bingalls-research/bsim/out/production/bsim-bristol/'
cd $RUNDIR
JAVACMD='java -cp .:/scratch/jchalatu/bingalls-research/bsim/lib/* BasicSimulation2D.BasicSimulation2D'
ARGS='-pop 1 -simt 5.0 -simdt 0.05 -export_time 0.05'
DATAPATHS="-input_data /scratch/jchalatu/bingalls-research/bsim/run/initialization_data/onecell-1800by1800.csv -export_path ${EXPORTPATH}/${FOLDER}"

cmd="${JAVACMD} ${ARGS} ${DATAPATHS}"

echo $cmd

$cmd
