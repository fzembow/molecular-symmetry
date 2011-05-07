
# check for correct arguments
if [ $# -neq 1 ]
then
echo "usage: ./run_hadoop.sh (input filenames file)"
exit 1
fi

# get reference to hadoop streaming jar
export SJAR=/usr/lib/hadoop/contrib/streaming/hadoop-streaming-0.20.2*.jar

# remove old output directory
rm -rf hadoop_output/

# spin off hadoop job
hadoop jar $SJAR -input $1 -output hadoop_output/ -mapper $(pwd)/mapper.py -numReduceTasks 0 -file mapper.py -file extract.py

