
# remove existing output, if any
hadoop fs -rmr hadoop_output/
rm -rf hadoop_output/

# get reference to hadoop streaming jar
export SJAR=/usr/lib/hadoop/contrib/streaming/hadoop-streaming-0.20.2*.jar

# remove old output directory
rm -rf hadoop_output/

# spin off hadoop job
hadoop jar $SJAR -input urls-full.txt -output hadoop_output/ -mapper $(pwd)/mapper.py -reducer $(pwd)/reducer.py -numReduceTasks 500 -file mapper.py -file extract.py -file reducer.py

# retrieve results
hadoop fs -copyToLocal hadoop_output/ hadoop_output/

# concatenate files
cat hadoop_output/part* ./results.txt

# change commas to tabs in the results
perl -pi -e 's/,/\t/g' results.txt 

# sort the results
sort -n -o results-sorted.txt -k 1 results.txt
rm results.txt
mv results-sorted.txt results.txt