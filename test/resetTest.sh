#!/bin/bash

# Reset the test case
rm testConfig.csv.summary*
cd testData
rm -r Final/ Intermediate/ Stats/ logs/
cd ../
echo "done"
