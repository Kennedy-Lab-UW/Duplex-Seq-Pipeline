#!/bin/bash

# Reset the test case
rm testConfig.csv.summary*
cd testData
rm -r Final/ Intermediate/ Stats/ logs/ testData_config.sh
cd ../
echo "done"
