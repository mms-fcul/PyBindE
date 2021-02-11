#!/bin/bash -e

dependencies=('pandas' 'numpy' 'scipy' 'freesasa')

for i in ${dependencies[*]};
do
  echo ${i}
  pip3 install ${i}
  pip3 install --upgrade ${i}
done