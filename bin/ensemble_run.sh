#!/bin/bash

printf "%%=================================================%%\n"
printf "%% An ensemble of aleatoric stochastic simulations %%\n"
printf "%%=================================================%%\n"

read -p "What's the ensemble size: " ensemble_size

for (( c=1; c<=$ensemble_size; c++ ))
do
    ./src/stoexsim/do_exsim_aleatoric.csh >/dev/null
done

echo "+++++ generating $ensemble_size simulations done +++++"
