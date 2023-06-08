#!/bin/bash

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

printf "%%=================================================%%\n"
printf "%% An ensemble of aleatoric stochastic simulations %%\n"
printf "%%=================================================%%\n"

read -p "What's the ensemble size: " ensemble_size

for (( c=1; c<=$ensemble_size; c++ ))
do
    ./src/stoexsim/do_exsim_aleatoric.csh >/dev/null
done

echo "+++++ generating $ensemble_size simulations done +++++"
