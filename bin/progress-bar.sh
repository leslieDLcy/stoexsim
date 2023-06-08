#!/usr/bin/env bash

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

progress-bar() {
  local duration=${1}


    already_done() { for ((done=0; done<$elapsed; done++)); do printf "▇"; done }
    remaining() { for ((remain=$elapsed; remain<$duration; remain++)); do printf " "; done }
    percentage() { printf "| %s%%" $(( (($elapsed)*100)/($duration)*100/100 )); }
    clean_line() { printf "\r"; }

  for (( elapsed=1; elapsed<=$duration; elapsed++ )); do
      already_done; remaining; percentage
      sleep 1
      clean_line
  done
  clean_line
}



# SLEEP_DURATION=${SLEEP_DURATION:=1}  # default to 1 second, use to speed up tests

# progress-bar() {
#   local duration
#   local columns
#   local space_available
#   local fit_to_screen  
#   local space_reserved

#   space_reserved=6   # reserved width for the percentage value
#   duration=${1}
#   columns=$(tput cols)
#   space_available=$(( columns-space_reserved ))

#   if (( duration < space_available )); then 
#   	fit_to_screen=1; 
#   else 
#     fit_to_screen=$(( duration / space_available )); 
#     fit_to_screen=$((fit_to_screen+1)); 
#   fi

#   already_done() { for ((done=0; done<(elapsed / fit_to_screen) ; done=done+1 )); do printf "▇"; done }
#   remaining() { for (( remain=(elapsed/fit_to_screen) ; remain<(duration/fit_to_screen) ; remain=remain+1 )); do printf " "; done }
#   percentage() { printf "| %s%%" $(( ((elapsed)*100)/(duration)*100/100 )); }
#   clean_line() { printf "\r"; }

#   for (( elapsed=1; elapsed<=duration; elapsed=elapsed+1 )); do
#       already_done; remaining; percentage
#       sleep "$SLEEP_DURATION"
#       clean_line
#   done
#   clean_line
# }