#!/bin/bash

#=================================================
input_path="./TestImg_Align"  # specify input image folder here
# output_path="./TestImg_Align/"$(date +%Y_%m_%d_%H_%M_%S)
output_path="./TestImg_Align/Output"$(date +%Y_%m_%d_%H_%M_%S)
#=================================================

if [ ! -d "${output_path}" ]; then
    mkdir ${output_path}
fi

run_all(){
    echo "----------------------------------------------------------"
    echo "$(date +%Y.%m.%d) $(date +%T)"
    group_end=25
    max_exp=2
    max_con=2
    max_sat=2
    min_exp=1
    min_con=1
    min_sat=1
    step_exp=1
    step_con=1
    step_sat=1
    group_start=24
    weight_exposure=${min_exp}
    until [ $weight_exposure -ge ${max_exp} ]
    do
      weight_contrast=${min_con}
      until [ $weight_contrast -ge ${max_con} ]
      do
        weight_saturation=${min_sat}
        until [ $weight_saturation -ge ${max_sat} ]
        do
          i=${group_start}
          until [ $i -ge ${group_end} ]
          do
            _path_group=""
            for file in ${input_path}/${i}_*g; do
                # _filename="${file##*/}"
                # _filename_no_ext="${_filename%.*}"
                _path_group="${_path_group} ${file}"
                echo "\n+++ $file +++"
            done
            outname="output"${i}"_sat"${weight_saturation}"_exp"${weight_exposure}"_con"${weight_contrast}
            # echo "${_path_group}"
            ./hdrFusion $_path_group "${output_path}/${outname}.jpg" ${weight_saturation} ${weight_exposure} ${weight_contrast}
            i=$((i+1))
          done
          weight_saturation=$((weight_saturation+step_sat))
        done
        weight_contrast=$((weight_contrast+step_con))
      done
      weight_exposure=$((weight_exposure+step_exp))
    done
    echo "All Done. Results are saved to ${output_path}"
}

run_all
