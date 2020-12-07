#!/usr/bin/env bash

COMPLETED_FILES_LIST="completed-files"

for url in $(ag --no-numbers -v '^#' dropbox-url-list) ; do

    filename=${url:69}

    if [[ -f "${COMPLETED_FILES_LIST}" && -n $( ag --no-numbers "${filename}" "${COMPLETED_FILES_LIST}" | ag -v '^#' ) ]] ; then
        echo "=> already processed file ${filename}"
        continue
    fi

    if [[ ! -f ${filename} ]] ; then
        echo "===== downloading ${filename} from ${url} ====="
        curl -OL $url
    fi

    files_name=${filename/res/files}
    source_dirname=${files_name/.zip/}
    if [[ ! -d ${source_dirname} ]] ; then
        mkdir ${source_dirname}
        pushd ${source_dirname}

        unzip -o ../${filename}

        if [[ -d ${source_dirname} ]] ; then
            echo "found embedded ${source_dirname} directory, removing"
            mv ${source_dirname}/* .
            rmdir ${source_dirname}
        fi

        popd
    fi

    if [[ ! -d ${source_dirname} ]] ; then
        echo "-> couldn't find dir ${source_dirname}, something's wrong"
        exit 1
    fi

    scen=${source_dirname#*_scen}
    echo "----- processing scenario ${scen} at $(date) -----"

    for iu in $(awk -F , '{print $2}'< $source_dirname/scen${scen}IUscountry.csv | sed 's/"//g' | grep -vw x) ; do

        country=${iu:0:3}
        found_code=${iu:3}
        padded_code=$( printf %05.0f ${found_code} )
        padded_iu=${country}${padded_code}

        iu_data_dir=data/${country}/${padded_iu}
        iu_files=$(ls ${source_dirname}/*${iu}*.csv ${source_dirname}/*${iu}*.p 2>/dev/null | tr '\n' ' ')

        for file in ${iu_files} ; do
            if [[ ! -d ${iu_data_dir} ]] ; then
                echo "-> making data directory ${iu_data_dir}"
                mkdir -p ${iu_data_dir}
            fi
            renamed_file=$( echo $file | sed -e "s/${found_code}/${padded_code}/g" -e "s/${source_dirname}/data\/${country}\/${padded_iu}/g" )
            cmd="mv $file $renamed_file"
            echo ${cmd}
            eval "${cmd}"
        done

    done

    echo "-> removing file ${filename}"
    rm -f ${filename}

    echo "-> removing dir ${dirname}"
    rm -rf "${dirname}"

    echo $filename >> ${COMPLETED_FILES_LIST}

    echo "----- rsync-ing scenario ${scen} -----"
    gsutil -m rsync -r ./data gs://ntd-disease-simulator-data/diseases/sch-mansoni/source-data

    echo "----- rm-ing data folders -----"
    rm -rf data
    rm -rf ${source_dirname}

    echo
    echo

done
