#!/bin/bash

export PATH="snap_path"/bin/snap/bin:$PATH
gptPath="gpt"

graphXmlPath="path_for_the_xml_file/BE-Bra_snap.xml"


removeExtension() {
    file="$1"
    echo "$(echo "$file" | sed -r 's/\.[^\.]*$//')"
}


sourceDirectory="sentinel2/L2A/BE-Bra"
targetDirectory="sentinel2/cube/BE-Bra"
targetFilePrefix="sub"

n_item=$(find * -maxdepth 0 | wc -l)
d=0
# A nice way to get notifications on Telegram
curl -s -X POST $URL -d chat_id=$CHAT_ID -d text="Starting to process BE-Bra (SNAP)"

for F in $(ls -1d "${sourceDirectory}"/S2*.SAFE); do
  d=$((d+1))
  sourceFile="$(realpath "$F")"
  targetFile="${targetDirectory}/${targetFilePrefix}_$(removeExtension "$(basename ${F})").nc"
  echo $sourceFile
  #echo $targetFile
  gpt ${graphXmlPath} -e -PtargetbasePath=$targetFile -PsourceFile=$sourceFile #-PgeoRegion=$polygon
  #curl -s -X POST $URL -d chat_id=$CHAT_ID -d text="SNAP processing for BE-Bra: $((100 * d / n_item)) %"
  rm -r /User/homes/dpabon/.snap/var/cache/s2tbx/l2a-reader/7.0.0/*
done


rename 's/sub_S2A_MSIL2A_//;' /Net/Groups/BGI/work_3/EFP_Trustee/sentinel2/cube/BE-Bra/*
rename 's/sub_S2B_MSIL2A_//;' /Net/Groups/BGI/work_3/EFP_Trustee/sentinel2/cube/BE-Bra/*

# A nice way to get notifications on Telegram
curl -s -X POST $URL -d chat_id=$CHAT_ID -d text="BE-Bra (SNAP) done"