#!/bin/bash

filename=$1
new_filename=${filename}_swapped

id_map_file="link_69328_19266_genotyped.txt"

declare -A id_map
echo "Reading id map file...$id_map_file"

while IFS=" " read -r key value; do
    id_map[$key]=$value
done < "$id_map_file"

echo "Done, performing sed ${#id_map[@]}"
value=$(<$filename)

for id in "${!id_map[@]}"
do
  value="${value/"$id"/"${id_map[$id]}"}"
done

echo "$value" > "$new_filename"

